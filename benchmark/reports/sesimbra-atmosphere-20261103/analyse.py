"""Score the recorded optical comparison; never read cover targets."""

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio import windows

from kelpobserve import config
from kelpobserve.local_cover import clean, digest
from kelpobserve.optical_corrections import GlintCalibrationError, fit_block_glint
from kelpobserve.pneo_optical_review import ROLES, NativeCorrection, block_ids, floor
from kelpobserve.water_constraints import WaterConstraintError, separate_iops

CACHE = config.load().require("cache") / "optical-atmosphere-20261103"
AUDIT = config.BENCHMARK / "reports/r7-optical-diagnostic-20261102"
S2 = "s2-20260702-S2B-R037"
PN = "pneo-20260708"
WAVES = {S2: [442, 492, 559, 665, 833], PN: [439, 486, 563, 655, 826]}


def dump(p, x):
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + "\n")


def read_candidate(native, directory, prefix, scene):
    bbox = rasterio.transform.array_bounds(*native.win.shape, native.win.transform)
    data = []
    used = []
    for wave in WAVES[scene]:
        p = next(directory.glob(f"*_{prefix}_{wave}.tif"))
        used.append(p)
        with rasterio.open(p) as src:
            w = windows.from_bounds(*bbox, transform=src.transform).round_offsets().round_lengths()
            assert np.allclose(tuple(src.window_transform(w)), tuple(native.win.transform), atol=1e-6, rtol=0)
            assert src.crs == native.win.crs
            data.append(native.rv.read_clipped(p, w))
    native.nir = data[-1]
    native.data = dict(zip(ROLES, data[:-1]))
    good = np.all(np.isfinite(data), axis=0)
    fraction = native.rv._average(
        np.where(native.valid, good.astype(float), np.nan),
        native.win,
        native.transform,
        native.shape,
        native.win.crs,
        np.nan,
    )[native.row, native.col]
    y = native.apply(np.zeros(4), 0)
    x = native.rv._average(
        np.where(native.valid, native.nir, np.nan), native.win, native.transform, native.shape, native.win.crs, np.nan
    )[native.row, native.col]
    complete = (fraction >= 1 - 1e-7) & np.isfinite(y).all(axis=1) & np.isfinite(x)
    y[~complete] = np.nan
    x[~complete] = np.nan
    return x, y, complete, used


def stats(frame):
    records = []
    for (scene, branch, fold, region), g in frame.groupby(["scene", "branch", "fold", "region"]):
        for band in ROLES:
            y = g["rhos_" + band].to_numpy()
            finite = np.isfinite(y)
            residuals = []
            for _, block in g[finite].groupby("block"):
                if len(block) >= 10:
                    residuals.extend(block["rhos_" + band] - block["rhos_" + band].median())
            records.append(
                {
                    "scene": scene,
                    "branch": branch,
                    "fold": fold,
                    "region": region,
                    "band": band,
                    "candidate_cells": len(g),
                    "valid_cells": int(finite.sum()),
                    "coverage": float(finite.mean()),
                    "negative_fraction": float((y[finite] < 0).mean()) if finite.any() else None,
                    "below_minus_0p001": float((y[finite] < -0.001).mean()) if finite.any() else None,
                    "median": float(np.median(y[finite])) if finite.any() else None,
                    "mad_sigma": float(1.4826 * np.median(np.abs(residuals))) if residuals else None,
                    "within_nir_support": float(g.within_nir_support.mean())
                    if g.within_nir_support.notna().any()
                    else None,
                }
            )
    return pd.DataFrame(records)


def run(output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    source = AUDIT / "paired-grid-cells.csv.gz"
    assert digest(source) == json.loads((AUDIT / "outputs.json").read_text())["sha256"][source.name]
    cols = ["scene", "region", "x", "y", "depth_m", "shore_m"] + [
        f"rhos_{b}_{v}" for b in (*ROLES, "nir") for v in ["V0", "V1"]
    ]
    g = pd.read_csv(source, usecols=cols)
    g = g[g.scene.isin([S2, PN])].copy()
    metas = json.loads((AUDIT / "scene-metadata.json").read_text())
    inputs = {source}
    all_rows = []
    calibrations = {}
    reproduction = {}
    support = {}
    for scene, part in g.groupby("scene"):
        print("Scoring", scene, flush=True)
        native = NativeCorrection(scene, part)
        inputs.update(native.inputs)
        block = block_ids(part)
        deep = (part.region == "deep_reference").to_numpy()
        bx = np.floor(part.x.to_numpy() / 100) * 100 + 50
        split = np.median(np.unique(bx[deep]))
        west = deep & (bx <= split - 50)
        east = deep & (bx >= split + 50)
        folds = [
            ("west_to_east", west, east),
            ("east_to_west", east, west),
            ("full_reference_preview", deep, np.ones(len(part), bool)),
        ]

        def append(branch, y, fold, test, support=None, part=part, block=block):
            row = part.loc[test, ["scene", "region", "x", "y", "depth_m", "shore_m"]].copy()
            row["block"] = block[test]
            row["branch"] = branch
            row["fold"] = fold
            row["within_nir_support"] = support[test] if support is not None else np.nan
            for i, b in enumerate(ROLES):
                row["rhos_" + b] = y[test, i]
            all_rows.append(row)

        # Earlier branches remain explicitly separate controls.
        for v in ["V0", "V1"]:
            y = part[[f"rhos_{b}_{v}" for b in ROLES]].to_numpy()
            for fold, _, test in folds:
                append("archived_" + v, y, fold, test)
        candidates = []
        if scene == S2:
            for branch in ["s2_radcor_fixed", "s2_radcor_auto"]:
                for glint in [False, True]:
                    for prefix in ["rhosu", "rhos"]:
                        candidates.append(
                            (
                                branch + ("_SWIR" if glint else "_no_glint") + "_" + prefix,
                                CACHE / branch / ("glint" if glint else ""),
                                prefix,
                            )
                        )
        else:
            candidates = [(b, CACHE / b, "rhos") for b in ["pneo_marine_control", "pneo_continental"]]
        for branch, directory, prefix in candidates:
            x, y, complete, files = read_candidate(native, directory, prefix, scene)
            inputs.update(files)
            support[branch] = {"cells": len(part), "complete_native_support_cells": int(complete.sum())}
            if branch == "pneo_marine_control":
                old = part[[f"rhos_{b}_V0" for b in ROLES]].to_numpy()
                reproduction[scene] = float(np.nanmax(np.abs(old - y)))
                assert reproduction[scene] < 2e-7, reproduction
            for fold, train, test in folds:
                append(branch, y, fold, test)
                if scene != PN:
                    continue
                fitmask = train & complete
                try:
                    cal = fit_block_glint(x[fitmask], y[fitmask], block[fitmask])
                except GlintCalibrationError as error:
                    calibrations[branch + "_" + fold] = {"status": "unavailable", "reason": str(error)}
                    continue
                inrange = (x >= cal.nir_p01) & (x <= cal.nir_p99) & complete
                calibrations[branch + "_" + fold] = {
                    "slopes": cal.slopes,
                    "nir_p05": cal.nir_baseline,
                    "nir_p01": cal.nir_p01,
                    "nir_p99": cal.nir_p99,
                    "informative_blocks": cal.block_count,
                    "training_cells": int(fitmask.sum()),
                    "test_cells": int(test.sum()),
                    "disjoint": not bool(set(block[train]) & set(block[test]))
                    if fold != "full_reference_preview"
                    else None,
                }
                for kind, base in [("empirical_zero", 0), ("empirical_p05", cal.nir_baseline)]:
                    corrected = native.apply(cal.slopes, base)
                    corrected[~complete] = np.nan
                    append(branch + "_" + kind, corrected, fold, test, inrange)
    frame = pd.concat(all_rows, ignore_index=True)
    frame["all_bands_valid"] = np.isfinite(frame[["rhos_" + b for b in ROLES]]).all(axis=1)
    frame["common_support"] = frame.groupby(["scene", "fold", "x", "y"]).all_bands_valid.transform("all")
    frame.to_csv(out / "predictions.csv.gz", index=False)
    evaluation = stats(frame)
    evaluation.to_csv(out / "evaluation.csv", index=False)
    stats(frame[frame.common_support]).to_csv(out / "evaluation-common.csv", index=False)
    changes = []
    pairs = []
    for state in ["fixed", "auto"]:
        for glint in ["no_glint", "SWIR"]:
            pairs.append((f"s2_radcor_{state}_{glint}_rhosu", f"s2_radcor_{state}_{glint}_rhos", "adjacency"))
        pairs.append((f"s2_radcor_{state}_no_glint_rhos", f"s2_radcor_{state}_SWIR_rhos", "SWIR_glint"))
    pairs.append(("s2_radcor_fixed_SWIR_rhos", "s2_radcor_auto_SWIR_rhos", "aerosol_solution"))
    full = frame[(frame.fold == "full_reference_preview") & frame.common_support]
    for before, after, effect in pairs:
        keys = ["scene", "region", "x", "y"]
        a = full[full.branch == before].set_index(keys)
        b = full[full.branch == after].set_index(keys)
        assert a.index.equals(b.index)
        for region in ["deep_reference", "habitat"]:
            mask = a.index.get_level_values("region") == region
            for band in ROLES:
                delta = b.loc[mask, "rhos_" + band] - a.loc[mask, "rhos_" + band]
                changes.append(
                    {
                        "effect": effect,
                        "before": before,
                        "after": after,
                        "region": region,
                        "band": band,
                        "n": len(delta),
                        "median_change": delta.median(),
                        "p05": delta.quantile(0.05),
                        "p95": delta.quantile(0.95),
                    }
                )
    pd.DataFrame(changes).to_csv(out / "paired-changes.csv", index=False)
    # Depth and shore-distance strata are descriptive, with fixed bins.
    preview = frame[frame.fold == "full_reference_preview"].copy()
    preview["shore_bin"] = pd.cut(
        preview.shore_m, [-np.inf, 100, 300, 1000, np.inf], labels=["<=100m", "100-300m", "300-1000m", ">1000m"]
    )
    strata = preview.groupby(["scene", "branch", "region", "shore_bin"], observed=True)[
        ["rhos_blue", "rhos_green"]
    ].agg(["count", "median"])
    strata.to_csv(out / "shore-strata.csv")
    spectra_path = config.OPTICS / "reports/e5s-kelp-spectra/results.json"
    products_path = config.BENCHMARK / "registration/evidence/cmems-optical-inputs-20261102.json"
    spectra = json.loads(spectra_path.read_text())
    products = {t["date"]: t for t in json.loads(products_path.read_text())["targets"]}
    states = {}
    closures = []
    inputs.update([spectra_path, products_path, AUDIT / "scene-metadata.json", AUDIT / "water-states.json"])
    for scene in [S2, PN]:
        meta = metas[scene]
        p = products[meta["date"]]["variables"]
        states[scene] = []
        for shape in ["median", "2025-12-24", "2026-02-22"]:
            e = copy.deepcopy(spectra)
            e["iops"]["median"] = e["iops"][shape]
            for ks in [-1, 0, 1]:
                for bs in [-1, 0, 1]:
                    kd = p["KD490"]["value"] * (1 + ks * p["KD490"]["uncertainty_pct"] / 100)
                    bb = p["BBP"]["value"] * (1 + bs * p["BBP"]["uncertainty_pct"] / 100)
                    record = {"state": f"{shape}_k{ks:+d}_b{bs:+d}", "shape": shape}
                    try:
                        record.update(
                            separate_iops(
                                kd,
                                bb,
                                meta["sza"],
                                [meta["wavelengths"][b] for b in ["blue", "green"]],
                                e,
                                cdm443=p["CDM"]["value"],
                                cdm_error_pct=p["CDM"]["uncertainty_pct"],
                            )
                        )
                    except WaterConstraintError as error:
                        record.update(status="infeasible", reason=str(error))
                    states[scene].append(record)
        old = json.loads((AUDIT / "water-states.json").read_text())[scene][0]
        states[scene].append(dict(old, state="r7_nominal", shape="r7"))
        for (branch, fold), part in frame[(frame.scene == scene) & (frame.region == "deep_reference")].groupby(
            ["branch", "fold"]
        ):
            obs = part[["rhos_blue", "rhos_green"]].to_numpy()
            ok = np.isfinite(obs).all(axis=1)
            z = part.depth_m.to_numpy() + meta["water_level_m"]
            for state in states[scene]:
                if "a" not in state:
                    continue
                lo = floor(state["a"], state["bb"], meta["sza"], z)
                hi = floor(state["a"], state["bb"], meta["sza"], z, 0.6)
                for tol in [0.001, 0.002]:
                    compatible = ok & ((obs >= lo - tol) & (obs <= hi + tol)).all(axis=1)
                    closures.append(
                        {
                            "scene": scene,
                            "branch": branch,
                            "fold": fold,
                            "state": state["state"],
                            "shape": state["shape"],
                            "cdm_conflict": state.get("cdm_conflict"),
                            "tolerance": tol,
                            "candidate_cells": len(obs),
                            "valid_cells": int(ok.sum()),
                            "compatible_fraction_all": float(compatible.mean()),
                            "compatible_fraction_valid": float(compatible[ok].mean()) if ok.any() else None,
                            "blue_residual": float(np.median((obs - lo)[ok, 0])) if ok.any() else None,
                            "green_residual": float(np.median((obs - lo)[ok, 1])) if ok.any() else None,
                        }
                    )
    pd.DataFrame(closures).to_csv(out / "water-closure.csv", index=False)
    # Common support recomputation retains the water grid and both tolerances.
    common_closures = []
    for (scene, branch, fold), part in frame[(frame.region == "deep_reference") & frame.common_support].groupby(
        ["scene", "branch", "fold"]
    ):
        meta = metas[scene]
        obs = part[["rhos_blue", "rhos_green"]].to_numpy()
        z = part.depth_m.to_numpy() + meta["water_level_m"]
        for state in states[scene]:
            if "a" not in state:
                continue
            lo = floor(state["a"], state["bb"], meta["sza"], z)
            hi = floor(state["a"], state["bb"], meta["sza"], z, 0.6)
            for tol in [0.001, 0.002]:
                match = ((obs >= lo - tol) & (obs <= hi + tol)).all(axis=1)
                common_closures.append(
                    {
                        "scene": scene,
                        "branch": branch,
                        "fold": fold,
                        "state": state["state"],
                        "shape": state["shape"],
                        "cdm_conflict": state.get("cdm_conflict"),
                        "tolerance": tol,
                        "cells": len(obs),
                        "compatible_fraction": float(match.mean()),
                        "blue_residual": float(np.median((obs - lo)[:, 0])),
                        "green_residual": float(np.median((obs - lo)[:, 1])),
                    }
                )
    pd.DataFrame(common_closures).to_csv(out / "water-closure-common.csv", index=False)
    dump(out / "water-states.json", states)
    dump(out / "calibration.json", calibrations)
    dump(
        out / "verification.json",
        {
            "pneo_control_max_absolute_rhos_error": reproduction,
            "support_by_branch": support,
            "raster_grid_matches_original": True,
            "cover_targets_read": False,
            "cover_scores_computed": False,
            "aerosol_estimation_not_spatially_held_out": True,
        },
    )
    dump(
        out / "provenance.json",
        {
            "inputs_sha256": {str(p): digest(p) for p in sorted(inputs)},
            "analysis_sha256": digest(Path(__file__)),
            "processing": "4 bounded ACOLITE branches; no source settings changed after scoring",
        },
    )
    print("Saved", out, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--cache-root", type=Path, default=CACHE)
    a = p.parse_args()
    CACHE = a.cache_root
    run(a.output)
