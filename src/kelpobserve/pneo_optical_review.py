"""Run the recorded PNeo glint and water review without reading habitat labels."""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio import windows
from rasterio.transform import from_origin

from kelpobserve import config
from kelpobserve.coastal.inversion import lee
from kelpobserve.coastal.optics import water
from kelpobserve.local_cover import clean, digest, legacy_readers
from kelpobserve.optical_corrections import GlintCalibrationError, fit_block_glint
from kelpobserve.pneo_optical_report import render
from kelpobserve.water_constraints import WaterConstraintError, separate_iops

ROLES = ("coastal", "blue", "green", "red")
AUDIT = config.BENCHMARK / "reports/r7-optical-diagnostic-20261102"
DESIGN = config.BENCHMARK / "registration/pneo-optical-repair-20261103.json"


class NativeCorrection:
    """Apply the nonlinearity before averaging, on the unchanged r1 water mask."""

    def __init__(self, scene, part):
        self.inputs = []
        x0, x1 = part.x.min() - 5, part.x.max() + 5
        y0, y1 = part.y.min() - 5, part.y.max() + 5
        self.transform = from_origin(x0, y1, 10, 10)
        self.shape = (round((y1 - y0) / 10), round((x1 - x0) / 10))
        self.row = np.rint((y1 - 5 - part.y.to_numpy()) / 10).astype(int)
        self.col = np.rint((part.x.to_numpy() - x0 - 5) / 10).astype(int)
        with legacy_readers(config.load()) as (rb, rv, spec):
            self.rv = rv
            self.win = rv.read_window(spec, scene, (x0 - 20, y0 - 20, x1 + 20, y1 + 20), "V0")
            self.valid = self.win.valid()
            self.nir = self.win.bands["nir"]
            self.data = {b: self.win.bands[b] for b in ROLES if b != "coastal"}
            rasters = rb.acolite_rasters(rv.scene_dir(spec, scene, "V0"), "rhos")
            path = rb.nearest_band(rasters, 443)[1]
            with rasterio.open(path) as src:
                bbox = rasterio.transform.array_bounds(*self.win.shape, self.win.transform)
                w = windows.from_bounds(*bbox, transform=src.transform).round_offsets().round_lengths()
                self.data["coastal"] = rv.read_clipped(path, w)
            for role in (*ROLES, "nir"):
                target = {"coastal": 443, "blue": 490, "green": 560, "red": 655, "nir": 825}[role]
                self.inputs.append(rb.nearest_band(rasters, target)[1])
            self.inputs.append(rv.mask_band(spec, rv.scene_dir(spec, scene, "V0"))[0])
            self.inputs.extend([Path(rb.__file__), Path(rv.__file__)])
        self.index = part.index

    def apply(self, slopes, baseline):
        result = []
        amplitude = np.maximum(self.nir - baseline, 0)
        # Preserve the pinned ACOLITE correction mask.
        amplitude = np.where(self.nir < 0.1, amplitude, 0)
        for band, slope in zip(ROLES, slopes):
            corrected = np.where(self.valid, self.data[band] - amplitude * slope, np.nan)
            grid = self.rv._average(corrected, self.win, self.transform, self.shape, self.win.crs, np.nan)
            result.append(grid[self.row, self.col])
        return np.column_stack(result)


def write(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def block_ids(g):
    return (np.floor(g.x / 100).astype(int).astype(str) + "_" + np.floor(g.y / 100).astype(int).astype(str)).to_numpy()


def within_block_slopes(x, y, ids):
    slopes = []
    for block in np.unique(ids):
        m = ids == block
        if m.sum() < 10 or np.quantile(x[m], 0.95) - np.quantile(x[m], 0.05) < 0.002:
            continue
        dx = x[m] - np.mean(x[m])
        dy = y[m] - np.mean(y[m], axis=0)
        slopes.append((dx[:, None] * dy).sum(axis=0) / (dx * dx).sum())
    return np.median(slopes, axis=0) if slopes else np.full(y.shape[1], np.nan)


def floor(a, bb, sza, z, rho=0):
    a, bb, z = np.asarray(a), np.asarray(bb), np.asarray(z)
    u = bb / (a + bb)
    ic = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(sza)))
    kc = (ic + lee._du_column(u)) * (a + bb)
    kb = (ic + lee._du_bottom(u)) * (a + bb)
    sub = (0.089 * u + 0.125 * u * u) * (1 - np.exp(-z[:, None] * kc)) + rho / np.pi * np.exp(-z[:, None] * kb)
    return np.pi * water.rrs_subsurface_to_above(sub)


def metrics(g, corrected, support, scene, fold, method, training_blocks):
    rows = []
    x = g.rhos_nir_V0.to_numpy()
    ids = block_ids(g)
    slopes = within_block_slopes(x, corrected, ids)
    for i, role in enumerate(ROLES):
        y = corrected[:, i]
        pieces = []
        for block in np.unique(ids):
            mask = ids == block
            if mask.sum() >= 10:
                pieces.append(y[mask] - np.median(y[mask]))
        residual = np.concatenate(pieces) if pieces else np.array([np.nan])
        rows.append(
            {
                "scene": scene,
                "fold": fold,
                "method": method,
                "band": role,
                "cells": len(g),
                "blocks": len(np.unique(ids)),
                "training_blocks": training_blocks,
                "within_support_fraction": np.mean(support),
                "negative_fraction": np.mean(y < 0),
                "below_minus_0p001_fraction": np.mean(y < -0.001),
                "p05": np.quantile(y, 0.05),
                "median": np.median(y),
                "p95": np.quantile(y, 0.95),
                "residual_nir_slope": slopes[i],
                "within_block_mad_sigma": 1.4826 * np.nanmedian(np.abs(residual)),
                "negative_fraction_within_support": np.mean(y[support] < 0) if support.any() else None,
            }
        )
    return rows


def run(output, design_path=DESIGN):
    historical = config.BENCHMARK / "reports/r7-local-cover-20261102-r3"
    old_files = {p: digest(p) for folder in (historical, AUDIT) for p in folder.rglob("*") if p.is_file()}
    frozen = json.loads((historical / "freeze.json").read_text())["code_sha256"]
    if any(digest(config.REPO / path) != value for path, value in frozen.items()):
        raise RuntimeError("Frozen r7 source has changed; review provenance before proceeding")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    design = json.loads(Path(design_path).read_text())
    write(output / "design.json", design)
    source = AUDIT / "paired-grid-cells.csv.gz"
    provenance = json.loads((AUDIT / "outputs.json").read_text())["sha256"]
    assert digest(source) == provenance[source.name]
    columns = ["scene", "region", "x", "y", "depth_m", "shore_m"]
    columns += [f"rhos_{b}_{v}" for b in (*ROLES, "nir") for v in ("V0", "V1")]
    g = pd.read_csv(source, usecols=columns)
    g = g[g.scene.isin(design["scenes"])].copy()
    scene_meta = json.loads((AUDIT / "scene-metadata.json").read_text())
    old_states = json.loads((AUDIT / "water-states.json").read_text())
    spectra_path = config.OPTICS / "reports/e5s-kelp-spectra/results.json"
    audit_path = config.BENCHMARK / "registration/evidence/cmems-optical-inputs-20261102.json"
    spectra = json.loads(spectra_path.read_text())
    products = {t["date"]: t for t in json.loads(audit_path.read_text())["targets"]}
    records, calibrations, assignments, closure, water_records, previews = [], {}, [], [], {}, []
    native_inputs = set()
    reconstruction_errors = {}
    for scene, part in g.groupby("scene"):
        print("Reviewing", scene, flush=True)
        deep = part[part.region == "deep_reference"].copy()
        bx = np.floor(deep.x / 100) * 100 + 50
        split = float(np.median(np.unique(bx)))
        west = deep[bx <= split - 50]
        east = deep[bx >= split + 50]
        folds = [("west_to_east", west, east), ("east_to_west", east, west)]
        theoretical = np.array(
            [np.median((deep[f"rhos_{b}_V0"] - deep[f"rhos_{b}_V1"]) / deep.rhos_nir_V0) for b in ROLES]
        )
        # The coefficient is constant in the pinned fixed-aerosol ACOLITE branch.
        native = NativeCorrection(scene, part)
        native_inputs.update(native.inputs)
        reconstructed = native.apply(theoretical, 0)
        observed = part[[f"rhos_{b}_V1" for b in ROLES]].to_numpy()
        reconstruction_errors[scene] = float(np.max(np.abs(reconstructed - observed)))
        if reconstruction_errors[scene] >= 2e-7:
            raise RuntimeError(f"Native V1 reproduction failed: {reconstruction_errors[scene]}")
        print("Native V1 reproduction error", reconstruction_errors[scene], flush=True)
        meta = scene_meta[scene]
        product = products[meta["date"]]["variables"]
        states = []
        for ks in (-1, 0, 1):
            for bs in (-1, 0, 1):
                kd = product["KD490"]["value"] * (1 + ks * product["KD490"]["uncertainty_pct"] / 100)
                bb = product["BBP"]["value"] * (1 + bs * product["BBP"]["uncertainty_pct"] / 100)
                name = f"KD_BBP_k{ks:+d}_b{bs:+d}"
                try:
                    state = separate_iops(
                        kd,
                        bb,
                        meta["sza"],
                        [meta["wavelengths"][b] for b in ("blue", "green")],
                        spectra,
                        cdm443=product["CDM"]["value"],
                        cdm_error_pct=product["CDM"]["uncertainty_pct"],
                    )
                    states.append({"state": name, **state})
                except WaterConstraintError as error:
                    states.append({"state": name, "status": "infeasible", "reason": str(error)})
        states.insert(0, {"state": "r7_nominal", **old_states[scene][0]})
        water_records[scene] = states
        calibrations[scene] = {"theoretical_slopes": theoretical, "split_easting_m": split, "folds": {}}
        # Full-reference fitting is a development preview, never a held-out score.
        folds.append(("full_reference_preview", deep, part))
        for fold, train, test in folds:
            train = train.copy()
            test = test.copy()
            ids = block_ids(train)
            if fold != "full_reference_preview":
                assert not set(ids) & set(block_ids(test))
            if len(train) < 1 or len(test) < 1:
                raise RuntimeError("Empty spatial split")
            baseline = float(np.quantile(train.rhos_nir_V0, 0.05))
            nir_range = np.quantile(train.rhos_nir_V0, [0.01, 0.99])
            info = {
                "training_cells": len(train),
                "test_cells": len(test),
                "training_blocks": len(np.unique(ids)),
                "test_blocks": len(np.unique(block_ids(test))),
                "nir_baseline_p05": baseline,
                "nir_support": nir_range,
            }
            for side, data in [("training", train), ("evaluation", test)]:
                for block in np.unique(block_ids(data)):
                    assignments.append({"scene": scene, "fold": fold, "side": side, "block": block})
            try:
                calibration = fit_block_glint(train.rhos_nir_V0, train[[f"rhos_{b}_V0" for b in ROLES]], ids)
                info["empirical"] = asdict(calibration)
                info["empirical_status"] = "available"
                empirical = calibration.slopes
            except GlintCalibrationError as error:
                empirical = None
                info["empirical_status"] = "unavailable"
                info["reason"] = str(error)
            calibrations[scene]["folds"][fold] = info
            x = test.rhos_nir_V0.to_numpy()
            y = test[[f"rhos_{b}_V0" for b in ROLES]].to_numpy()
            positions = part.index.get_indexer(test.index)
            if (positions < 0).any():
                raise RuntimeError("Evaluation rows are absent from native sampling grid")
            support = (x >= nir_range[0]) & (x <= nir_range[1])
            alternatives = {
                "V0": y,
                "V1_LUT_zero": test[[f"rhos_{b}_V1" for b in ROLES]].to_numpy(),
                "LUT_p05": native.apply(theoretical, baseline)[positions],
            }
            if empirical is not None:
                alternatives["empirical_zero"] = native.apply(empirical, 0)[positions]
                alternatives["empirical_p05"] = native.apply(empirical, baseline)[positions]
                for q in (0.01, 0.10):
                    alternatives[f"empirical_p{int(q * 100):02d}_sensitivity"] = native.apply(
                        empirical, float(np.quantile(train.rhos_nir_V0, q))
                    )[positions]
            print(scene, fold, info["empirical_status"], flush=True)
            for method, corrected in alternatives.items():
                for region in test.region.unique():
                    mask = (test.region == region).to_numpy()
                    rows = metrics(test[mask], corrected[mask], support[mask], scene, fold, method, len(np.unique(ids)))
                    records.extend([{**row, "region": region} for row in rows])
                if fold == "full_reference_preview" and method == "empirical_p05":
                    preview = test[["scene", "region", "x", "y", "depth_m", "shore_m"]].copy()
                    for i, role in enumerate(ROLES):
                        preview[f"rhos_{role}_candidate"] = corrected[:, i]
                    preview["within_nir_support"] = support
                    preview["absolute_radiometry_validated"] = False
                    preview["method"] = "empirical_p05"
                    preview["evaluation_role"] = "full_reference_development_preview"
                    preview["qualified_for_cover"] = False
                    previews.append(preview)
                ref = (test.region == "deep_reference").to_numpy()
                if not ref.any():
                    continue
                z = test.depth_m.to_numpy()[ref] + meta["water_level_m"]
                obs = corrected[ref, 1:3]
                for state in states:
                    if "a" not in state:
                        continue
                    low = floor(state["a"], state["bb"], meta["sza"], z)
                    high = floor(state["a"], state["bb"], meta["sza"], z, 0.6)
                    both = ((obs >= low - 0.001) & (obs <= high + 0.001)).all(axis=1)
                    closure.append(
                        {
                            "scene": scene,
                            "fold": fold,
                            "method": method,
                            "water_state": state["state"],
                            "cdm_conflict": state.get("cdm_conflict"),
                            "cells": len(obs),
                            "both_bands_within_fraction": float(both.mean()),
                            "median_blue_residual": np.median(obs[:, 0] - low[:, 0]),
                            "median_green_residual": np.median(obs[:, 1] - low[:, 1]),
                        }
                    )
    pd.DataFrame(records).to_csv(output / "evaluation.csv", index=False)
    pd.DataFrame(closure).to_csv(output / "water-closure.csv", index=False)
    pd.DataFrame(assignments).to_csv(output / "spatial-assignments.csv", index=False)
    if previews:
        pd.concat(previews).to_csv(output / "candidate-grid-diagnostics.csv.gz", index=False)
    write(output / "calibration.json", calibrations)
    write(output / "water-states.json", water_records)
    unchanged = all(digest(path) == value for path, value in old_files.items())
    if not unchanged:
        raise RuntimeError("A source evidence file changed during the optical review")
    write(
        output / "provenance.json",
        {
            "created_at": datetime.now(UTC).isoformat(),
            "inputs_sha256": {
                str(p): digest(p)
                for p in [
                    source,
                    spectra_path,
                    audit_path,
                    Path(design_path),
                    AUDIT / "scene-metadata.json",
                    AUDIT / "water-states.json",
                ]
                + sorted(native_inputs)
            },
            "code_sha256": {
                p.name: digest(p)
                for p in [
                    Path(__file__),
                    Path(__file__).with_name("optical_corrections.py"),
                    Path(__file__).with_name("water_constraints.py"),
                    Path(__file__).with_name("pneo_optical_report.py"),
                ]
            },
            "cover_targets_read": False,
            "cover_scores_computed": False,
            "spatial_splits_disjoint": True,
            "legacy_glint_reproduced_within_rhos": 2e-7,
            "legacy_glint_reproduction_max_abs_rhos": reconstruction_errors,
            "absolute_radiometry_validated": False,
            "historical_outputs_unchanged": unchanged,
            "frozen_r7_code_hashes_verified": len(frozen),
        },
    )
    write(
        output / "verification.json",
        {
            "historical_file_hashes_verified": len(old_files),
            "historical_outputs_unchanged": unchanged,
            "frozen_r7_code_hashes_verified": len(frozen),
            "legacy_native_correction_max_abs_rhos": reconstruction_errors,
            "cover_targets_read": False,
            "cover_scores_computed": False,
        },
    )
    write(
        output / "decision.json",
        {
            "status": "not_qualified_for_cover",
            "formal_binary_radiometric_gate_registered": False,
            "basis": "Development review: spatial baseline transfer and absolute optical closure remain unresolved.",
            "adopt_as_default": False,
            "existing_r7_verdict_changed": False,
        },
    )
    render(output)
    write(output / "outputs.json", {"sha256": {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}})
    print("Completed", output, flush=True)
    return output
