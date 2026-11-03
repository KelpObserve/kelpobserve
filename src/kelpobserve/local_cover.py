"""Campaign-matched, fixed-physics kelp-cover development experiment (r7).

Run with ``kelpobserve local-cover --output NEW_DIRECTORY``. The frozen protocol,
input and code hashes, all units and diagnostics accompany each run. Earlier
benchmark outputs and image products are read only.
"""

from __future__ import annotations

import contextlib
import copy
import hashlib
import importlib.metadata
import itertools
import json
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from kelpobserve import config, replay
from kelpobserve import cover_retrieval as cr
from kelpobserve.attenuation import bands

REGISTRATION = config.BENCHMARK / "registration/revision-r7.json"
OPTICAL_AUDIT = config.BENCHMARK / "registration/evidence/cmems-optical-inputs-20261102.json"
E5S = config.OPTICS / "reports/e5s-kelp-spectra/results.json"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


@contextlib.contextmanager
def legacy_readers(paths):
    """Reuse registered extraction with relocated read-only archive paths."""
    sys.path.insert(0, str(config.BENCHMARK))
    from benchmark import revision as rv
    from benchmark import run_benchmark as rb

    relocated = replay.relocation(replay.EXPERIMENTS["r1"], paths, config.BENCHMARK, archive=True)
    spec = relocated.json(json.loads((config.BENCHMARK / "registration/experiment.json").read_text()))
    original_scenes = rb.scene_directories

    def scenes(spec):
        return {
            key: {**value, "output_dir": str(relocated(value["output_dir"]))}
            for key, value in original_scenes(spec).items()
        }

    with (
        patch.object(rb, "OUTPUT", config.ARCHIVE),
        patch.object(rb, "scene_directories", scenes),
        patch.object(rv, "ACOLITE_R1", paths.require("cache") / "acolite-r1"),
    ):
        yield rb, rv, spec


def members(e5s, sensor, background, quantile="median"):
    table = e5s["class_bands_sentinel2"]["pooled"] if sensor == "sentinel2" else e5s["class_bands_pleiades_neo"]

    def spectrum(name):
        return np.array([table[name]["bands"][b][quantile] for b in ("B2", "B3")])

    kelp = spectrum("kelp")
    base = (
        np.mean([spectrum(n) for n in ("turf", "coralline_crust", "red_foliose")], axis=0)
        if background == "understory"
        else spectrum("bare")
    )
    return kelp, base


def optical_state(e5s, wavelengths, kd490, solar_zenith, shape="median"):
    e = copy.deepcopy(e5s)
    e["iops"]["median"] = e["iops"][shape]
    result = bands.at_bands(kd490, wavelengths, solar_zenith, e)
    return np.array(result["a"]), np.array(result["bb"])


def read_units(protocol, paths):
    """Extract support and targets; no physical fit or model score here."""
    import rasterio

    from kelpobserve.coastal.provenance import resolve_run

    all_units, scene_info, sources = (
        [],
        {},
        {REGISTRATION, E5S, OPTICAL_AUDIT, config.BENCHMARK / "registration/experiment.json"},
    )
    with legacy_readers(paths) as (rb, rv, spec):
        levels = rb.acquisition_levels(spec)
        sources.add(config.BENCHMARK / "registration" / spec["water_levels"]["acquisition_path"])
        sources.update(Path(spec["sources"][k]) for k in ("lidar", "dtm"))
        inventory = resolve_run(config.ARCHIVE / "inventory")
        sources.update(inventory / name for name in ("events.gpkg", "targets_development.csv", "targets_september.csv"))
        sources.update(
            [
                config.ARCHIVE / "inventory/latest.json",
                config.ARCHIVE / "prepare/latest.json",
                resolve_run(config.ARCHIVE / "prepare") / "manifest.json",
            ]
        )
        for scene in protocol["scenes"]:
            sensor, campaign, scene_id = scene["sensor"], scene["campaign"], scene["id"]
            events = rv.campaign_events(spec, campaign)
            targets = rv.targets_for(campaign)
            for variant in protocol["radiometry_variants"]:
                window = rv.read_window(spec, scene_id, rv.event_bounds(events, 100), variant)
                directory = rv.scene_dir(spec, scene_id, variant)
                rasters = rb.acolite_rasters(directory, "rhos")
                roles = spec["scene_qc"]["bands_nm"][sensor]
                names = {role: rb.nearest_band(rasters, roles[role]) for role in rv.ROLES}
                sources.update(path for _, path in names.values())
                sources.add(rv.mask_band(spec, directory)[0])
                key = scene_id + ":" + variant
                scene_info[key] = {
                    **scene,
                    "variant": variant,
                    "solar_zenith_deg": window.solar_zenith_deg,
                    "wavelengths_nm": [names[r][0] for r in ("blue", "green")],
                    "acquisition_level_m": levels[scene_id],
                    "eligible_events": len(events),
                }
                grid_transform, _, _ = rb.common_grid(spec, rv.event_bounds(events, 100), spec["selected_grid_scene"])
                cells = {
                    event.event_id: tuple(
                        map(int, rasterio.transform.rowcol(grid_transform, event.geometry.x, event.geometry.y))
                    )
                    for event in events.itertuples()
                }
                for scale in ["10m"] if sensor == "sentinel2" else ["native", "10m"]:
                    extract = rv.native_units if scale == "native" else rv.grid_units
                    units = rb.attach_targets(extract(spec, events, window, levels[scene_id], campaign), targets)
                    units["mixed_transects"] = units.get("mixed_transects", False)
                    units["mixed_strata"] = units.get("mixed_strata", False)
                    units["cell_ids"] = units.event_ids.map(
                        lambda text, cells=cells: sorted({cells[e] for e in json.loads(text)})
                    )
                    units["scene"] = scene_id
                    units["variant"] = variant
                    units["sensor"] = sensor
                    units["scale"] = scale
                    units["scene_key"] = key
                    # Explicit exclusions are retained in the unit table.
                    units["eligible"] = (
                        (units.reason == "ok")
                        & (units.depth > 0)
                        & np.isfinite(units[["blue", "green", "red", "depth", "cover_pct"]]).all(axis=1)
                        & ~units.mixed_transects
                        & ~units.mixed_strata
                    )
                    units["analysis_exclusion"] = [
                        (
                            "mixed_transects"
                            if u.mixed_transects
                            else "mixed_strata"
                            if u.mixed_strata
                            else "missing_cover"
                            if not np.isfinite(u.cover_pct)
                            else "missing_depth"
                            if not np.isfinite(u.depth)
                            else "nonpositive_depth"
                            if u.depth <= 0
                            else u.reason
                        )
                        for u in units.itertuples()
                    ]
                    all_units.append(units)
    return pd.concat(all_units, ignore_index=True), scene_info, sources


def score(predicted, observed):
    p, y = np.asarray(predicted, dtype=float), np.asarray(observed, dtype=float)
    ok = np.isfinite(p) & np.isfinite(y)
    return {
        "n": int(ok.sum()),
        "mae_pp": float(np.mean(np.abs(p[ok] - y[ok]))) if ok.any() else None,
        "bias_pp": float(np.mean(p[ok] - y[ok])) if ok.any() else None,
        "rmse_pp": float(np.sqrt(np.mean((p[ok] - y[ok]) ** 2))) if ok.any() else None,
    }


def baselines(units, protocol):
    from sklearn.linear_model import Ridge

    data = units[units.eligible].copy()
    families = {"mean": [], "depth": ["depth"], "raw": ["blue", "green"], "raw_depth": ["blue", "green", "depth"]}
    for name in families:
        units["baseline_" + name] = np.nan
    folds = cr.purged_folds(data, protocol["min_training_units"], protocol["min_training_per_stratum"])
    for fold in folds:
        if not fold["sufficient"]:
            continue
        train, test = units.loc[fold["train_indices"]], units.loc[fold["held_indices"]]
        for name, columns in families.items():
            if not columns:
                predictions = np.full(len(test), train.cover_pct.mean())
            else:
                x, y = train[columns].to_numpy(), train.cover_pct.to_numpy()
                mean, sd = x.mean(axis=0), x.std(axis=0)
                sd = np.where(sd > 1e-12, sd, 1.0)
                model = Ridge(alpha=protocol["ridge_alpha"]).fit((x - mean) / sd, y)
                predictions = model.predict((test[columns].to_numpy() - mean) / sd)
            units.loc[fold["held_indices"], "baseline_" + name] = np.clip(predictions, 0, 100)
    return units, folds


def paired_summary(frame, column, seed, draws):
    shared = frame[np.isfinite(frame[column].astype(float)) & np.isfinite(frame.baseline_depth.astype(float))]
    if shared.empty:
        return {"n": 0}
    errors = np.abs(shared[column] - shared.cover_pct) - np.abs(shared.baseline_depth - shared.cover_pct)
    by_group = pd.DataFrame({"transect": shared.transect, "delta": errors}).groupby("transect").delta.mean()
    rng = np.random.default_rng(seed)
    samples = rng.choice(by_group.to_numpy(), size=(draws, len(by_group)), replace=True).mean(axis=1)
    return {
        "n": len(shared),
        "transects": len(by_group),
        "mae_delta_physics_minus_depth_pp": float(errors.mean()),
        "equal_transect_delta_pp": float(by_group.mean()),
        "descriptive_cluster_interval_95": np.quantile(samples, [0.025, 0.975]).tolist()
        if len(by_group) >= 2
        else None,
        "interval_status": "descriptive" if len(by_group) >= 2 else "not_estimable_single_transect",
        "leave_one_transect_out_delta_pp": [
            (float(by_group.drop(g).mean()) if len(by_group) > 1 else None) for g in by_group.index
        ],
        "limitation": "At most three transect-number groups; bootstrap interval is descriptive, not confirmatory.",
    }


def estimate_units(units, info, water_match, e5s, protocol):
    """No field cover enters any physical estimate or compatibility decision."""
    output, scenarios = [], []
    fractions = np.linspace(0, 1, protocol["fraction_grid_size"])
    tolerances = protocol["rhos_tolerances"]
    main_tolerance = str(protocol["primary_rhos_tolerance"])
    kd0 = water_match["value"]
    pct = water_match["uncertainty_pct"] / 100
    wavelengths, sza = info["wavelengths_nm"], info["solar_zenith_deg"]
    states = {}
    for water_step, shape in itertools.product([-1, 0, 1], protocol["water_shapes"]):
        kd = kd0 * (1 + water_step * pct)
        if kd > 0:
            states[water_step, shape] = optical_state(e5s, wavelengths, kd, sza, shape)
    a0, bb0 = states[0, "median"]
    for unit in units.itertuples():
        if not unit.eligible:
            continue
        observed = np.array([unit.blue, unit.green])
        kelp, background = members(e5s, info["sensor"], "understory")
        curve = cr.curves(fractions, unit.depth, a0, bb0, sza, kelp, background)
        result = cr.retrieve(observed, curve, fractions, tolerances)
        row = {
            "index": unit.Index,
            "physical_cover_pct": result["cover_pct"],
            "boundary": result["boundary"],
            "residual_rhos": result["max_abs_residual_rhos"],
            "rho_blue": cr.benthic_reflectance(observed, unit.depth, a0, bb0, sza)[0],
            "rho_green": cr.benthic_reflectance(observed, unit.depth, a0, bb0, sza)[1],
        }
        for tol, interval in result["compatibility"].items():
            label = tol.replace(".", "p")
            row["compatible_" + label] = interval["any"]
            row["cover_low_" + label] = interval["min_pct"]
            row["cover_high_" + label] = interval["max_pct"]
        row["nominal_compatible"] = result["compatibility"][main_tolerance]["any"]
        row["compatible_cover_pct"] = result["compatibility"][main_tolerance]["best_pct"]
        lows, highs, bests = [], [], []
        for (water_step, shape), depth_shift, background_name, quantile in itertools.product(
            states, protocol["depth_shifts_m"], protocol["backgrounds"], protocol["spectral_quantiles"]
        ):
            depth = unit.depth + depth_shift
            if depth <= 0:
                continue
            a, bb = states[water_step, shape]
            kelp, background = members(e5s, info["sensor"], background_name, quantile)
            fit = cr.retrieve(
                observed, cr.curves(fractions, depth, a, bb, sza, kelp, background), fractions, [float(main_tolerance)]
            )
            compatible = fit["compatibility"][main_tolerance]
            if compatible["any"]:
                lows.append(compatible["min_pct"])
                highs.append(compatible["max_pct"])
            bests.append(fit["cover_pct"])
            scenarios.append(
                {
                    "index": unit.Index,
                    "water_step": water_step,
                    "shape": shape,
                    "depth_shift_m": depth_shift,
                    "background": background_name,
                    "quantile": quantile,
                    "cover_pct": fit["cover_pct"],
                    "compatible": compatible["any"],
                    "residual_rhos": fit["max_abs_residual_rhos"],
                }
            )
        row.update(
            sensitivity_low_pct=min(lows) if lows else None,
            sensitivity_high_pct=max(highs) if highs else None,
            sensitivity_compatible_scenarios=len(lows),
            sensitivity_scenarios=len(bests),
            all_scenarios_best_fit_min_pct=min(bests),
            all_scenarios_best_fit_max_pct=max(bests),
        )
        output.append(row)
    return pd.DataFrame(output).set_index("index") if output else pd.DataFrame(), scenarios


def positive_control(e5s, protocol, scene_info, water_by_date):
    errors = []
    fractions = np.linspace(0, 1, protocol["fraction_grid_size"])
    for info in scene_info.values():
        match = water_by_date[info["date"]]["KD490"]
        a, bb = optical_state(e5s, info["wavelengths_nm"], match["value"], info["solar_zenith_deg"])
        for depth, background_name in itertools.product([5.0, 10.0, 17.0], protocol["backgrounds"]):
            kelp, background = members(e5s, info["sensor"], background_name)
            curve = cr.curves(fractions, depth, a, bb, info["solar_zenith_deg"], kelp, background)
            for cover in [0.0, 0.1, 0.35, 0.7, 1.0]:
                observed = cr.curves([cover], depth, a, bb, info["solar_zenith_deg"], kelp, background)[0]
                fit = cr.retrieve(observed, curve, fractions, protocol["rhos_tolerances"])
                errors.append(abs(fit["cover_pct"] - 100 * cover))
    assert max(errors) < 0.11, "Noise-free recovery control failed"
    return {
        "n": len(errors),
        "maximum_cover_error_pp": max(errors),
        "scope": "Own-model, exact-water/depth/endmember numerical control; no independent physical validation.",
    }


def summarise_group(frame, protocol):
    eligible = frame[frame.eligible]
    result = {
        "total_units": len(frame),
        "eligible_units": len(eligible),
        "exclusions": frame.loc[~frame.eligible, "analysis_exclusion"].value_counts().to_dict(),
    }
    if eligible.empty or "physical_cover_pct" not in frame:
        return result
    nominal = eligible[eligible.nominal_compatible.eq(True)]
    result.update(
        best_fit_diagnostic=score(eligible.physical_cover_pct, eligible.cover_pct),
        compatible_units=len(nominal),
        compatible_share_of_eligible=len(nominal) / len(eligible),
        compatible_diagnostic=score(nominal.compatible_cover_pct, nominal.cover_pct),
        boundary_units=int(eligible.boundary.sum()),
        median_max_abs_residual_rhos=float(eligible.residual_rhos.median()),
        paired_vs_depth=paired_summary(eligible, "physical_cover_pct", protocol["seed"], protocol["bootstrap_draws"]),
        compatible_paired_vs_depth=paired_summary(
            nominal, "compatible_cover_pct", protocol["seed"], protocol["bootstrap_draws"]
        ),
    )
    paired = eligible[np.isfinite(eligible.baseline_depth)]
    result["matched_scores"] = {
        name: score(paired[column], paired.cover_pct)
        for name, column in [
            ("physical_best_fit", "physical_cover_pct"),
            ("depth", "baseline_depth"),
            ("mean", "baseline_mean"),
            ("raw", "baseline_raw"),
            ("raw_depth", "baseline_raw_depth"),
        ]
    }
    result["by_stratum"] = {
        name: {
            "units": len(g),
            "best_fit": score(g.physical_cover_pct, g.cover_pct),
            "compatible_units": int(g.nominal_compatible.sum()),
        }
        for name, g in eligible.groupby("stratum")
    }
    result["compatibility_by_tolerance"] = {
        str(t): int(eligible["compatible_" + str(t).replace(".", "p")].sum()) for t in protocol["rhos_tolerances"]
    }
    return result


def plots(frame, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    primary = frame[(frame.variant == "V1") & frame.eligible]
    panels = list(primary.groupby(["sensor", "campaign", "scale"], sort=False))
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
    for ax, ((sensor, campaign, scale), g) in zip(axes.flat, panels):
        for stratum, colour in [("shallow", "#128b89"), ("deep", "#e59b31")]:
            s = g[g.stratum == stratum]
            ax.scatter(s.cover_pct, s.physical_cover_pct, c=colour, label=stratum, alpha=0.8)
        ax.plot([0, 100], [0, 100], "--", color="grey", lw=1)
        ax.set_title(f"{sensor} · {campaign} · {scale}", fontsize=10)
        ax.set_xlim(-3, 103)
        ax.set_ylim(-3, 103)
        ax.set_xlabel("Diver cover (%)")
        ax.set_ylabel("Closest physical fit (%)")
        ax.grid(alpha=0.15)
    axes.flat[0].legend()
    fig.suptitle(
        "Nominal closest fits, including incompatible retrievals\nDevelopment diagnostics; V1 glint-corrected imagery"
    )
    fig.tight_layout()
    fig.savefig(output / "cover-diagnostics.png", dpi=180)
    plt.close(fig)


def report(record):
    rows = []
    for key, value in record["groups"].items():
        sensor, campaign, scale, variant = key.split("|")
        if variant != "V1":
            continue
        m = value.get("matched_scores", {})

        def number(x):
            return "—" if x is None else f"{x:.1f}"

        rows.append(
            f"| {sensor} | {campaign} | {scale} | {value['eligible_units']}/{value['total_units']} | "
            f"{value.get('compatible_units', 0)} | {m.get('physical_best_fit', {}).get('n', 0)} | {number(m.get('physical_best_fit', {}).get('mae_pp'))} | "
            f"{number(m.get('depth', {}).get('mae_pp'))} | {number(m.get('raw_depth', {}).get('mae_pp'))} | "
            f"{number(m.get('physical_best_fit', {}).get('bias_pp'))} |"
        )
    text = [
        "# Sesimbra r7: campaign-matched physical kelp-cover retrieval",
        "",
        f"Run: {record['started_at']}. Protocol frozen: {record['registration_frozen_at']}.",
        "",
        "The experiment estimates a kelp spectral mixture fraction with fixed Lee physics, depth and spectral references. July and September are evaluated against their own diver observations. No field-cover fitting enters the physical estimate. All field labels were used previously; this is development evidence.",
        "",
        (
            "No primary group produced a spectrally compatible cover retrieval. The closest-fit errors below describe rejected estimates; they do not establish kelp detection."
            if all(v.get("compatible_units", 0) == 0 for k, v in record["groups"].items() if k.endswith("|V1"))
            else "Primary retrieval availability is reported separately from closest-fit cover error below."
        ),
        "",
        "## Primary V1 results",
        "",
        "Errors below use the identical units with valid within-campaign baseline predictions. The physical column includes the closest fit even when its spectrum is incompatible; it is a diagnostic of failure as well as accuracy. Compatible counts cover all eligible units.",
        "",
        "| Sensor | Campaign | Footprint | Eligible / candidate units | Spectrally compatible | Matched units | Physical MAE (pp) | Depth MAE (pp) | Raw + depth MAE (pp) | Physical bias (pp) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
        *rows,
        "",
        "Compatibility means that at least one cover fraction reproduces both blue and green rhos within ±0.001. This is a declared residual tolerance based on the earlier clean-scene variability scale; it is not an accuracy guarantee or a confidence interval. Results at ±0.0005 and ±0.002 are retained in results.json. The absence of an application error tolerance precludes a map-ready verdict.",
        "",
        "## Retrieval and comparison",
        "",
        "The primary scene products use the existing V1 glint correction; unchanged V0 products are a fixed sensitivity. The 2 July S2 image precedes the 7–8 July dives by 5–6 days. The 31 August S2 and 3 September PNeo images are 1–2 days from the plausible September survey dates. PNeo 8 July is the shallow-survey day and one day after the deep survey. June retains its historical unavailable-scene outcome.",
        "",
        "PNeo uses its native multispectral imagery: the native sampling footprint is a 3 × 3 pixel mean (about 3.6 m), alongside 10 m aggregation. S2 uses 10 m pixels. The unchanged r1 pixel masks and fractional valid-area rules apply. Mixed-transect or mixed-stratum cells are excluded from this matched comparison, and every exclusion remains in units.csv. A quadrat mean samples part of a satellite footprint; it does not measure the complete-pixel cover.",
        "",
        "Water magnitude comes from the frozen, product-qualified same-day KD490 neighbour. A single factor scales the winter E5s non-water absorption and backscatter shape to that KD490. This remains an assumed summer spectrum. BBP443/CDM443 are retained in the source audit and are not combined into this retrieval before consistency checks.",
        "",
        "The fixed borrowed spectral references are the E5s sensor-convolved kelp medians and the mean of turf, coralline-crust and red-foliose medians. The forward model uses the actual blue/green raster wavelengths with centre-band water coefficients; nonlinear band integration remains a model approximation. Fractions are searched on 0–100% in 0.1-point steps. Their relationship to planar diver cover is conditional on the canopy and background model.",
        "",
        "Depth-only, mean, raw blue/green and raw-plus-depth baselines are fitted separately within each campaign, with ridge alpha 1 and training-only standardisation. Each fold holds one transect number across both strata and purges shared 10 m cells. Folds require eight training units and three per stratum. Three transect-number groups provide limited independence; cluster-bootstrap intervals and leave-one-group-out changes are descriptive. An interval is not estimated when only one transect contributes. In September the 10 m comparisons retain only four matched units from one transect; their errors do not summarise all 12 eligible cells.",
        "",
        "## Sensitivities and limits",
        "",
        "Each eligible unit carries scenarios combining KD490 at its central value and ± the provider percentage error bar; the median and each winter spectral shape; effective depth shifts −2/0/+2 m; understory/bare backgrounds; and coherent per-class spectral quartiles p25/median/p75. These are bounded stress scenarios without assigned probabilities. Quartile combinations do not reproduce a measured joint spectral distribution. Provider error bars have no assigned Gaussian confidence level. The depth scenarios also expose sensitivity to canopy height, but do not implement canopy radiative transfer or resolve horizontal positional errors.",
        "",
        "Each scenario preserves the two water coefficients as a shared state. Per-band reflectance compatibility uses a box tolerance and makes no independence claim. units.csv records the union of compatible fractions across scenarios, all-scenario closest-fit range, nominal residuals and unclipped effective benthic reflectance. sensitivity.csv.gz records every scenario. Wide ranges or rejected fits must not be interpreted as zero kelp.",
        "",
        "No empirical offset, endmember or optical parameter was selected against diver cover. Physics-based fractions cannot distinguish kelp species from spectrally similar vegetation without additional evidence. Water representativeness, borrowed spectra, canopy geometry, geolocation and image residuals remain limits.",
        "",
        f"The own-physics positive control recovered {record['positive_control']['n']} synthetic combinations with maximum cover error {record['positive_control']['maximum_cover_error_pp']:.3g} percentage points. This checks implementation with exact inputs. It does not validate the physical assumptions.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "kelpobserve local-cover --output out/local-cover/new-run",
        "```",
        "",
        "registration.json and freeze.json were written before real-cover scoring. input-lock.json pins every raster and reference read; code and package versions are recorded in freeze.json. Existing archived runs remain unchanged. results.json contains all V0/V1 groups, matched metrics, fold support, strata and paired diagnostics.",
        "",
        "![Nominal closest-fit cover diagnostics](cover-diagnostics.png)",
        "",
    ]
    return "\n".join(text)


def run(output: Path, registration: Path = REGISTRATION):
    registration_record = json.loads(registration.read_text())
    if registration_record["status"] != "frozen_development":
        raise ValueError("r7 must be frozen for development before scoring")
    protocol = registration_record["execution_protocol"]
    for relative, expected in registration_record["execution_code_sha256"].items():
        if digest(config.REPO / relative) != expected:
            raise ValueError(f"Code differs from the pre-score freeze: {relative}")
    for relative, expected in registration_record["execution_input_sha256"].items():
        if digest(config.REPO / relative) != expected:
            raise ValueError(f"Input differs from the pre-score freeze: {relative}")
    output.mkdir(parents=True, exist_ok=False)
    started = datetime.now(UTC).isoformat(timespec="seconds")
    write(output / "registration.json", registration_record)
    code = [
        Path(__file__),
        Path(cr.__file__),
        Path(bands.__file__),
        config.BENCHMARK / "benchmark/revision.py",
        config.BENCHMARK / "benchmark/run_benchmark.py",
        config.LIBRARY / "coastal/inversion/lee.py",
        config.LIBRARY / "coastal/optics/water.py",
        config.LIBRARY / "coastal/optics/qaa.py",
    ]
    freeze = {
        "created_at": started,
        "registration_sha256": digest(registration),
        "code_sha256": {str(p.relative_to(config.REPO)): digest(p) for p in code},
        "python": platform.python_version(),
        "packages": sorted(f"{d.metadata['Name']}=={d.version}" for d in importlib.metadata.distributions()),
    }
    write(output / "freeze.json", freeze)
    paths = config.load()
    units, info, sources = read_units(protocol, paths)
    audit = json.loads(OPTICAL_AUDIT.read_text())
    water_by_date = {t["date"]: t["variables"] for t in audit["targets"]}
    for item in info.values():
        assert water_by_date[item["date"]]["KD490"]["product_checks_passed"]
    bundle = config.REPO / ".cache/cmems-optics/20261102-v202603"
    for file in audit["manifest"]["files"]:
        p = bundle / file["file"]
        assert digest(p) == file["sha256"], p
        sources.add(p)
    sources.add(bundle / "manifest.json")
    assert digest(bundle / "manifest.json") == audit["manifest_sha256"]
    write(
        output / "input-lock.json",
        {
            "frozen_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "files": {str(p): digest(p) for p in sorted(sources)},
        },
    )
    units.to_csv(output / "support-before-scoring.csv", index=False)
    write(output / "scene-info.json", info)
    e5s = json.loads(E5S.read_text())
    control = positive_control(e5s, protocol, info, water_by_date)
    write(output / "positive-control.json", control)
    frames = []
    scenario_rows = []
    folds = {}
    groups = {}
    for key, part in units.groupby(["sensor", "campaign", "scale", "variant"], sort=False):
        label = "|".join(key)
        print("Scoring", label, flush=True)
        part, folds[label] = baselines(part.copy(), protocol)
        scene = info[part.scene_key.iloc[0]]
        physical, sensitivity = estimate_units(part, scene, water_by_date[scene["date"]]["KD490"], e5s, protocol)
        part = part.join(physical)
        scenario_rows.extend([{"group": label, **r} for r in sensitivity])
        frames.append(part)
        groups[label] = summarise_group(part, protocol)
    full = pd.concat(frames).sort_index()
    full.to_csv(output / "units.csv", index_label="index")
    pd.DataFrame(scenario_rows).to_csv(output / "sensitivity.csv.gz", index=False, compression="gzip")
    record = {
        "started_at": started,
        "completed_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "registration_frozen_at": registration_record["frozen_at"],
        "protocol": protocol,
        "groups": groups,
        "folds": folds,
        "positive_control": control,
        "verdict_scope": "Continuous development evidence; no application release gate or transfer claim.",
    }
    write(output / "results.json", record)
    plots(full, output)
    (output / "report.md").write_text(report(record))
    write(output / "outputs.json", {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()})
    return output
