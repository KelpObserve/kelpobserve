"""Export existing r7 predictions and optical diagnostics at their field locations.

No fitting to field labels, retrieval changes, or edits to the frozen run.
Run with the kelpobserve virtualenv and --repo /path/to/kelpobserve --output NEW_DIR.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from kelpobserve import cover_retrieval as cr
from kelpobserve.attenuation import bands
from kelpobserve.coastal.inversion import lee
from kelpobserve.local_cover import clean, members

RUN = "benchmark/reports/r7-local-cover-20261102-r3"
VALID_COVER = {"measured", "absence_derived_zero", "absence_recorded_zero"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def export(repo, output):
    run = repo / RUN
    freeze = json.loads((run / "freeze.json").read_text())
    locked = json.loads((run / "input-lock.json").read_text())["files"]
    for rel, expected in freeze["code_sha256"].items():
        assert sha(repo / rel) == expected, f"Code changed since scoring: {rel}"
    inventory = next(Path(p).parent for p in locked if p.endswith("/events.gpkg"))
    audit_path = repo / "benchmark/registration/evidence/cmems-optical-inputs-20261102.json"
    spectra_path = repo / "optics/reports/e5s-kelp-spectra/results.json"
    inputs = [
        run / "units.csv",
        run / "scene-info.json",
        run / "registration.json",
        run / "freeze.json",
        run / "input-lock.json",
        inventory / "events.gpkg",
        inventory / "targets_development.csv",
        inventory / "targets_september.csv",
        audit_path,
        spectra_path,
    ]
    for path in inputs:
        if str(path) in locked:
            assert sha(path) == locked[str(path)], f"Input changed since scoring: {path}"
    before = {str(p): sha(p) for p in inputs}
    run_before = {str(p): sha(p) for p in run.iterdir() if p.is_file()}
    protocol = json.loads((run / "registration.json").read_text())["execution_protocol"]
    tolerance = protocol["primary_rhos_tolerance"]
    info = json.loads((run / "scene-info.json").read_text())
    audit = json.loads(audit_path.read_text())
    water_by_date = {r["date"]: r for r in audit["targets"]}
    spectra = json.loads(spectra_path.read_text())
    units = pd.read_csv(run / "units.csv")
    events = gpd.read_file(inventory / "events.gpkg", layer="events").to_crs(4326)
    assert events.event_id.is_unique
    events = events.set_index("event_id")
    targets = pd.concat([pd.read_csv(inventory / f) for f in ("targets_development.csv", "targets_september.csv")])
    assert targets.event_id.is_unique
    targets = targets.set_index("event_id")
    scene_rows, states, context_rows = [], {}, []
    for key, scene in info.items():
        match = water_by_date[scene["date"]]
        kd = match["variables"]["KD490"]
        state = bands.at_bands(kd["value"], scene["wavelengths_nm"], scene["solar_zenith_deg"], spectra)
        a, bb = np.array(state["a"]), np.array(state["bb"])
        kelp, background = members(spectra, scene["sensor"], protocol["primary_background"])
        states[key] = a, bb, kelp, background
        row = {
            "scene_key": key,
            "image_date": scene["date"],
            "sensor": scene["sensor"],
            "variant": scene["variant"],
            "solar_zenith_deg": scene["solar_zenith_deg"],
            "acquisition_level_m": scene["acquisition_level_m"],
            "kd490_m_inv": kd["value"],
            "kd490_provider_error_pct": kd["uncertainty_pct"],
            "water_input_status": "scene_assigned_KD490_scaled_winter_IOP_shape",
            "water_source_longitude": match["longitude"],
            "water_source_latitude": match["latitude"],
            "water_source_distance_to_reef_km": match["distance_km"],
            "scale_on_winter_iop_shape": state["scale_on_site_shape"],
        }
        for i, band in enumerate(("blue", "green")):
            row.update(
                {
                    f"{band}_nm": scene["wavelengths_nm"][i],
                    f"a_{band}_m_inv": a[i],
                    f"bb_{band}_m_inv": bb[i],
                    f"kd_{band}_m_inv": state["kd"][i],
                    f"kelp_endmember_{band}": kelp[i],
                    f"background_endmember_{band}": background[i],
                }
            )
        scene_rows.append(row)
        if scene["variant"] == protocol["primary_radiometry"]:
            for variable, data in match["variables"].items():
                context_rows.append(
                    {
                        "scene": scene["id"],
                        "image_date": scene["date"],
                        "variable": variable,
                        **data,
                        "longitude": match["longitude"],
                        "latitude": match["latitude"],
                        "distance_to_reef_km": match["distance_km"],
                        "used_in_r7": variable == "KD490",
                        "status": "satellite_product_candidate",
                        "reef_representativeness": match["reef_representativeness"],
                    }
                )
    scene_table = pd.DataFrame(scene_rows)
    optical = scene_table.drop(columns=["sensor", "variant"])
    enriched = units.merge(optical, on="scene_key", how="left", validate="many_to_one")
    enriched["source_unit_index"] = enriched["index"]
    enriched["fit_status"] = np.where(
        ~enriched.eligible,
        "not_scored",
        np.where(enriched.nominal_compatible == True, "compatible_diagnostic", "incompatible_diagnostic"),
    )
    enriched["primary_sensor_support"] = (enriched.variant == "V1") & (
        ((enriched.sensor == "sentinel2") & (enriched.scale == "10m"))
        | ((enriched.sensor == "pleiades_neo") & (enriched.scale == "native"))
    )
    enriched["footprint"] = np.where(enriched.scale == "10m", "10m_grid_cell", "3x3_native_pixels_approx_3p6m")
    enriched["compatibility_tolerance_rhos"] = tolerance
    checks = 0
    fractions = np.linspace(0, 1, protocol["fraction_grid_size"])
    for idx, unit in enriched[enriched.eligible].iterrows():
        a, bb, kelp, background = states[unit.scene_key]
        curve = cr.curves(fractions, unit.depth, a, bb, unit.solar_zenith_deg, kelp, background)
        observed = np.array([unit.blue, unit.green])
        fit = cr.retrieve(observed, curve, fractions, [tolerance])
        assert np.isclose(fit["cover_pct"], unit.physical_cover_pct, atol=1e-9)
        assert np.isclose(fit["max_abs_residual_rhos"], unit.residual_rhos, atol=1e-12)
        assert fit["compatibility"][str(tolerance)]["any"] == unit.nominal_compatible
        bottom = cr.benthic_reflectance(observed, unit.depth, a, bb, unit.solar_zenith_deg)
        assert np.allclose(bottom, [unit.rho_blue, unit.rho_green], atol=1e-10)
        best = cr.curves([unit.physical_cover_pct / 100], unit.depth, a, bb, unit.solar_zenith_deg, kelp, background)[0]
        floor = cr.curves([0], unit.depth, a, bb, unit.solar_zenith_deg, np.zeros(2), np.zeros(2))[0]
        inv_cos = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(unit.solar_zenith_deg)))
        kb = (inv_cos + lee._du_bottom(bb / (a + bb))) * (a + bb)
        for i, band in enumerate(("blue", "green")):
            diagnostics = {
                f"predicted_{band}_rhos": best[i],
                f"model_minus_observed_{band}_rhos": best[i] - observed[i],
                f"model_min_{band}_rhos": curve[:, i].min(),
                f"model_max_{band}_rhos": curve[:, i].max(),
                f"black_bottom_floor_{band}_rhos": floor[i],
                f"observed_below_floor_{band}": bool(observed[i] < floor[i] - tolerance),
                f"bottom_attenuation_factor_{band}": np.exp(-kb[i] * unit.depth),
                f"kb_{band}_m_inv": kb[i],
            }
            for name, value in diagnostics.items():
                if name not in enriched:
                    enriched[name] = pd.Series(
                        index=enriched.index, dtype="object" if isinstance(value, bool) else "float64"
                    )
                enriched.at[idx, name] = value
        checks += 1
    quadrat_rows = []
    for _, unit in enriched.iterrows():
        ids = json.loads(unit.event_ids)
        assert len(ids) == unit.n_quadrats and len(set(ids)) == len(ids)
        contributing = []
        for event_id in ids:
            target = targets.loc[event_id]
            valid_cover = target.cover_status in VALID_COVER and pd.notna(target.cover_pct)
            if valid_cover:
                contributing.append(float(target.cover_pct))
            event = events.loc[event_id]
            assert event.campaign == unit.campaign
            quadrat_rows.append(
                {
                    **unit.to_dict(),
                    "event_id": event_id,
                    "quadrat_longitude": event.geometry.x,
                    "quadrat_latitude": event.geometry.y,
                    "field_cover_pct": target.cover_pct,
                    "field_cover_status": target.cover_status,
                    "contributed_to_scored_cover": bool(unit.eligible and valid_cover),
                    "field_transect": event.transect,
                    "field_stratum": event.stratum,
                    "quadrat_number": event.quadrat,
                    "field_date_workbook": event.workbook_date,
                    "field_date_transects": event.transects_date,
                    "field_date_discrepancy": event.date_discrepancy,
                    "field_depth_label_m": event.workbook_depth_label_m,
                    "quadrat_lidar_depth_m": event.lidar_depth_m,
                    "field_depth_qc_status": event.depth_qc_status,
                    "field_join_class": event.join_class,
                    "field_eligible_registered": event.eligible_registered,
                }
            )
        assert len(contributing) == unit.n_scored_quadrats
        if contributing:
            assert np.isclose(np.mean(contributing), unit.cover_pct)
    quadrats = pd.DataFrame(quadrat_rows)
    assert not quadrats.duplicated(["source_unit_index", "event_id"]).any()
    primary = quadrats[quadrats.primary_sensor_support].copy()
    scored = primary[primary.contributed_to_scored_cover].copy()
    excluded = primary[~primary.contributed_to_scored_cover].copy()
    # Coordinates identify the quadrat; the prediction's support stays explicit.
    headline_columns = [
        "event_id",
        "quadrat_longitude",
        "quadrat_latitude",
        "campaign",
        "field_transect",
        "field_stratum",
        "field_date_workbook",
        "field_date_transects",
        "field_date_discrepancy",
        "sensor",
        "image_date",
        "scene",
        "variant",
        "unit_id",
        "source_unit_index",
        "footprint",
        "n_quadrats",
        "n_scored_quadrats",
        "field_cover_pct",
        "field_cover_status",
        "cover_pct",
        "physical_cover_pct",
        "fit_status",
        "compatible_cover_pct",
        "nominal_compatible",
        "boundary",
        "depth",
        "blue_nm",
        "green_nm",
        "blue",
        "green",
        "red",
        "nir",
        "predicted_blue_rhos",
        "predicted_green_rhos",
        "model_minus_observed_blue_rhos",
        "model_minus_observed_green_rhos",
        "residual_rhos",
        "compatibility_tolerance_rhos",
        "rho_blue",
        "rho_green",
        "a_blue_m_inv",
        "a_green_m_inv",
        "bb_blue_m_inv",
        "bb_green_m_inv",
        "kd490_m_inv",
        "kd_blue_m_inv",
        "kd_green_m_inv",
        "water_input_status",
        "bottom_attenuation_factor_blue",
        "bottom_attenuation_factor_green",
        "black_bottom_floor_blue_rhos",
        "black_bottom_floor_green_rhos",
        "observed_below_floor_blue",
        "observed_below_floor_green",
        "sensitivity_low_pct",
        "sensitivity_high_pct",
        "sensitivity_compatible_scenarios",
        "sensitivity_scenarios",
        "all_scenarios_best_fit_min_pct",
        "all_scenarios_best_fit_max_pct",
        "baseline_depth",
        "baseline_raw",
        "baseline_raw_depth",
        "baseline_mean",
    ]
    rename = {
        "cover_pct": "footprint_observed_mean_cover_pct",
        "physical_cover_pct": "closest_fit_cover_pct",
        "compatible_cover_pct": "compatible_fit_cover_pct",
        "depth": "footprint_water_depth_m",
        "blue": "observed_blue_rhos",
        "green": "observed_green_rhos",
        "red": "observed_red_rhos",
        "nir": "observed_nir_rhos",
        "rho_blue": "effective_bottom_blue_unclipped",
        "rho_green": "effective_bottom_green_unclipped",
        "baseline_depth": "depth_baseline_cover_pct",
        "baseline_raw": "raw_baseline_cover_pct",
        "baseline_raw_depth": "raw_depth_baseline_cover_pct",
        "baseline_mean": "mean_baseline_cover_pct",
    }
    scored = scored[headline_columns].rename(columns=rename)
    output.mkdir(parents=True, exist_ok=False)
    scored.to_csv(output / "scored-quadrats-primary.csv", index=False)
    quadrats.rename(columns=rename).to_csv(output / "all-quadrat-associations.csv", index=False)
    enriched.rename(columns=rename).to_csv(output / "footprints-all.csv", index=False)
    excluded.rename(columns=rename).to_csv(output / "unscored-quadrats-primary.csv", index=False)
    scene_table.to_csv(output / "optical-inputs-by-scene.csv", index=False)
    pd.DataFrame(context_rows).to_csv(output / "copernicus-context-by-scene.csv", index=False)
    features = [
        {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [r["quadrat_longitude"], r["quadrat_latitude"]]},
            "properties": r,
        }
        for r in scored.to_dict("records")
    ]
    write_json(output / "scored-quadrats-primary.geojson", {"type": "FeatureCollection", "features": features})
    groups = enriched[enriched.primary_sensor_support].groupby(["sensor", "campaign", "scale"])
    counts = []
    for (sensor, campaign, scale), group in groups:
        subset = scored[(scored.sensor == sensor) & (scored.campaign == campaign)]
        counts.append(
            {
                "sensor": sensor,
                "campaign": campaign,
                "scale": scale,
                "candidate_footprints": len(group),
                "scored_footprints": int(group.eligible.sum()),
                "scored_quadrats": len(subset),
                "compatible_footprints": int((group.nominal_compatible == True).sum()),
            }
        )
    assert scored.closest_fit_cover_pct.notna().all()
    assert scored.compatible_fit_cover_pct.isna().all()
    assert (scored.fit_status == "incompatible_diagnostic").all()
    assert {str(p): sha(p) for p in inputs} == before
    assert {str(p): sha(p) for p in run.iterdir() if p.is_file()} == run_before
    verification = {
        "source_run": str(run),
        "created_at": datetime.now(UTC).isoformat(),
        "source_and_run_files_unchanged": True,
        "frozen_code_hashes_verified": len(freeze["code_sha256"]),
        "predictions_residuals_compatibility_and_benthic_reflectance_recomputed": checks,
        "unit_cover_means_and_event_memberships_verified": len(units),
        "primary_scored_quadrat_scene_rows": len(scored),
        "unique_primary_scored_quadrats": scored.event_id.nunique(),
        "primary_unscored_quadrat_scene_rows": len(excluded),
        "all_variant_quadrat_scene_scale_rows": len(quadrats),
        "primary_groups": counts,
        "inputs_sha256": before,
        "export_script_sha256": sha(__file__),
    }
    write_json(output / "verification.json", verification)
    print(json.dumps({k: v for k, v in verification.items() if k != "inputs_sha256"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.repo.resolve(), args.output.resolve())
