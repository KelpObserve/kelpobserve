"""Fixed-model, separate-month leave-one-quadrat-out cover evaluation.

No model selection or tuning occurs here. All numerical preprocessing is fitted
inside each outer fold. Scene-wide, unlabelled Lyzenga calibration is frozen.
"""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

BENCHMARK = Path(__file__).resolve().parents[3] / "benchmark"
INPUTS = BENCHMARK / "inputs" / "monthly"
RESULTS = BENCHMARK / "results" / "monthly"
BANDS = ["R", "G", "B", "NIR", "RE", "DB"]
PAIRS = [("DB", "B"), ("DB", "G"), ("B", "G"), ("G", "R"), ("B", "R"), ("DB", "R")]
MEANS = [f"{band}_mean" for band in BANDS]
STATS = [f"{band}_{stat}" for band in BANDS for stat in ["mean", "std", "range"]]
LYZENGA = [f"DII_{a}_{b}" for a, b in PAIRS]
FEATURES = {
    "band_means": MEANS,
    "lyzenga": LYZENGA,
    "means_lyzenga": MEANS + LYZENGA,
    "spectral_stats": STATS,
    "stats_lyzenga": STATS + LYZENGA,
    "depth": ["depth"],
    "none": [],
}
KEYS = ["campaign", "subset", "experiment"]


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


@dataclass(frozen=True)
class Experiment:
    name: str
    regressor: str
    features: str = "none"
    window_px: int = 0
    background: str = "raw"


def experiments(subset: str) -> list[Experiment]:
    """The combinations explored in this study, not an exhaustive grid search."""
    result = [Experiment("mean", "mean"), Experiment("median", "median")]
    if subset != "all_cover":
        result.append(Experiment("depth_ridge", "ridge", "depth"))
    if subset == "matched_depth":
        return result + [
            Experiment("spectral_stats_5_ridge", "ridge", "spectral_stats", 5),
            Experiment("band_means_1_ridge", "ridge", "band_means", 1),
        ]
    for window in [1, 3, 5]:
        for feature in ["band_means", "lyzenga", "means_lyzenga", "spectral_stats", "stats_lyzenga"]:
            if window == 1 and feature in ["spectral_stats", "stats_lyzenga"]:
                continue  # Single-pixel std/range are zero; these would duplicate mean-only fits.
            result.append(Experiment(f"{feature}_{window}_ridge", "ridge", feature, window))
    if subset == "all_cover":
        for regressor in ["random_forest", "extra_trees", "svr", "pls"]:
            result.append(Experiment(f"spectral_stats_5_{regressor}", regressor, "spectral_stats", 5))
        for background in ["global", "local"]:
            for regressor in ["ridge", "extra_trees"]:
                result.append(
                    Experiment(f"spectral_stats_5_{regressor}_{background}", regressor, "spectral_stats", 5, background)
                )
    return result


def make_model(name: str):
    if name == "ridge":
        return make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    if name in ["random_forest", "extra_trees"]:
        cls = RandomForestRegressor if name == "random_forest" else ExtraTreesRegressor
        return cls(n_estimators=300, max_depth=3, min_samples_leaf=5, max_features=1.0, random_state=42, n_jobs=1)
    if name == "svr":
        return make_pipeline(StandardScaler(), SVR(kernel="rbf", C=10.0, epsilon=5.0, gamma="scale"))
    if name == "pls":
        return make_pipeline(StandardScaler(), PLSRegression(n_components=2, scale=False, max_iter=1000, tol=1e-6))
    raise ValueError(f"Unknown regressor: {name}")


def load_inputs(folder: Path = INPUTS) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = json.loads((folder / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if digest(folder / name) != expected:
            raise ValueError(f"Input checksum mismatch: {name}")
    observations = pd.read_csv(folder / "observations.csv")
    features = pd.read_csv(folder / "features.csv")
    if observations.event_id.duplicated().any() or features.duplicated(["event_id", "window_px"]).any():
        raise ValueError("Duplicate quadrat or window identifiers")
    if set(observations.campaign) != {"july", "september"}:
        raise ValueError("Expected July and September campaigns only")
    if not observations.quality_eligible.isin([True, False]).all():
        raise ValueError("Quality eligibility must be boolean")
    if not observations.cover_pct.dropna().between(0, 100).all():
        raise ValueError("Cover must be in [0,100] or missing")
    if not np.array_equal(observations.cover_pct.eq(0), observations.presence.eq(0)):
        raise ValueError("Zero cover and explicit absence disagree")
    for window in [1, 3, 5]:
        q = features[features.window_px == window]
        if set(q.event_id) != set(observations.event_id):
            raise ValueError(f"Incomplete quadrat set for {window}x{window}")
    if not np.isfinite(features[STATS + LYZENGA + ["row", "column"]]).all().all():
        raise ValueError("Non-finite satellite features")
    if observations.depth.isin([np.inf, -np.inf]).any():
        raise ValueError("Non-finite depth")
    return observations, features


def select_samples(observations: pd.DataFrame, campaign: str, subset: str) -> pd.DataFrame:
    q = observations[(observations.campaign == campaign) & observations.cover_pct.notna()]
    if subset in ["matched_depth", "quality_filtered"]:
        q = q[q.depth.notna()]
    if subset == "quality_filtered":
        q = q[q.quality_eligible]
    return q.sort_values("event_id").reset_index(drop=True)


def folds(q: pd.DataFrame):
    if q.campaign.nunique() != 1 or q.event_id.duplicated().any() or len(q) < 3:
        raise ValueError("LOO requires unique quadrats from exactly one month, with at least three samples")
    for held_out in range(len(q)):
        yield np.delete(np.arange(len(q)), held_out), held_out


def subtract_background(X, y_train, groups, train, test, kind):
    """Only y_train is accepted: held-out cover cannot influence background selection."""
    controls = train[y_train == 0]
    if not len(controls):
        raise ValueError("No empty training quadrats for background subtraction")
    mean_columns = [STATS.index(name) for name in MEANS]
    a, b = X[train].copy(), X[[test]].copy()
    audit = {"global_controls": controls.tolist(), "training_references": {}}
    if kind == "global":
        bg = X[controls][:, mean_columns].mean(axis=0)
        a[:, mean_columns] -= bg
        b[:, mean_columns] -= bg
        audit.update(test_controls=controls.tolist(), local_used=False)
    elif kind == "local":
        for j, index in enumerate([*train, test]):
            pool = controls[controls != index]
            local = pool[groups[pool] == groups[index]]
            chosen = local if len(local) else pool
            if not len(chosen):
                raise ValueError("No leave-self-out empty reference for a training row")
            bg = X[chosen][:, mean_columns].mean(axis=0)
            if index == test:
                b[0, mean_columns] -= bg
                audit.update(test_controls=chosen.tolist(), local_used=bool(len(local)))
            else:
                a[j, mean_columns] -= bg
                audit["training_references"][int(index)] = chosen.tolist()
    else:
        raise ValueError(f"Unknown background: {kind}")
    return a, b, audit


def predict_experiment(q: pd.DataFrame, features: pd.DataFrame, experiment: Experiment):
    cols = FEATURES[experiment.features]
    if experiment.window_px:
        table = features[features.window_px == experiment.window_px].set_index("event_id").loc[q.event_id]
        X = table[cols].to_numpy()
    else:
        X = q[cols].to_numpy()
    if not np.isfinite(X).all():
        raise ValueError(f"Missing features for {experiment.name}")
    y = q.cover_pct.to_numpy()
    groups = (q.stratum.astype(str) + "_T" + q.transect.astype(str)).to_numpy()
    predictions = np.empty(len(q))
    audits = []
    for train, test in folds(q):
        if experiment.regressor in ["mean", "median"]:
            statistic = np.mean if experiment.regressor == "mean" else np.median
            predictions[test] = statistic(y[train])
            continue
        a, b = X[train], X[[test]]
        if experiment.background != "raw":
            a, b, audit = subtract_background(X, y[train], groups, train, test, experiment.background)
            audits.append(
                {
                    "held_out": q.event_id.iloc[test],
                    "global_controls": q.event_id.iloc[audit["global_controls"]].tolist(),
                    "test_controls": q.event_id.iloc[audit["test_controls"]].tolist(),
                    "local_used": audit["local_used"],
                    "training_references": {
                        q.event_id.iloc[i]: q.event_id.iloc[ids].tolist()
                        for i, ids in audit["training_references"].items()
                    },
                }
            )
        model = make_model(experiment.regressor)
        model.fit(a, y[train])
        predictions[test] = np.asarray(model.predict(b)).reshape(-1)[0]
    if not np.isfinite(predictions).all():
        raise ValueError(f"Non-finite predictions from {experiment.name}")
    return predictions, audits


def score(y, predicted) -> dict:
    return {
        "n": len(y),
        "mae_pp": mean_absolute_error(y, predicted),
        "rmse_pp": float(np.sqrt(mean_squared_error(y, predicted))),
        "r2": r2_score(y, predicted),
        "bias_pp": float(np.mean(predicted - y)),
        "clipped_mae_pp": mean_absolute_error(y, np.clip(predicted, 0, 100)),
        "out_of_bounds": int(((predicted < 0) | (predicted > 100)).sum()),
    }


def evaluate(inputs: Path = INPUTS, output: Path = RESULTS) -> pd.DataFrame:
    observations, features = load_inputs(inputs)
    output.mkdir(parents=True, exist_ok=True)
    metrics, predictions, split_records, background_records, overlaps, selections = [], [], [], [], [], []
    for campaign in ["july", "september"]:
        subsets = ["all_cover", "matched_depth"] + (["quality_filtered"] if campaign == "september" else [])
        for subset in subsets:
            q = select_samples(observations, campaign, subset)
            y = q.cover_pct.to_numpy()
            for row in observations[observations.campaign == campaign].itertuples():
                reasons = []
                if pd.isna(row.cover_pct):
                    reasons.append("missing_cover")
                if subset != "all_cover" and pd.isna(row.depth):
                    reasons.append("missing_depth")
                if subset == "quality_filtered" and not row.quality_eligible:
                    reasons.append("quality_flag")
                selections.append(
                    {
                        "campaign": campaign,
                        "subset": subset,
                        "event_id": row.event_id,
                        "included": not reasons,
                        "exclusion": ";".join(reasons),
                    }
                )
            for train, test in folds(q):
                split_records.append(
                    {
                        "campaign": campaign,
                        "subset": subset,
                        "held_out": q.event_id.iloc[test],
                        "training": q.event_id.iloc[train].tolist(),
                    }
                )
            for window in [1, 3, 5]:
                pixels = features[features.window_px == window].set_index("event_id").loc[q.event_id]
                for event, row in pixels.iterrows():
                    other = pixels.drop(event)
                    count = (
                        (other.row - row.row).abs().lt(window) & (other.column - row.column).abs().lt(window)
                    ).sum()
                    overlaps.append(
                        {
                            "campaign": campaign,
                            "subset": subset,
                            "window_px": window,
                            "event_id": event,
                            "overlapping_training_quadrats": int(count),
                        }
                    )
            for experiment in experiments(subset):
                predicted, audits = predict_experiment(q, features, experiment)
                info = {"campaign": campaign, "subset": subset, "experiment": experiment.name}
                metrics.append({**info, **asdict(experiment), **score(y, predicted)})
                predictions.extend(
                    {
                        **info,
                        "event_id": event,
                        "observed": float(actual),
                        "predicted": float(pred),
                        "train_n": len(q) - 1,
                    }
                    for event, actual, pred in zip(q.event_id, y, predicted)
                )
                background_records.extend({**info, **record} for record in audits)
            print(f"{campaign}: {subset}, {len(q)} quadrats, {len(experiments(subset))} experiments", flush=True)
    m, p = pd.DataFrame(metrics), pd.DataFrame(predictions)
    if p.duplicated(KEYS + ["event_id"]).any():
        raise ValueError("Duplicate held-out predictions")
    m.to_csv(output / "metrics.csv", index=False)
    p.to_csv(output / "predictions.csv", index=False)
    pd.DataFrame(selections).to_csv(output / "sample_selection.csv", index=False)
    pd.DataFrame(overlaps).to_csv(output / "window_overlap.csv", index=False)
    write_json(output / "folds.json", split_records)
    write_json(output / "background_audit.json", background_records)
    write_json(
        output / "protocol.json",
        {
            "evaluation": "Separate July and September leave-one-quadrat-out; no cross-month training or tuning.",
            "target": "CCMAR quadrat cover percent, not a fraction of presence points.",
            "input_manifest_sha256": digest(inputs / "manifest.json"),
            "models": {
                name: repr(make_model(name)) for name in ["ridge", "random_forest", "extra_trees", "svr", "pls"]
            },
            "baselines": "Training-fold mean and median, recalculated for each held-out quadrat.",
            "preprocessing": "StandardScaler inside each fold; trees unscaled. No clipping for primary metrics.",
            "calibration": "Frozen unlabelled full-scene low-NIR calibration includes held-out image areas, not cover labels.",
            "subsets": {
                "all_cover": "All observed covers; satellite models and baselines.",
                "matched_depth": "Common cover+depth observations for mean, median, depth, 5x5 stats and single pixel.",
                "quality_filtered": "September sensitivity: measured cover+depth and existing inventory quality eligibility.",
            },
            "background": "Empty training covers only. Global per-band mean subtraction; local stratum/transect reference with global fallback and leave-self-out training references. Only band means are shifted.",
            "limitations": [
                "Small spatially dependent samples; overlapping windows; no new-location or new-date validation.",
                "Multiple models examined on repeatedly inspected outcomes; exploratory comparisons, no validated winner.",
                "Depth/position QC flags retained in primary samples and exposed by quality-filtered sensitivity.",
                "Depth stratum and transect indirectly inform local background preprocessing.",
                "Native 1.2m pixel exceeds 0.5m quadrat; 3x3 and 5x5 span 3.6m and 6m.",
                "July REFLECTANCE/RAYLEIGH and September BASIC/NONE products differ radiometrically.",
            ],
            "versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "sklearn": sklearn.__version__,
            },
        },
    )
    return m
