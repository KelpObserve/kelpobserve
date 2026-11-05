"""Generate the matched-sample monthly comparison figures from held-out predictions."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kelpobserve.monthly.evaluation import score

PANELS = [
    ("median", "Median baseline", "#687581"),
    ("depth_ridge", "Depth-only ridge", "#D36B35"),
    ("spectral_stats_5_ridge", "Spectral-statistics ridge", "#087E8B"),
]


def plot_months(results: Path, output: Path) -> list[Path]:
    predictions = pd.read_csv(results / "predictions.csv")
    p = predictions[
        (predictions.subset == "matched_depth") & predictions.experiment.isin([panel[0] for panel in PANELS])
    ]
    if p.duplicated(["campaign", "experiment", "event_id"]).any():
        raise ValueError("Duplicate predictions in figure inputs")
    output.mkdir(parents=True, exist_ok=True)
    low = min(-5, np.floor(p.predicted.min() / 10) * 10 - 5)
    high = max(105, np.ceil(p.predicted.max() / 10) * 10 + 5)
    paths = []
    with plt.rc_context(
        {"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "svg.hashsalt": "kelpobserve-monthly"}
    ):
        for campaign in ["july", "september"]:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), layout="constrained", sharex=True, sharey=True)
            reference = None
            for ax, (name, title, color) in zip(axes, PANELS):
                q = p[(p.campaign == campaign) & (p.experiment == name)].sort_values("event_id")
                observed = q.set_index("event_id").observed
                if q.empty:
                    raise ValueError(f"Missing {campaign} {name} predictions")
                if reference is None:
                    reference = observed
                else:
                    pd.testing.assert_series_equal(reference, observed)
                metrics = score(q.observed.to_numpy(), q.predicted.to_numpy())
                ax.scatter(
                    q.observed, q.predicted, s=62, color=color, alpha=0.75, edgecolors="white", linewidth=0.6, zorder=3
                )
                ax.plot([0, 100], [0, 100], "--", color="#747d86", lw=1.2, label="Perfect prediction")
                ax.axhline(0, color="#bcc2c9", lw=0.7)
                ax.set(
                    xlim=(-5, 105),
                    ylim=(low, high),
                    xticks=range(0, 101, 20),
                    xlabel="Observed kelp cover (%)",
                    title=f"{title}\nMAE {metrics['mae_pp']:.2f} pp · RMSE {metrics['rmse_pp']:.2f} pp",
                )
                ax.grid(alpha=0.18)
                ax.legend(loc="upper left", frameon=False, fontsize=9)
            axes[0].set_ylabel("Leave-one-out predicted kelp cover (%)")
            fig.suptitle(
                f"{campaign.title()} · {len(reference)} quadrats · separate monthly leave-one-out evaluation",
                fontsize=16,
                fontweight="bold",
            )
            fig.text(
                0.01,
                -0.025,
                "Identical quadrats across panels; each prediction excludes that quadrat from training. Predictions are not clipped.\n"
                "Spectral statistics: original 18 features (six-band mean, standard deviation and range), 5×5 PNEO window. Errors in percentage points (pp).",
                fontsize=10,
                color="#52606d",
                va="top",
            )
            for extension in ["png", "svg"]:
                path = output / f"{campaign}_median_depth_spectral.{extension}"
                fig.savefig(path, dpi=180, bbox_inches="tight", metadata={"Date": None} if extension == "svg" else {})
                if extension == "svg":
                    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
                paths.append(path)
            plt.close(fig)
    return paths
