"""Commands for the separate monthly experiment; registered replays are untouched."""

from pathlib import Path


def register(commands):
    from kelpobserve.monthly.evaluation import BENCHMARK, INPUTS, RESULTS

    monthly = commands.add_parser("monthly", help="separate monthly leave-one-quadrat-out evaluation")
    actions = monthly.add_subparsers(dest="monthly_action", required=True)
    evaluate = actions.add_parser("evaluate", help="evaluate July and September independently")
    evaluate.add_argument("--inputs", type=Path, default=INPUTS)
    evaluate.add_argument("--output", type=Path, default=RESULTS)
    plot = actions.add_parser("plot", help="draw the monthly comparison figures")
    plot.add_argument("--results", type=Path, default=RESULTS)
    plot.add_argument("--output", type=Path, default=BENCHMARK / "docs" / "figures" / "monthly")
    extract = actions.add_parser("extract", help="extract native PNEO features from licensed imagery")
    extract.add_argument("--pneo-root", type=Path, required=True)
    extract.add_argument("--inputs", type=Path, default=INPUTS)
    extract.add_argument("--output", type=Path, default=BENCHMARK.parent / "out" / "monthly-extraction")
    extract.add_argument("--recalibrate", action="store_true", help="recompute unlabelled full-scene calibration")


def run(args):
    if args.monthly_action == "evaluate":
        from kelpobserve.monthly.evaluation import evaluate

        metrics = evaluate(args.inputs, args.output)
        print(
            metrics[metrics.subset == "matched_depth"][["campaign", "experiment", "n", "mae_pp", "rmse_pp"]].to_string(
                index=False
            )
        )
    elif args.monthly_action == "plot":
        from kelpobserve.monthly.plots import plot_months

        for path in plot_months(args.results, args.output):
            print(path)
    else:
        from kelpobserve.monthly.pneo import extract

        print(extract(args.pneo_root, args.inputs, args.output, args.recalibrate))
    return 0
