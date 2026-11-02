"""Command line of the KelpObserve reproduction archive.

    kelpobserve experiments              # the experiments that can be replayed
    kelpobserve replay v1                # replay one and compare it with the archived runs
    kelpobserve replay r1 --stages "revision models"
    kelpobserve recompare out/replay/r1-20261018T203922Z
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from kelpobserve import config, replay


def summary(report: dict) -> dict:
    return {
        "experiment": report["experiment"],
        "checks": report["checks"],
        "stages": {label: {"reproduced": result["reproduced"],
                           "files": {name: entry["status"] for name, entry in result.get("files", {}).items()}}
                   for label, result in report["stages"].items()},
        "not_replayed": report["not_replayed"],
        "stopped": report.get("stopped"),
        "report": f"{report['workspace']}/replay.json",
    }


def kd490_command(args) -> int:
    from kelpobserve.attenuation import cmems, table

    paths = config.load()
    if args.action == "fetch":
        if args.env_file:
            for line in args.env_file.read_text().splitlines():
                key, sep, value = line.partition("=")
                if sep and key.strip() in cmems.CREDENTIALS:
                    os.environ[key.strip()] = value.strip().strip("'\"")
        recorded = json.loads(table.EVIDENCE.read_text()) if table.EVIDENCE.is_file() else None
        local = table.SUBSETS / recorded["file"] if recorded else None
        if local and local.is_file():
            ok = cmems.matches(local, recorded)
            print(f"{local} {'is' if ok else 'is not'} the recorded subset; nothing was downloaded.")
            return 0 if ok else 1
        record = cmems.fetch(table.SUBSETS)
        if recorded is None:
            cmems.write_provenance(record, table.EVIDENCE)
            print(json.dumps(record, indent=2))
            return 0
        if record["data_sha256"] == recorded.get("data_sha256") or record["sha256"] == recorded["sha256"]:
            print(f"The download matches the recorded subset; {table.EVIDENCE} is kept unchanged.")
            return 0
        other = table.EVIDENCE.with_name(f"cmems-kd490-subset.{record['fetched_at'].replace(':', '')}.json")
        cmems.write_provenance(record, other)
        print(f"The download differs from the recorded subset; its record is {other}, and {table.EVIDENCE} is kept.")
        return 1
    record = table.build(paths)
    print(json.dumps({"written": table.write(record), "scenes": record["scenes"]}, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="kelpobserve", description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("experiments", help="list the replayable experiments and their stages")
    replayer = commands.add_parser("replay", help="replay an experiment in a new workspace and compare it")
    replayer.add_argument("experiment", choices=sorted(replay.EXPERIMENTS))
    replayer.add_argument("--stages", nargs="+", help="only these stages, by label or output folder")
    again = commands.add_parser("recompare", help="compare a finished replay again, without rerunning it")
    again.add_argument("workspace", type=Path)
    kd490 = commands.add_parser("kd490", help="CMEMS KD490 as an external attenuation source")
    kd490_commands = kd490.add_subparsers(dest="action", required=True)
    fetcher = kd490_commands.add_parser("fetch", help="freeze the KD490 subset for the Sesimbra reef")
    fetcher.add_argument("--env-file", type=Path, help="file holding COPERNICUSMARINE_SERVICE_USERNAME and _PASSWORD")
    kd490_commands.add_parser("table", help="in-scene QAA, nearest CMEMS pixel and implied IOPs per scene")
    args = parser.parse_args(argv)

    if args.command == "kd490":
        return kd490_command(args)

    if args.command == "experiments":
        for name, experiment in replay.EXPERIMENTS.items():
            print(f"{name}: {experiment.title}")
            for stage in experiment.stages:
                print(f"    {stage.label}")
            for label, reason in experiment.not_replayed.items():
                print(f"    not replayed: {label} ({reason})")
        return 0
    if args.command == "recompare":
        report = replay.recompare(args.workspace.resolve(), config.load())
    else:
        report = replay.run(args.experiment, config.load(), tuple(args.stages) if args.stages else None)
    print(json.dumps(summary(report), indent=2))
    reproduced = all(result["reproduced"] for result in report["stages"].values())
    return 0 if reproduced and not report.get("stopped") else 1


if __name__ == "__main__":
    sys.exit(main())
