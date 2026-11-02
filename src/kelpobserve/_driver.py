"""Run one stage inside a replay: python -m kelpobserve._driver MODULE CALL ARGS...

CALL is main for modules with their own command line, which then receives ARGS, or the name of
the function to call. The replay hook (replay_site/sitecustomize.py) has already rebound the
experiment's anchors when this starts.
"""

from __future__ import annotations

import importlib
import os
import sys

from kelpobserve._patch import installed


def main() -> int:
    module, call, argv = sys.argv[1], sys.argv[2], sys.argv[3:]
    if not installed(os.environ.get("KELPOBSERVE_REPLAY_JOB", "")):
        raise SystemExit("The replay hook did not run; start replays with kelpobserve replay.")
    target = importlib.import_module(module)
    if call == "main":
        sys.argv = [module.rsplit(".", 1)[-1], *argv]
    result = getattr(target, call)()
    print(f"{module} {call}: {result}", flush=True)
    return result if isinstance(result, int) else 0


if __name__ == "__main__":
    sys.exit(main())
