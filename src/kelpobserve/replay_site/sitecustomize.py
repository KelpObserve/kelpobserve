"""Replay hook. kelpobserve.replay puts this folder first on PYTHONPATH and names the job in
KELPOBSERVE_REPLAY_JOB, so that every interpreter of a replay, worker processes included,
rebinds the benchmark runner before any stage code runs."""

import os

if os.environ.get("KELPOBSERVE_REPLAY_JOB"):
    try:
        from kelpobserve._patch import install
    except ImportError:
        pass
    else:
        install(os.environ["KELPOBSERVE_REPLAY_JOB"])
