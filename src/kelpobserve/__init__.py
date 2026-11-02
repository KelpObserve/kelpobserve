"""Reproduction archive of the KelpObserve Sesimbra analyses (AQUASERV TA project PID 42082)."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("kelpobserve")
except PackageNotFoundError:
    __version__ = "0+unknown"
