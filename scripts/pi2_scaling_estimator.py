#!/usr/bin/env python3
"""
This script has moved to: scripts/experimental/pi2_scaling_estimator.py

It is maintained as a DC/adiabatic diagnostic only (not for universal G calibration).
This stub forwards execution to the new location for backward compatibility.
"""
import os
import runpy


def _forward():
    here = os.path.dirname(__file__)
    new_path = os.path.join(here, "experimental", "pi2_scaling_estimator.py")
    if not os.path.exists(new_path):
        raise FileNotFoundError(f"Cannot locate experimental estimator at {new_path}")
    runpy.run_path(new_path, run_name="__main__")


if __name__ == "__main__":
    _forward()