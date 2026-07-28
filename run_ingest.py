"""
=========================================================
BAH 2026 — PS-6 :: Track A runner
=========================================================
Run order:

  1. python run_ingest.py orbits      <- FIRST, always
  2. python run_ingest.py inventory
  3. python run_ingest.py static
  4. python run_ingest.py meteo
  5. python run_ingest.py stack

Step 1 tells you which Sentinel-1 relative orbit to pin.
Do NOT skip it — mixing orbits corrupts the SAR time series
and the damage is invisible until classification underperforms.
=========================================================
"""

import sys

import gee_collections as gc
import gee_export as ge
from gee_config import SEASONS, DEFAULT_SEASON


def cmd_orbits(aoi, season):
    start, end = SEASONS[season]
    gc.inspect_s1_orbits(aoi, start, end, "DESCENDING")
    gc.inspect_s1_orbits(aoi, start, end, "ASCENDING")
    print("Pin the winner as S1_RELATIVE_ORBIT in gee_config.py,")
    print("then re-run everything else.\n")


def cmd_inventory(aoi, season):
    import build_inventory
    build_inventory.main([season])


def cmd_static(aoi, season):
    tasks = ge.export_static(aoi)
    ge.monitor(tasks)


def cmd_meteo(aoi, season):
    tasks = ge.export_meteorology(aoi, season)
    ge.monitor(tasks)


def cmd_stack(aoi, season):
    tasks = ge.export_season_stack(aoi, season, tiled=True)
    print(f"\n{len(tasks)} tasks queued. Monitor at "
          f"https://code.earthengine.google.com/tasks")
    ge.monitor(tasks, interval=60)


COMMANDS = {
    "orbits":    cmd_orbits,
    "inventory": cmd_inventory,
    "static":    cmd_static,
    "meteo":     cmd_meteo,
    "stack":     cmd_stack,
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print(__doc__)
        print("commands:", ", ".join(COMMANDS))
        return

    cmd = sys.argv[1]
    season = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_SEASON

    if season not in SEASONS:
        print(f"unknown season '{season}'. options: {list(SEASONS)}")
        return

    gc.initialize()
    aoi = gc.load_aoi()
    print(f"AOI loaded :: Mandya + Mysuru")

    COMMANDS[cmd](aoi, season)


if __name__ == "__main__":
    main()
