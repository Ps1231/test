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
from gee_config import SEASONS, DEFAULT_SEASON, AOI_FILE


def describe_aoi():
    """Report which AOI actually loaded.

    The old version printed a hardcoded label, so the log said nothing about
    which of the two geojson files was in use. Read it from the file.
    """
    import json
    try:
        with open(AOI_FILE) as f:
            props = json.load(f)["features"][0]["properties"]
        return (f"{props.get('districts', '?')} "
                f"({props.get('area_sq_km', 0):.0f} km2) <- {AOI_FILE.name}")
    except Exception as exc:
        return f"{AOI_FILE.name} (could not read properties: {exc})"


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
    # Track B: stack goes STRAIGHT TO LOCAL DISK — no Drive.
    # (static + meteo still use Drive; they're unchanged.)
    # Optional 3rd arg overrides the auto grid:
    #     python run_ingest.py stack rabi_2023_24 16
    from gee_config import LOCAL_STACK_DIR
    grid = int(sys.argv[3]) if len(sys.argv) > 3 else None
    written = ge.export_season_stack_local(aoi, season, str(LOCAL_STACK_DIR),
                                           grid_rows=grid, grid_cols=grid)
    print(f"\n{len(written)} stack tiles saved locally to {LOCAL_STACK_DIR}")
    print("Next: run liss3_process.py, then merge_cube.py")


def cmd_monthly(aoi, season):
    # monthly composites straight to local disk (window_days=30);
    # files are tagged "monthly" so they never clash with 8-day.
    from gee_config import LOCAL_STACK_DIR
    written = ge.export_season_stack_local(aoi, season, str(LOCAL_STACK_DIR),
                                           window_days=30)
    print(f"\n{len(written)} monthly tiles saved locally to {LOCAL_STACK_DIR}")


COMMANDS = {
    "orbits":    cmd_orbits,
    "inventory": cmd_inventory,
    "static":    cmd_static,
    "meteo":     cmd_meteo,
    "stack":     cmd_stack,
    "monthly":   cmd_monthly,
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
    print(f"AOI loaded :: {describe_aoi()}")

    COMMANDS[cmd](aoi, season)


if __name__ == "__main__":
    main()