"""
=========================================================
BAH 2026 — PS-6 :: Track B :: Drive -> local downloader
=========================================================
GEE can only export to Google Drive or GEE Assets — it cannot write to
your local disk. So the stack is exported to Drive as usual, and THIS
module pulls those tiles down into a local folder afterwards.

Scope: STACK ONLY. Static and meteo exports are left on Drive, untouched.
Track B needs the stack on disk because merge_cube.py combines it with the
locally-processed LISS-III output.

Two ways to authenticate to Drive:
  1. PyDrive2 with the same Google account you use for GEE (interactive).
  2. A service-account JSON that has access to the Drive folder.

If you'd rather not use the Drive API at all, you can always just download
the tiles by hand from drive.google.com/BAH2026_GEE into LOCAL_STACK_DIR —
merge_cube.py only cares that the files are there.
"""

import os
import glob
from gee_config import DRIVE_FOLDER, LOCAL_STACK_DIR


def download_stack_tiles(season, mode="8day", local_dir=None,
                         drive_folder=DRIVE_FOLDER):
    """
    Download this season's stack tiles (<season>_<mode>_r*c*.tif) from the
    Drive export folder into local_dir. Returns the list of local paths.

    Uses PyDrive2. Install once:  pip install pydrive2
    """
    local_dir = str(local_dir or LOCAL_STACK_DIR)
    os.makedirs(local_dir, exist_ok=True)

    try:
        from pydrive2.auth import GoogleAuth
        from pydrive2.drive import GoogleDrive
    except ImportError:
        raise SystemExit(
            "pydrive2 not installed. Either run  pip install pydrive2\n"
            "or download the tiles manually from Drive into:\n  " + local_dir)

    gauth = GoogleAuth()
    gauth.LocalWebserverAuth()          # opens a browser once, caches creds
    drive = GoogleDrive(gauth)

    # find the Drive folder id
    folders = drive.ListFile({
        "q": f"title='{drive_folder}' and mimeType='application/vnd.google-apps.folder' and trashed=false"
    }).GetList()
    if not folders:
        raise SystemExit(f"Drive folder '{drive_folder}' not found.")
    folder_id = folders[0]["id"]

    prefix = f"{season}_{mode}_"
    files = drive.ListFile({
        "q": f"'{folder_id}' in parents and trashed=false"
    }).GetList()
    tiles = [f for f in files
             if f["title"].startswith(prefix) and f["title"].endswith(".tif")]

    if not tiles:
        print(f"No tiles matching {prefix}*.tif in Drive/{drive_folder} yet.")
        print("Has the export finished? Check the GEE Tasks page.")
        return []

    local_paths = []
    for f in tiles:
        dest = os.path.join(local_dir, f["title"])
        print(f"  downloading {f['title']} ...")
        f.GetContentFile(dest)
        local_paths.append(dest)

    print(f"Downloaded {len(local_paths)} stack tile(s) -> {local_dir}")
    return local_paths


def local_tiles_present(season, mode="8day", local_dir=None):
    """Quick check: are the stack tiles already on disk?"""
    local_dir = str(local_dir or LOCAL_STACK_DIR)
    return sorted(glob.glob(os.path.join(local_dir, f"{season}_{mode}_r*c*.tif")))


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("usage: python download_stack.py <season> [mode]")
        raise SystemExit(1)
    season = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "8day"
    download_stack_tiles(season, mode)
