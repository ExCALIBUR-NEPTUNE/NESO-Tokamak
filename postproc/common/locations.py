import os.path
from pathlib import Path

# Some standard locations in the repo. Assumes this file is two directories below the repo root.
LOCATIONS = {}
LOCATIONS["repo_root"] = os.path.abspath(os.path.dirname(__file__) + "/../..")
LOCATIONS["runs"] = os.path.normpath(
    os.path.join(LOCATIONS["repo_root"], "example-runs")
)
LOCATIONS["postproc"] = os.path.normpath(
    os.path.join(LOCATIONS["repo_root"], "postproc")
)


def get_run_dir(run_lbl):
    """Returns the directory path for a simulation run. Throws FileNotFoundError if the directory doesn't exist."""
    run_dir = os.path.join(LOCATIONS["runs"], run_lbl)
    if os.path.isdir(run_dir):
        return run_dir
    else:
        raise FileNotFoundError(f"No run directory at {run_dir}")


def get_postproc_dir(run_lbl):
    """Generates a directory path for postprocessing a simulation run. Creates the directory if it doesn't exist already."""
    pp_dir = os.path.join(LOCATIONS["postproc"], run_lbl)
    Path(pp_dir).mkdir(parents=True, exist_ok=True)
    if os.path.isdir(pp_dir):
        return pp_dir
    else:
        raise RuntimeError(f"Failed to create output dir at {pp_dir}")
