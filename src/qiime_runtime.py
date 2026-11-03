"""Use the interpreter selected by conda run even if a shell prepends system Python."""
import os
from pathlib import Path
import sys


def selected_python():
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix and Path(prefix).name == "qiime2-amplicon-2025.7":
        selected = Path(prefix).resolve()
        # A login shell can prepend system R as well as system Python.
        # Keep subprocess executables in the already selected conda environment.
        bindir = str(selected / "bin")
        path = os.environ.get("PATH", "")
        if path.split(os.pathsep)[0] != bindir:
            os.environ["PATH"] = bindir + os.pathsep + path
        if Path(sys.prefix).resolve() != selected:
            python = selected / "bin/python3"
            os.execv(str(python), [str(python), "-B", *sys.argv])
