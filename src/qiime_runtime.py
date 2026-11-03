"""Use the interpreter selected by conda run even if a shell prepends system Python."""
import os
from pathlib import Path
import sys


def selected_python():
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix and Path(prefix).name == "qiime2-amplicon-2025.7":
        selected = Path(prefix).resolve()
        if Path(sys.prefix).resolve() != selected:
            python = selected / "bin/python3"
            os.execv(str(python), [str(python), "-B", *sys.argv])
