"""
Shared helpers for parsing OpenFOAM field and summary files from
validate scripts. Imported by individual `tests/<case>/validate`
scripts via a sys.path insert (see any in-tree validate for the
two-line idiom).

No external dependencies beyond the Python standard library; runs
under the system /usr/bin/python3 provided by the CI image without
needing the OpenFOAM Python bindings.
"""

import os
import re
import sys


# Matches OpenFOAM nonuniform-list scalar internalField blocks:
#
#     internalField   nonuniform List<scalar>
#     1234
#     (
#     1.23 4.56 ...
#     )
#
# Group 1 is the (textual) cell count; group 2 is the contents
# between the parentheses, ready to be split on whitespace.
_NONUNIFORM_SCALAR_RE = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(([^)]+)\)",
    re.S,
)

_UNIFORM_SCALAR_RE = re.compile(
    r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;",
)


def read_internal_scalar_field(path):
    """Parse a scalar internalField from an OF field file.

    Returns the per-cell values as a list of floats. For a non-uniform
    field that's the full per-cell list in mesh order. For a uniform
    field it's a one-element list `[value]` -- callers that compute
    mean/min/max get the right answer with no special-casing. Exits
    the validate with a clear error if the file contains neither.
    """
    with open(path) as f:
        txt = f.read()
    m = _NONUNIFORM_SCALAR_RE.search(txt)
    if m:
        return list(map(float, m.group(2).split()))
    m = _UNIFORM_SCALAR_RE.search(txt)
    if m:
        return [float(m.group(1))]
    sys.exit(f"validate: could not parse scalar internalField at {path}")


def read_uniform_scalar(path):
    """Parse a uniform scalar internalField from an OF field file.

    Returns the single scalar value. Exits the validate with a clear
    error if the file does not contain a uniform scalar internalField.
    """
    with open(path) as f:
        txt = f.read()
    m = _UNIFORM_SCALAR_RE.search(txt)
    if not m:
        sys.exit(f"validate: could not parse uniform scalar field at {path}")
    return float(m.group(1))


def find_latest_time(base_dir):
    """Find the lexically-largest time-named subdir of `base_dir`.

    Time dirs are subdirectories whose name parses as a positive
    float ("0", "1", "0.5", etc.). Returns the absolute path of the
    largest, or None if no time subdirs exist.

    Used by validates that read function-object output from
    `postProcessing/<name>/<time>/...` where the function object's
    write step may run at any case time.
    """
    if not os.path.isdir(base_dir):
        return None
    candidates = []
    for entry in os.listdir(base_dir):
        try:
            t = float(entry)
        except ValueError:
            continue
        if t < 0:
            continue
        candidates.append((t, entry))
    if not candidates:
        return None
    candidates.sort()
    return os.path.join(base_dir, candidates[-1][1])


def parse_summary(path):
    """Parse an OpenFOAM-style key-value summary.dat file.

    File shape:
        # comments
        key1   value1
        key2   value2
        ...

    Returns a dict mapping keys to floats. Lines starting with '#'
    or empty lines are skipped; lines whose value column does not
    parse as a float are skipped silently (keeps the parser permissive
    for entries that mix in non-numeric tokens).
    """
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2 or parts[0].startswith("#"):
                continue
            try:
                out[parts[0]] = float(parts[-1])
            except ValueError:
                pass
    return out
