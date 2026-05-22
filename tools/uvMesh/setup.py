"""Compatibility shim for `pip install` with older pip versions.

The package metadata canonically lives in `pyproject.toml` (PEP 621). pip's
build-isolation venv pulls in a pinned-old setuptools that doesn't read
PEP 621, leading to an "UNKNOWN" wheel name. Repeating the metadata here
makes the install work regardless of build-time setuptools version.
"""
from setuptools import find_packages, setup

setup(
    name="uvmesh",
    version="0.2.0",
    description="Hybrid O-grid-annulus + polyhedral-bulk mesh generator for UV reactor cases",
    packages=find_packages(include=["uvmesh", "uvmesh.*"]),
    python_requires=">=3.8",
    install_requires=["numpy"],
    extras_require={
        "tests": ["pytest>=7"],
    },
)
