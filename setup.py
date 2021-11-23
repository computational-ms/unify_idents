#!/usr/bin/env python3
from setuptools import setup
import os

version_path = os.path.join(os.path.dirname(__file__), "unify_idents", "version.txt")
with open(version_path, "r") as version_file:
    unify_idents_version = version_file.read().strip()

setup(
    name="unify_idents",
    version=unify_idents_version,
    packages=[
        "unify_idents",
        "unify_idents.engine_parsers",
        "unify_idents.engine_parsers.ident",
        "unify_idents.engine_parsers.quant",
    ],
    package_dir={"unify_idents": "unify_idents"},
    #    package_data={
    #        "unify_idents": ["version.txt", "engine_parsers"],
    #        "": ["tests", "example_scripts"],
    #    },
    python_requires=">=3.8.0",
    install_requires=["pytest", "uparma"],
    description="Unify PSM idents",
    long_description="Unifying PSM identifications of different peptide search engines",
    author="M. Koesters, A. Vlasov, C. Fufezan",
    url="http://computational-ms.github.com",
    license="The MIT license",
    platforms="any that supports python 3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: SunOS/Solaris",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
)
