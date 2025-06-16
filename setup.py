#!/usr/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Blümel, Julian Dosch
#
# This file is part of Spice.
#
#  Spice is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Spice is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="spice",
    version="1.0",
    python_requires='>=3.9.0',
    description="Splicing-based protein isoform comparison estimator.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Christian Blümel",
    author_email="christianbluemel@stud.uni-frankfurt.de",
    url="",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # Make sure these are on PyPI
        'greedyFAS',         # <-- If not on PyPI, consider using a GitHub URL
        'pyranges',
        'requests',
        'plotly',
        'tqdm',
        'numpy',
        'matplotlib',
        'scipy',
        'pyyaml',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            "spice.library = spice.spice_library:main",
            "spice.result = spice.spice_result:main",
            "spice.novel = spice.spice_novel:main",
            "spice.makejobs = spice.FASJobAssistant:main",
            "spice.parse_domain_out = spice.parse_domain_out:main"
        ],
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
