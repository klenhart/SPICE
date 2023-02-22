#!/bin/env python

#######################################################################
# Copyright (C) 2022 Christian, Blümel, Julian Dosch
#
# This file is part of grand-trumpet.
#
#  grand-trumpet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  grand-trumpet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with grand-trumpet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="grand-trumpet",
    version="0.1",
    python_requires='>=3.9.0',
    description="Genome wide appliance and visualization of the FAS algorithm for Alternative Splicing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Christian Blümel",
    author_email="christianbluemel@stud.uni-frankfurt.de",
    url="",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'greedyFAS',
        'pyranges',
        'requests',
        'sys',
        'json',
        'os',
        'argparse',
        'itertools',
        'plotly',
        'pathlib',
        'tqdm'
    ],
    entry_points={
        'console_scripts': ["gt.run = grand-trumpet.fas_lib:main",
                            "gt.run = grand-trumpet.fas_compare:main",
                            "gt.run = grand-trumpet.fas_bashAssist:main",
                            "gt.run = grand-trumpet.get_domain_importance:main",
                            "gt.run = grand-trumpet.fas_handler:main",
                            "gt.run = grand-trumpet.fas_movement:main"],
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