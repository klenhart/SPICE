#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
"""
Created on Tue Aug  2 11:13:06 2022
@author: chrisbl
"""

import plotly.graph_objects as go
import math
from valves.fas_utility import calc_rmsd
    

def scale_list(scale_factor, float_list):
    output_list = []
    for value in float_list:
        output_list.append(value * scale_factor)
    return output_list


def make_graph(fas_lib, gene_id, sample_names, categories, sigma_list, rmsd, filepath):
    fig = go.Figure()
    for i, name in enumerate(sample_names):
        fig.add_trace(go.Scatterpolar(
            r=sigma_list[i] + [sigma_list[i][0]],
            theta=categories + [categories[0]],
            fill="toself",
            name=name
            ))
    fig.add_annotation(x=0, y=-0.2,
                text="RMSD=" + str(rmsd),
                showarrow=False)
    fig.add_annotation(x=0, y=-0.14,
                text=gene_id,
                showarrow=False)
    fig.update_layout(
        autosize=False,
        width=604,
        height=387,
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True,
        paper_bgcolor="LightSteelBlue"
        )
    #fig.show()
    fig.write_image(file=filepath)
    

def main():
    pass


if __name__ == "__main__":
    main()
    