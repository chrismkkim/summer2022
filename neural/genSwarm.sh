#!/bin/bash
dirname="/home/kimchm/balanced/analysis/alm/lickright/"
fname="script_genu.swarm"
touch "${dirname}${fname}"
for repi in {1..8000};
    do
        echo "julia main_genuStep1.jl ${repi}" >> "${dirname}${fname}"
done