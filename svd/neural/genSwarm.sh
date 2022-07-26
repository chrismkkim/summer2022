#!/bin/bash
dirname="/home/kimchm/balanced/analysis/alm/svd/"
fname="script_genubal.swarm"
touch "${dirname}${fname}"
for repi in {1..8000};
    do
        echo "julia main_bal_genuStep2.jl ${repi}" >> "${dirname}${fname}"
done