#!/bin/bash
dirname="/home/kimchm/balanced/dale/fitdata/5k/almTrained/"
fname="script_run.swarm"
touch "${dirname}${fname}"
for lami in {1..5};
do
    for Li in {1..5};
    do
        echo "julia bnet.jl ${lami} ${Li}" >> "${dirname}${fname}"
    done
done
