function runperformance(xtarg,xsum, almOrd,matchedCells)

    pcor = zeros(length(almOrd))    
    for nid = 1:length(almOrd)
        ci_alm = almOrd[nid] # alm neuron
        ci = matchedCells[nid] # model neuron

        xtarg_slice = @view xtarg[1:end-1,ci_alm]
        xsum_slice = @view xsum[:,ci]
        pcor[nid] = cor(xtarg_slice, xsum_slice)
    end

    return pcor

end