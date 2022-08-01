function runperformance_nrep(xtarg,xsum, almOrd,matchedCells)

    nrep = 20
    pcor = zeros(length(almOrd))    
    for nid = 1:length(almOrd)
        ci_alm = almOrd[nid] # alm neuron
        ci = matchedCells[nid] # model neuron

        ptmp = zeros(nrep)
        for repi = 1:nrep
            xtarg_slice = @view xtarg[1:end-1,ci_alm]
            xsum_slice = @view xsum[:,ci,repi]
            ptmp[repi] = cor(xtarg_slice, xsum_slice)
        end
        pcor[nid] = mean(ptmp)
    end

    return pcor

end