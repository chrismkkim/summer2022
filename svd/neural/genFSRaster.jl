function genFSRaster(rate_R, rate_L, cell_R, cell_L, units_R, units_L, spiketime_R, spiketime_L, eventtime_R, eventtime_L)
    
    start_time = -2.0
    
    # #---------- remove neurons with rate < 1Hz ----------#
    # meanrate_R = mean(rate_R, dims=2)[:]
    # meanrate_L = mean(rate_L, dims=2)[:]
    # idx1Hz_R = meanrate_R .> 1.0
    # idx1Hz_L = meanrate_L .> 1.0
    # idx1Hz = idx1Hz_R .* idx1Hz_L

    # rate_R = rate_R[idx1Hz,:]
    # rate_L = rate_L[idx1Hz,:]
    # cells = cell_R[idx1Hz]

    cells = copy(cell_R)

    #---------- order neurons by selectivity ----------#
    rate_diff = mean(rate_L - rate_R, dims=2)[:]
    cellsOrd = cells[sortperm(rate_diff)]

    #---------- spike raster: lick right / left ----------#
    ns_R = zeros(Int, length(cellsOrd))
    times_R = zeros(length(cellsOrd), 300)
    ns_L = zeros(Int, length(cellsOrd))
    times_L = zeros(length(cellsOrd), 300)
    triali = 3

    for celli = 1:length(cellsOrd)
        println(celli, " / ", length(cellsOrd))

        ci = cellsOrd[celli]

        # lick right - raster
        ind_R = findall(x->x==ci, units_R)
        spiketime_ci_R = spiketime_R[ind_R]
        eventtime_ci_R = eventtime_R[ind_R]        

        if length(spiketime_ci_R) >= triali # number of recorded trials > triali
            spktime = spiketime_ci_R[triali] .- eventtime_ci_R[triali] .- start_time
            idx = (spktime .> 0).*(spktime .< 2.0)
            nspk = Int(sum(idx))
            if nspk > 0
                # println(nspk)
                # println(spktime[idx])

                ns_R[celli] = nspk
                times_R[celli, 1:nspk] .= spktime[idx]
            end
        end

        # lick left - raster
        ind_L = findall(x->x==ci, units_L)
        spiketime_ci_L = spiketime_L[ind_L]
        eventtime_ci_L = eventtime_L[ind_L]

        if length(spiketime_ci_L) >= triali # number of recorded trials > triali
            spktime = spiketime_ci_L[triali] .- eventtime_ci_L[triali] .- start_time
            idx = (spktime .> 0).*(spktime .< 2.0)
            nspk = Int(sum(idx))
            if nspk > 0
                ns_L[celli] = nspk
                times_L[celli, 1:nspk] .= spktime[idx]
            end
        end

    end


    return ns_R, times_R, ns_L, times_L, cellsOrd, rate_diff



end


