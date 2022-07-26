function load_fs_1trials(rate_R, rate_L, cell_R, cell_L, units_R, units_L, spiketime_R, spiketime_L, eventtime_R, eventtime_L)
    
    start_time = -2.0
    
    ntrials = 30

    cells = copy(cell_R)

    #---------- order neurons by selectivity ----------#
    rate_diff = mean(rate_R - rate_L, dims=2)[:]
    cellsOrd = cells[sortperm(rate_diff)]

    #---------- spike raster: lick right / left ----------#
    ns_R = zeros(Int, length(cellsOrd), ntrials)
    ns_L = zeros(Int, length(cellsOrd), ntrials)

    for triali = 1:ntrials
        println(triali, " / ", length(ntrials))

        for celli = 1:length(cellsOrd)

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

                    ns_R[celli, triali] = nspk
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
                    ns_L[celli, triali] = nspk
                end
            end

        end
    end

    # correlation: 1 trial vs. psth
    pearcorr = zeros(ntrials)
    rate_diff_ord = rate_diff[sortperm(rate_diff)]
    rate_diff_1trial = (ns_R - ns_L) / 2
    for ii = 1:ntrials        
        pearcorr[ii] = cor(rate_diff_ord, rate_diff_1trial[:,ii])
    end


    return pearcorr

end




function load_pyr_1trials(rate_R, rate_L, cell_R, cell_L, units_R, units_L, spiketime_R, spiketime_L, eventtime_R, eventtime_L)
    
    start_time = -2.0
    
    ntrials = 30

    #---------- remove neurons with rate < 1Hz ----------#
    meanrate_R = mean(rate_R, dims=2)[:]
    meanrate_L = mean(rate_L, dims=2)[:]
    idx1Hz_R = meanrate_R .> 1.0
    idx1Hz_L = meanrate_L .> 1.0
    idx1Hz = idx1Hz_R .* idx1Hz_L

    rate_R = rate_R[idx1Hz,:]
    rate_L = rate_L[idx1Hz,:]
    cells = cell_R[idx1Hz]
    
    #---------- order neurons by selectivity ----------#
    rate_diff = mean(rate_R - rate_L, dims=2)[:]
    cellsOrd = cells[sortperm(rate_diff)]

    #---------- spike raster: lick right / left ----------#
    ns_R = zeros(Int, length(cellsOrd), ntrials)
    ns_L = zeros(Int, length(cellsOrd), ntrials)

    for triali = 1:ntrials
        println(triali, " / ", length(ntrials))

        for celli = 1:length(cellsOrd)

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

                    ns_R[celli, triali] = nspk
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
                    ns_L[celli, triali] = nspk
                end
            end

        end
    end

    # correlation: 1 trial vs. psth
    pearcorr = zeros(ntrials)
    rate_diff_ord = rate_diff[sortperm(rate_diff)]
    rate_diff_1trial = (ns_R - ns_L) / 2
    for ii = 1:ntrials        
        pearcorr[ii] = cor(rate_diff_ord, rate_diff_1trial[:,ii])
    end


    return pearcorr



end



