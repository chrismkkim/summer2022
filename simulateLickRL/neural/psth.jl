function psth(dirsim, navg, nrep, lickRL)

    twin = 20
    Ncells = 5000
    stim_off = 1000
    train_time = 3000
    nstep = Int(train_time / twin) + 1
    rate = zeros(nstep, Ncells, nrep)

    # spikes to rate
    for repi = 1:nrep
        for avgi = 1:navg
            ii = (repi-1)*navg + avgi
            if lickRL == "right"
                times = load(dirsim * "times_R_$(ii).jld", "times")
                ns = load(dirsim * "ns_R_$(ii).jld", "ns")        
            elseif lickRL =="left"
                times = load(dirsim * "times_L_$(ii).jld", "times")
                ns = load(dirsim * "ns_L_$(ii).jld", "ns")        
            end
            for ci = 1:Ncells
                for spk = 1:ns[ci]
                    idx = Int(div(times[ci,spk], twin) + 1)
                    rate[idx, ci, repi] += 1 / (twin/1000) / navg
                end
            end
        end
    end

    timev = twin*collect(0:nstep-1)
    tidx = (timev .>= stim_off).*(timev .< train_time)
    rate = rate[tidx,:,:]

    return rate

end
