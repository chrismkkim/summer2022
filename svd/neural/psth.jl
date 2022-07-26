function psth(dirsim, nrep)

    twin = 20
    Ncells = 5000
    stim_off = 1000
    train_time = 3000
    nstep = Int(train_time / twin) + 1
    rate = zeros(nstep, Ncells, nrep)

    # load data
    tmp = load(dirsim * "times_rep1.jld", "times")
    times = zeros(size(tmp)[1], size(tmp)[2], nrep)
    ns = zeros(Int,Ncells, nrep)
    for repi = 1:nrep
        println(repi)
        times[:,:,repi] = load(dirsim * "times_rep$(repi).jld", "times")
        ns[:,repi] = load(dirsim * "ns_rep$(repi).jld", "ns")
    end

    # spikes to rate
    for repi = 1:nrep
        for ci = 1:Ncells
            for spk = 1:ns[ci, repi]
                idx = Int(div(times[ci,spk,repi], twin) + 1)
                rate[idx, ci, repi] += 1 / (twin/1000)
            end
        end
    end

    meanrate = mean(rate, dims=3)[:,:]

    timev = twin*collect(0:nstep-1)
    tidx = (timev .>= stim_off).*(timev .< train_time)
    meanrate = meanrate[tidx,:]

    return meanrate

end
