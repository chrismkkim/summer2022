function plt_inh_raster(dirfig, dirsim, rtargFS, pcorMatched, cellsMatched)

    Ne = 2500
    Ncells = 5000

    # load rate_navg
    rate_navg50 = load(dirsim * "rate_R_navg50.jld", "rate")[:,Ne+1:Ncells,:]
    rate_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,Ne+1:Ncells,:]
    
    # performance: Pyr - Exc
    pcorSorted = reverse(sortperm(pcorMatched))
    pcor_high = pcorSorted[pcorMatched[pcorSorted] .> 0.7]
    pcor_med = pcorSorted[pcorMatched[pcorSorted] .> 0.5]
    pcor_low = pcorSorted[pcorMatched[pcorSorted] .> 0.2]

    # load spikes
    navg = 50
    Ncells = 5000
    tmp = load(dirsim * "times_R_1.jld", "times")
    times = zeros(size(tmp)[1], size(tmp)[2], navg)
    ns = zeros(Int,Ncells, navg)
    for avgi = 1:navg
        times[:,:,avgi] = load(dirsim * "times_R_$(avgi).jld", "times")
        ns[:,avgi] = load(dirsim * "ns_R_$(avgi).jld", "ns")
    end
    times = times[Ne+1:Ncells,:,:]
    ns = ns[Ne+1:Ncells,:]

    #--------------------------------------------#
    #----- spike raster of selected neurons -----#
    #--------------------------------------------#
    for kk = 1:1
        nid = pcor_high[kk]
        # nid = pcor_med[end-kk+1]
        # nid = pcor_low[end-kk+1]

        ci = Int(cellsMatched[nid])
        ci_alm = nid
        fig = figure(figsize=(1.5,2.0))
        gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0)

        ax = fig.add_subplot(gs[1])
        for ii = 1:navg
            tidx = times[ci,:,ii] .> 1000
            plot(times[ci,tidx,ii], ii*ones(sum(tidx)), color="black", marker="o", ms=0.5, linestyle="")        
            xlim([1000,3000])
        end
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["bottom"].set_visible(false)
        tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
        # fig.subplots_adjust(bottom = 0.001)

        ax = fig.add_subplot(gs[0])
        timev = 20*collect(1:100)
        plot(timev, rtargFS[:,ci_alm], color="red", linewidth=0.8, alpha=1.0)    
        plot(timev, rate_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)
        plot(timev, rate_navg400[:,ci,1], color="black", linewidth=0.8, alpha=1.0)    
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        xlim([0,2000])
        xticks([0,1000,2000],[0,1,2],fontsize=8)
        yticks(fontsize=8)
        xlabel("time (s)", fontsize=8)
        ylabel("spk / s ", fontsize=8)

        tight_layout()

        savefig(dirfig * "fs_raster_pcor$(round(pcorMatched[nid], digits=2))_$(ci).png", dpi=300)
        close()
    end





end