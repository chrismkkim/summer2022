function plt_fs_raster(dirfig, dirsim, cells_R, units_R, spikeTime_R, eventTime_R, rate_R, pcor_R, mse_R,
    cells_L, units_L, spikeTime_L, eventTime_L, rate_L, pcor_L, mse_L, cells_matched)

# load rate_navg
Ncells = 5000
Ne = 2500
rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,Ne+1:Ncells,:]
rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[:,Ne+1:Ncells,:]

start_time = -2.0

# pcor = (pcor_R + pcor_L)/2
# pcorSorted = reverse(sortperm(pcor))
# cells_sorted = cells_R[pcorSorted]

mse = (mse_R + mse_L)/2
mseSorted = reverse(sortperm(mse))
cells_sorted = cells_R[mseSorted]

for celli = 1:length(mseSorted)
# for celli = 1:10
    println(celli)

    nid = mseSorted[celli]
    ci = cells_sorted[celli]
    ci_model = Int(cells_matched[nid])

    ind_R = findall(x->x==ci, units_R)
    ntrial_ci_R = size(ind_R)[1]
    spikeTime_ci_R = spikeTime_R[ind_R]
    eventTime_ci_R = eventTime_R[ind_R]

    ind_L = findall(x->x==ci, units_L)
    ntrial_ci_L = size(ind_L)[1]
    spikeTime_ci_L = spikeTime_L[ind_L]
    eventTime_ci_L = eventTime_L[ind_L]

    # maximum of model rate (rate_R_navg) and actual rate (rate_R)
    ratemax = maximum([rate_R_navg400[:,ci_model,1][:]; rate_R[nid,:][:]; rate_L_navg400[:,ci_model,1][:]; rate_L[nid,:][:]])

    fig=figure(figsize=(0.8,1.3))
    gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0.001, left=0.3, bottom=0.2)

    ax = fig.add_subplot(gs[1,1])
    for triali = 1:ntrial_ci_R
        spktime = spikeTime_ci_R[triali] .- eventTime_ci_R[triali] .- start_time
        idx = (spktime .> 0).*(spktime .< 2.0)
        if sum(idx) > 0
            spktime = spktime[idx]
            plot(spktime, (triali+ntrial_ci_L)*ones(length(spktime)), color="blue", marker="o", ms=0.5, mec="None", linestyle="")        
        end
    end
    for triali = 1:ntrial_ci_L
        spktime = spikeTime_ci_L[triali] .- eventTime_ci_L[triali] .- start_time
        idx = (spktime .> 0).*(spktime .< 2.0)
        if sum(idx) > 0
            spktime = spktime[idx]
            plot(spktime, triali*ones(length(spktime)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
        end
    end
    xlim([0, 2])
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.spines["left"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    # title("R", fontsize=8)


    ax = fig.add_subplot(gs[2,1])
    timev = 20*collect(1:100)
    plot(timev, rate_R[nid,:][1:100], color="blue", linewidth=0.8, alpha=1.0)    
    plot(timev, rate_L[nid,:][1:100], color="red", linewidth=0.8, alpha=1.0)    
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xlim([0,2000])
    ylim([0, ratemax])
    xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
    yticks(fontsize=8)

    # tick_params(labelleft=false)
    # ax.set_xlabel("time (s)", fontsize=8)


    # ax = fig.add_subplot(gs[1,1])
    # for triali = 1:ntrial_ci_L
    #     spktime = spikeTime_ci_L[triali] .- eventTime_ci_L[triali] .- start_time
    #     idx = (spktime .> 0).*(spktime .< 2.0)
    #     if sum(idx) > 0
    #         spktime = spktime[idx]
    #         plot(spktime, triali*ones(length(spktime)), color="red", marker="o", ms=1.0, mec="None", linestyle="")        
    #     end
    # end
    # xlim([0, 2])
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # ax.spines["left"].set_visible(false)
    # ax.spines["bottom"].set_visible(false)
    # tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    # title("L", fontsize=8)


    # ax = fig.add_subplot(gs[0,1])
    # timev = 20*collect(1:100)
    # plot(timev, rate_L[nid,:][1:100], color="red", linewidth=0.8, alpha=1.0)    
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([0,2000])
    # ylim([0, ratemax])
    # xticks([0,1000,2000],[0,1,2],fontsize=8)
    # yticks(fontsize=8)
    # xlabel("time (s)", fontsize=8)
    # ylabel("spk / s ", fontsize=8)

    # gs.update(left=0.22, bottom=0.2)

    savefig(dirfig * "raster/fs/fs_$(celli).png", dpi=600)
    savefig(dirfig * "raster/fs/fs_$(celli).pdf", dpi=600)

    close()

end


    # Ne = 2500
    # Ncells = 5000

    # # load rate_navg
    # rate_navg50 = load(dirsim * "rate_R_navg50.jld", "rate")[:,Ne+1:Ncells,:]
    # rate_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,Ne+1:Ncells,:]
    
    # # performance: Pyr - Exc
    # pcorSorted = reverse(sortperm(pcorMatched))
    # pcor_high = pcorSorted[pcorMatched[pcorSorted] .> 0.7]
    # pcor_med = pcorSorted[pcorMatched[pcorSorted] .> 0.5]
    # pcor_low = pcorSorted[pcorMatched[pcorSorted] .> 0.2]

    # # load spikes
    # navg = 50
    # Ncells = 5000
    # tmp = load(dirsim * "times_R_1.jld", "times")
    # times = zeros(size(tmp)[1], size(tmp)[2], navg)
    # ns = zeros(Int,Ncells, navg)
    # for avgi = 1:navg
    #     times[:,:,avgi] = load(dirsim * "times_R_$(avgi).jld", "times")
    #     ns[:,avgi] = load(dirsim * "ns_R_$(avgi).jld", "ns")
    # end
    # times = times[Ne+1:Ncells,:,:]
    # ns = ns[Ne+1:Ncells,:]

    # #--------------------------------------------#
    # #----- spike raster of selected neurons -----#
    # #--------------------------------------------#
    # for kk = 1:1
    #     nid = pcor_high[kk]
    #     # nid = pcor_med[end-kk+1]
    #     # nid = pcor_low[end-kk+1]

    #     ci = Int(cellsMatched[nid])
    #     ci_alm = nid
    #     fig = figure(figsize=(1.5,2.0))
    #     gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0)

    #     ax = fig.add_subplot(gs[1])
    #     for ii = 1:navg
    #         tidx = times[ci,:,ii] .> 1000
    #         plot(times[ci,tidx,ii], ii*ones(sum(tidx)), color="black", marker="o", ms=0.5, linestyle="")        
    #         xlim([1000,3000])
    #     end
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     ax.spines["left"].set_visible(false)
    #     ax.spines["bottom"].set_visible(false)
    #     tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    #     # fig.subplots_adjust(bottom = 0.001)

    #     ax = fig.add_subplot(gs[0])
    #     timev = 20*collect(1:100)
    #     plot(timev, rtargFS[:,ci_alm], color="red", linewidth=0.8, alpha=1.0)    
    #     plot(timev, rate_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)
    #     plot(timev, rate_navg400[:,ci,1], color="black", linewidth=0.8, alpha=1.0)    
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     xlim([0,2000])
    #     xticks([0,1000,2000],[0,1,2],fontsize=8)
    #     yticks(fontsize=8)
    #     xlabel("time (s)", fontsize=8)
    #     ylabel("spk / s ", fontsize=8)

    #     tight_layout()

    #     savefig(dirfig * "inh$(kk)_raster_pcor$(round(pcorMatched[nid], digits=2))_$(ci).png", dpi=300)
    #     close()
    # end





end