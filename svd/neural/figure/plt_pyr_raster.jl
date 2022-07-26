function plt_pyr_raster(dirfig, dirsim, cells_R, units_R, spikeTime_R, eventTime_R, rate_R, pcor_R,
    cells_L, units_L, spikeTime_L, eventTime_L, rate_L, pcor_L, almOrd, matchedCells)

#---------- load rate_navg ----------#
Ncells = 5000
Ne = 2500
rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,1:Ne,:]
rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[:,1:Ne,:]

#---------- remove neurons with rate < 1Hz ----------#
meanrate_R = mean(rate_R, dims=2)[:]
meanrate_L = mean(rate_L, dims=2)[:]
idx1Hz_R = meanrate_R .> 1.0
idx1Hz_L = meanrate_L .> 1.0
idx1Hz = idx1Hz_R .* idx1Hz_L

rate_R = rate_R[idx1Hz,:]
rate_L = rate_L[idx1Hz,:]
cells_R = cells_R[idx1Hz]

#---------- sort neurons by pcor ----------#
pcor = (pcor_R + pcor_L)/2
pcorSorted = reverse(sortperm(pcor))

rate_R = rate_R[almOrd[pcorSorted],:]
rate_L = rate_L[almOrd[pcorSorted],:]
cells_sorted = cells_R[almOrd[pcorSorted]]


start_time = -2.0

for celli = 1:300
# for celli = 1:length(pcorSorted)
    println(celli)

    ci = cells_sorted[celli]
    nid = pcorSorted[celli]
    ci_model = Int(matchedCells[nid])

    ind_R = findall(x->x==ci, units_R)
    ntrial_ci_R = size(ind_R)[1]
    spikeTime_ci_R = spikeTime_R[ind_R]
    eventTime_ci_R = eventTime_R[ind_R]

    ind_L = findall(x->x==ci, units_L)
    ntrial_ci_L = size(ind_L)[1]
    spikeTime_ci_L = spikeTime_L[ind_L]
    eventTime_ci_L = eventTime_L[ind_L]

    # maximum of model rate (rate_R_navg) and actual rate (rate_R)
    ratemax = maximum([rate_R_navg400[:,ci_model,1][:]; rate_R[celli,:][:]; rate_L_navg400[:,ci_model,1][:]; rate_L[celli,:][:]])

    fig = figure(figsize=(0.8,1.3))
    gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0.001, left=0.3, bottom=0.2)

    ax = fig.add_subplot(gs[1,1])
    for triali = 1:ntrial_ci_R
        spktime = spikeTime_ci_R[triali] .- eventTime_ci_R[triali] .- start_time
        idx = (spktime .> 0).*(spktime .< 2.0)
        if sum(idx) > 0
            spktime = spktime[idx]
            plot(spktime, (triali + ntrial_ci_L)*ones(length(spktime)), color="blue", marker="o", ms=0.5, mec="None", linestyle="")        
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
    plot(timev, rate_R[celli,:][1:100], color="blue", linewidth=0.8, alpha=1.0)    
    plot(timev, rate_L[celli,:][1:100], color="red", linewidth=0.8, alpha=1.0)    
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
    #         plot(spktime, triali*ones(length(spktime)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
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
    # plot(timev, rate_L[celli,:][1:100], color="red", linewidth=0.8, alpha=1.0)    
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([0,2000])
    # ylim([0, ratemax])
    # xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
    # yticks(fontsize=8)
    # xlabel("time (s)", fontsize=8)
    # ylabel("spk / s ", fontsize=8)

    # gs.update(left=0.22, bottom=0.2)

    savefig(dirfig * "raster/pyr/pyr_$(celli).png", dpi=600)
    savefig(dirfig * "raster/pyr/pyr_$(celli).pdf", dpi=600)

    close()

end


end






# fig=figure(figsize=(1.5,1.5))
# gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.25, hspace=0.001, left=0.3, bottom=0.3)

# ax = fig.add_subplot(gs[1,0])
# for triali = 1:ntrial_ci_R
#     spktime = spikeTime_ci_R[triali] .- eventTime_ci_R[triali] .- start_time
#     idx = (spktime .> 0).*(spktime .< 2.0)
#     if sum(idx) > 0
#         spktime = spktime[idx]
#         plot(spktime, triali*ones(length(spktime)), color="blue", marker="o", ms=0.5, mec="None", linestyle="")        
#     end
# end
# xlim([0, 2])
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# ax.spines["left"].set_visible(false)
# ax.spines["bottom"].set_visible(false)
# tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
# title("R", fontsize=8)


# ax = fig.add_subplot(gs[0,0])
# timev = 20*collect(1:100)
# plot(timev, rate_R[celli,:][1:100], color="blue", linewidth=0.8, alpha=1.0)    
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# xlim([0,2000])
# ylim([0, ratemax])
# xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
# tick_params(labelleft=false)
# ax.set_xlabel("time (s)", fontsize=8)


# ax = fig.add_subplot(gs[1,1])
# for triali = 1:ntrial_ci_L
#     spktime = spikeTime_ci_L[triali] .- eventTime_ci_L[triali] .- start_time
#     idx = (spktime .> 0).*(spktime .< 2.0)
#     if sum(idx) > 0
#         spktime = spktime[idx]
#         plot(spktime, triali*ones(length(spktime)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
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
# plot(timev, rate_L[celli,:][1:100], color="red", linewidth=0.8, alpha=1.0)    
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# xlim([0,2000])
# ylim([0, ratemax])
# xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
# yticks(fontsize=8)
# xlabel("time (s)", fontsize=8)
# ylabel("spk / s ", fontsize=8)
