matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

clrGray = [cmap.Greys(x) for x in [0.4, 0.8]]
clrGreen = [cmap.Greens(x) for x in [0.4, 0.8]]

function plt_selectivity_pop(dirfig, pyr_R, pyr_L, fs_R, fs_L, excDiff_1trial, inhDiff_1trial, matchedCells)
    
    # load trial averaged spiking rates
    # rate_R : [time x neurons x batch_size]
    # rate_L : [time x neurons x batch_size]
    fname_rate_R = dirsim * "rate_R_navg400.jld"
    fname_rate_L = dirsim * "rate_L_navg400.jld"
    rate_R = load(fname_rate_R, "rate")
    rate_L = load(fname_rate_L, "rate")

    # rate difference in Pyr and FS cells
    pyrDiff = mean(pyr_L - pyr_R, dims=1)[:] # mean in time
    fsDiff = mean(fs_L - fs_R, dims=1)[:] # mean in time

    # rate difference in Exc and Inh neuron models    
    # --- trial averaged ---
    # matchedCells: Exc neurons matched with ALM Pyr cells
    excDiff = mean(rate_L[:,matchedCells,1] - rate_R[:,matchedCells,1], dims=1)[:] # mean in time
    inhDiff = mean(rate_L[:,2501:end,1] - rate_R[:,2501:end,1], dims=1)[:] # mean in time

    # --- single trial ---
    excDiff_1trial = excDiff_1trial[matchedCells]
    excDiff_1trial = excDiff_1trial[sortperm(excDiff)] # sort by trial averaged rate diff
    inhDiff_1trial = inhDiff_1trial[sortperm(inhDiff)] # sort by trial averaged rate diff

    # --- neuron averaged ---
    excDiff_1trial_avg = funMovAvg(excDiff_1trial, 10) # average 20 adjacent neurons
    inhDiff_1trial_avg = funMovAvg(inhDiff_1trial, 10) # average 20 adjacent neurons


    figure(figsize=(3.0,1.8))
    # # ax = subplot(211)    
    plot(excDiff_1trial, color="blue", marker=".", ms=1, linestyle="", alpha=0.5)            
    plot(excDiff_1trial_avg, color="cyan", linestyle="-", lw=0.8, alpha=1)
    plot(excDiff[sortperm(excDiff)], color="black", linestyle="--", lw=0.8)
    # plot(collect(1:length(excDiff)), zeros(length(excDiff)), color="magenta", linestyle="--")
    legend_elements = [Line2D([0], [0], color="blue", lw=1, label="one trial"),
                       Line2D([0], [0], color="black", lw=1, label="avg trial"),
                       Line2D([0], [0], color="cyan", lw=1, label="avg neuron")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
    # legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, loc=2)     
    xlabel("exc neuron", fontsize=8)
    ylabel(L"$\Delta rate$" * " (Hz)", fontsize=8)
    xticks([0,1000,2000], fontsize=8)
    yticks([-20,0,20], fontsize=8)
    ylim([-23, 23])
    tight_layout()

    savefig(dirfig * "selectivity_e.png", dpi=600)



    figure(figsize=(3.0,1.8))
    plot(inhDiff_1trial, color="red", marker=".", ms=1, linestyle="", alpha=0.5)        
    plot(inhDiff_1trial_avg, color="yellow", linestyle="-", lw=0.8, alpha=1)    
    plot(inhDiff[sortperm(inhDiff)], color="black", linestyle="--", lw=0.8)    
    legend_elements = [Line2D([0], [0], color="red", lw=1, label="one trial"),
                       Line2D([0], [0], color="black", lw=1, label="avg trial"),
                       Line2D([0], [0], color="gold", lw=1, label="avg neuron")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
    # legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, loc=2)     
    xlabel("inh neuron", fontsize=8)
    ylabel(L"$\Delta rate$" * " (Hz)", fontsize=8)
    xticks([0,1000,2000], fontsize=8)
    yticks([-20,0,20], fontsize=8)
    ylim([-23, 23])
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)

    tight_layout()

    savefig(dirfig * "selectivity_i.png", dpi=600)




    figure(figsize=(2.0,1.8))
    # ax = subplot(211)
    plot(pyrDiff[sortperm(pyrDiff)], color="limegreen", marker="o", ms=1, linestyle="")
    plot(collect(1:length(pyrDiff)), zeros(length(pyrDiff)), color="gray", lw=0.8, linestyle="--")
    # plot(fsDiff, color="darkorange", marker="o", ms=2, linestyle="")
    xlabel("Pyr cell", fontsize=8)
    ylabel(L"$\Delta rate (Hz)$", fontsize=8)
    xticks(fontsize=8)
    yticks([-20,0,20], fontsize=8)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    tight_layout()

    savefig(dirfig * "selectivity_pyr.png", dpi=600)



    figure(figsize=(2.0,1.8))
    # ax = subplot(212)
    # plot(pyrDiff, color="limegreen", marker="o", ms=2, linestyle="")
    plot(fsDiff[sortperm(fsDiff)], color="darkorange", marker="o", ms=1, linestyle="")
    plot(collect(1:length(fsDiff)), zeros(length(fsDiff)), color="gray", lw=0.8, linestyle="--")
    xlabel("FS cell", fontsize=8)
    ylabel(L"$\Delta rate (Hz)$", fontsize=8)
    xticks(fontsize=8)
    yticks([-20,0,20], fontsize=8)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    tight_layout()

    savefig(dirfig * "selectivity_fs.png", dpi=600)




    # for ii = 1:length(navg_ind)
    #     ax = subplot(2,1,ii)
    #     cnt_R, bins_R = hist(pcor_R[:,ii], bins=100, range=(-1,1), color="blue", histtype="step")
    #     cnt_L, bins_L = hist(pcor_L[:,ii], bins=100, range=(-1,1), color="red", histtype="step")
    #     pmax_R = maximum(cnt_R)
    #     pmax_L = maximum(cnt_L)
    #     ptmpR = pcor_R[:,ii]
    #     ptmpL = pcor_L[:,ii]
    #     medR = median(ptmpR[.!isnan.(ptmpR)])
    #     medL = median(ptmpL[.!isnan.(ptmpL)])
    #     plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", linestyle="--")
    #     plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", linestyle="--")
    #     title(L"$N_{trial}=$" * "$(navg_list[navg_ind[ii]])", fontsize=12)

    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     if ii == 1
    #         legend_elements = [Line2D([0], [0], color="blue", lw=2, label="right"),
    #                            Line2D([0], [0], color="red", lw=2, label="left")]
    #         legend(loc=2, handles=legend_elements, frameon=false, fontsize=10, handlelength=1)     
    #         tick_params(labelbottom=false)
    #         yticks(fontsize=12)
    #     else
    #         xticks(fontsize=12)
    #         yticks(fontsize=12)
    #     end

    # end
    # xlabel("correlation", fontsize=12)
    # ylabel("count", fontsize=12)
    # tight_layout()

    # savefig(dirfig * "pyr_pcor.png", dpi=300)
     
end
    