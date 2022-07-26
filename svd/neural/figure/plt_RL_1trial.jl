function plt_RL_1trial(dirfig, pcor_exc, pcor_inh, pcor_null, pcor_pyr, pcor_fs)
        

    ntrials = 20

    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    hist(pcor_pyr[1:ntrials], bins=25, range=(0,1), color="limegreen", lw=0.8, histtype="step")
    hist(pcor_fs[1:ntrials], bins=25, range=(0,1), color="darkorange", lw=0.8, histtype="step")
    legend_elements = [Line2D([0], [0], color="limegreen", lw=1.5, label="Pyr"),
                       Line2D([0], [0], color="darkorange", lw=1.5, label="FS")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xticks([0,0.5,1], fontsize=8)
    yticks(fontsize=8)
    xlabel("single trial selectivity", fontsize=8)
    ylabel("count", fontsize=8)

    tight_layout()

    savefig(dirfig * "selectivity_1trial_data.png", dpi=600)
    savefig(dirfig * "selectivity_1trial_data.pdf", dpi=600)


    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    hist(pcor_exc, bins=50, range=(0,1), color="limegreen", lw=0.8, histtype="step")
    hist(pcor_inh, bins=50, range=(0,1), color="darkorange", lw=0.8, histtype="step")
    hist(pcor_null, bins=50, range=(0,1), color="gray", alpha=0.5, density=true)
    legend_elements = [Line2D([0], [0], color="limegreen", lw=1.5, label="Exc"),
                       Line2D([0], [0], color="darkorange", lw=1.5, label="Inh"),
                       Line2D([0], [0], color="gray", lw=1.5, label="Null")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xticks([0,0.5,1], fontsize=8)
    yticks(fontsize=8)
    xlabel("single trial selectivity", fontsize=8)
    ylabel("count", fontsize=8)

    tight_layout()

    savefig(dirfig * "selectivity_1trial_model.png", dpi=600)
    savefig(dirfig * "selectivity_1trial_model.pdf", dpi=600)    



    # figure(figsize=(1.7,1.5))
    # ax = subplot(111)
    # ii = 2 # navg = 400 (number of trial averages)
    # cnt_R, bins_R = hist(pcor_R[:,ii], bins=50, range=(-1,1), lw=0.8, color="blue", histtype="step")
    # cnt_L, bins_L = hist(pcor_L[:,ii], bins=50, range=(-1,1), lw=0.8, color="red", histtype="step")
    # pmax_R = maximum(cnt_R)
    # pmax_L = maximum(cnt_L)
    # ptmpR = pcor_R[:,ii]
    # ptmpL = pcor_L[:,ii]
    # medR = median(ptmpR[.!isnan.(ptmpR)])
    # medL = median(ptmpL[.!isnan.(ptmpL)])
    # plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", lw=0.8, linestyle="--")
    # plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", lw=0.8, linestyle="--")
    # # title(L"$N_{trial}=$" * "$(navg_list[navg_ind[ii]])", fontsize=12)

    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="right"),
    #                     Line2D([0], [0], color="red", lw=1.5, label="left")]
    # legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     
    # xticks([-1,-0.5,0,0.5,1], ["-1","-0.5","0","0.5","1"],fontsize=8)
    # yticks(fontsize=8)
    
    # xlabel("performance", fontsize=8)
    # ylabel("count", fontsize=8)
    # tight_layout()

    # savefig(dirfig * "pyr_pcor.png", dpi=600)
    # savefig(dirfig * "pyr_pcor.pdf", dpi=600)


    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.1, bottom=0.3)
    
    # ax = fig.add_subplot(gs[1])
    # for ci = 1:length(cellsOrd)
    #     nid = cellsOrd[ci]
    #     plot(times_L[nid,1:ns_L[nid]], ci*ones(ns_L[nid]), color="red", marker="o", ms=0.1, mec="None", linestyle="")
    # end
    # xlim([1000,p.train_time])
    # ylim([0,length(cellsOrd)])
    # xticks([1000, 2000, 3000],[-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # ax.set_ylabel("inh neuron", fontsize=8)
    # title("lick left", fontsize=8)
    
    # ax = fig.add_subplot(gs[2])
    # for ci = 1:length(cellsOrd)
    #     nid = cellsOrd[ci]
    #     plot(times_R[nid,1:ns_R[nid]], ci*ones(ns_R[nid]), color="blue", marker="o", ms=0.1, mec="None", linestyle="")
    # end
    # xlim([1000,p.train_time])
    # ylim([0,length(cellsOrd)])
    # # xticks([0, 1000, 2000, 3000],[0, 1,2,3], fontsize=8)
    # xticks([1000, 2000, 3000],[-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # title("lick right", fontsize=8)
    
    # ax = fig.add_subplot(gs[3])
    # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff_1trial[cellsOrd], collect(1:length(cellsOrd)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # plot(rate_diff[cellsOrd], collect(1:length(cellsOrd)), color="magenta", lw=0.8, linestyle="--", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # title("selectivity", fontsize=8)
    
    # savefig(dirfig * "inh_RL_raster.png", dpi=600)
    # savefig(dirfig * "inh_RL_raster.pdf", dpi=600)
    
    

    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=2, figure=fig, wspace=0.2, hspace=0.5, left=0.15, bottom=0.3, right=0.9, width_ratios=[1,2])

    # ax = fig.add_subplot(gs[1,1])
    # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff[cellsOrd], collect(1:length(cellsOrd)), color="black", lw=0.8, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-15,15])
    # # ylim([0,length(cellsOrd)])
    # xticks([-15,0,15], ["-15", "0","15"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # ax.set_ylabel("Inh neuron", fontsize=8)
    # # title("selectivity", fontsize=8)
        
    # ax = fig.add_subplot(gs[2])
    # imshow(transpose(Inhrate_R - Inhrate_L)[reverse(cellsOrd),:], vmin=-15, vmax=15, cmap="bwr_r", aspect="auto")
    # # xlim([0,2])
    # # ylim([0,length(cellsOrd_fs)])
    # xticks([0,50,100], [-2,-1,0], fontsize=8)
    # yticks([])
    # # yticks([100, 200, 300], fontsize=8)
    # ax.set_xlabel("time (s)", fontsize=8)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # cbar = colorbar(ticks=[-15,0,15])
    # cbar.ax.set_yticklabels(["-15", "0", "15"])
    # cbar.ax.tick_params(labelsize=8)

    # savefig(dirfig * "inh_RL_psth.png", dpi=600)
    # savefig(dirfig * "inh_RL_psth.pdf", dpi=600)

end
