function plt_pyr_selectivity(dirfig, rate_diff_pyr, rate_R_pyr, rate_L_pyr)
        

    rate_mean_pyr = mean((rate_R_pyr + rate_L_pyr)/2, dims=1)[:]
    rate_diffnorm_pyr = rate_diff_pyr ./ rate_mean_pyr
    sortcells = sortperm(rate_diffnorm_pyr)

    figure(figsize=(2.0,1.5))
    imshow((transpose(rate_R_pyr - rate_L_pyr)./rate_mean_pyr)[reverse(sortcells),:], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # xlim([0,2])
    # ylim([0,length(cellsOrd_fs)])
    xticks([0,50,100], [-2,-1,0], fontsize=7)
    yticks([0,1000], fontsize=7)
    # yticks([100, 200, 300], fontsize=8)
    xlabel("time (s)", fontsize=7)
    ylabel("Pyr neuron", fontsize=7)
    # ax.set_ylabel("FS cell", fontsize=8)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.set_yticklabels(["-1", "0", "1"])
    cbar.ax.tick_params(labelsize=7)
    tight_layout()

    # savefig("fig_selectivity_pyr.png", dpi=600)
    savefig(dirfig * "fig_selectivity_pyr.pdf", dpi=600)




    # sortcells = sortperm(rate_diff_pyr)

    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.15, bottom=0.3, right=0.9, width_ratios=[1,1,2])

    # ax = fig.add_subplot(gs[1])
    # rate_diff_pyr1trial = (ns_R_pyr - ns_L_pyr) / 2
    # plot(rate_diff_pyr1trial, collect(1:length(cellsOrd_pyr)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # # plot(rate_diff_fs[sortcells], collect(1:length(cellsOrd_fs)), color="magenta", lw=0.5, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # ax.set_ylabel("Pyr cell", fontsize=8)
    # # title("selectivity", fontsize=8)


    # ax = fig.add_subplot(gs[2])
    # # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff_pyr[sortcells], collect(1:length(cellsOrd_pyr)), color="black", lw=0.8, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # # title("selectivity", fontsize=8)
        

    # ax = fig.add_subplot(gs[3])
    # imshow(transpose(rate_R_pyr - rate_L_pyr)[reverse(sortcells),:], vmin=-25, vmax=25, cmap="bwr_r", aspect="auto")
    # # xlim([0,2])
    # # ylim([0,length(cellsOrd_fs)])
    # xticks([0,50,100], [-2,-1,0], fontsize=8)
    # yticks([])
    # # yticks([100, 200, 300], fontsize=8)
    # ax.set_xlabel("time (s)", fontsize=8)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # cbar = colorbar(ticks=[-25,0,25])
    # cbar.ax.set_yticklabels(["-25", "0", "25"])
    # cbar.ax.tick_params(labelsize=8)


    # savefig(dirfig * "pyr_RL_psth.png", dpi=600)
    # savefig(dirfig * "pyr_RL_psth.pdf", dpi=600)






    # ax = fig.add_subplot(gs[1])
    # for ci = 1:length(cellsOrd_pyr)
    #     plot(times_L_pyr[ci,1:ns_L_pyr[ci]], ci*ones(ns_L_pyr[ci]), color="red", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_pyr)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # ax.set_ylabel("pyr cell", fontsize=8)
    # title("lick left", fontsize=8)

    # ax = fig.add_subplot(gs[2])
    # for ci = 1:length(cellsOrd_pyr)
    #     plot(times_R_pyr[ci,1:ns_R_pyr[ci]], ci*ones(ns_R_pyr[ci]), color="blue", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_pyr)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # title("lick right", fontsize=8)

    # ax = fig.add_subplot(gs[3])
    # sortcells = sortperm(rate_diff_pyr)
    # rate_diff_pyr1trial = (ns_R_pyr - ns_L_pyr) / 2
    # plot(rate_diff_pyr1trial, collect(1:length(cellsOrd_pyr)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # plot(rate_diff_pyr[sortcells], collect(1:length(cellsOrd_pyr)), color="magenta", lw=0.5, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # title("selectivity", fontsize=8)

    # savefig(dirfig * "pyr_RL_raster.png", dpi=600)
    # savefig(dirfig * "pyr_RL_raster.pdf", dpi=600)




    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.1, bottom=0.3)

    # ax = fig.add_subplot(gs[1])
    # for ci = 1:length(cellsOrd_pyr)
    #     plot(times_L_pyr[ci,1:ns_L_pyr[ci]], ci*ones(ns_L_pyr[ci]), color="red", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_pyr)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # ax.set_ylabel("pyr cell", fontsize=8)
    # title("lick left", fontsize=8)

    # ax = fig.add_subplot(gs[2])
    # for ci = 1:length(cellsOrd_pyr)
    #     plot(times_R_pyr[ci,1:ns_R_pyr[ci]], ci*ones(ns_R_pyr[ci]), color="blue", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_pyr)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # title("lick right", fontsize=8)

    # ax = fig.add_subplot(gs[3])
    # sortcells = sortperm(rate_diff_pyr)
    # rate_diff_pyr1trial = (ns_R_pyr - ns_L_pyr) / 2
    # plot(rate_diff_pyr1trial, collect(1:length(cellsOrd_pyr)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # plot(rate_diff_pyr[sortcells], collect(1:length(cellsOrd_pyr)), color="magenta", lw=0.5, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # title("selectivity", fontsize=8)

    # savefig(dirfig * "pyr_RL_raster.png", dpi=600)
    # savefig(dirfig * "pyr_RL_raster.pdf", dpi=600)

end
