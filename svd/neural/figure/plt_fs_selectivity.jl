function plt_fs_selectivity(dirfig, rate_diff_fs, rate_R_fs, rate_L_fs)
        
    rate_mean_fs = mean((rate_R_fs + rate_L_fs)/2, dims=1)[:]
    rate_diffnorm_fs = rate_diff_fs ./ rate_mean_fs
    sortcells = sortperm(rate_diffnorm_fs)


    figure(figsize=(2.0,1.5))
    imshow((transpose(rate_R_fs - rate_L_fs)./rate_mean_fs)[reverse(sortcells),:], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # xlim([0,2])
    # ylim([0,length(cellsOrd_fs)])
    xticks([0,50,100], [-2,-1,0], fontsize=7)
    yticks([0,200], fontsize=7)
    # yticks([100, 200, 300], fontsize=8)
    xlabel("time (s)", fontsize=7)
    ylabel("FS neuron", fontsize=7)
    # ax.set_ylabel("FS cell", fontsize=8)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.set_yticklabels(["-1", "0", "1"])
    cbar.ax.tick_params(labelsize=7)
    tight_layout()

    # savefig("fig_selectivity_fs.png", dpi=600)
    savefig(dirfig * "fig_selectivity_fs.pdf", dpi=600)



    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.1, bottom=0.3)

    # ax = fig.add_subplot(gs[1])
    # for ci = 1:length(cellsOrd_fs)
    #     plot(times_L_fs[ci,1:ns_L_fs[ci]], ci*ones(ns_L_fs[ci]), color="red", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_fs)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # ax.set_ylabel("fs cell", fontsize=8)
    # title("lick left", fontsize=8)

    # ax = fig.add_subplot(gs[2])
    # for ci = 1:length(cellsOrd_fs)
    #     plot(times_R_fs[ci,1:ns_R_fs[ci]], ci*ones(ns_R_fs[ci]), color="blue", marker="o", ms=0.4, mec="None", linestyle="")
    # end
    # xlim([0,2])
    # ylim([0,length(cellsOrd_fs)])
    # xticks([0,1,2], [-2,-1,0], fontsize=8)
    # yticks([])
    # ax.set_xlabel("time (s)", fontsize=8)
    # title("lick right", fontsize=8)

    # ax = fig.add_subplot(gs[3])
    # rate_diff_pyr1trial = (ns_R_fs - ns_L_fs) / 2
    # plot(rate_diff_pyr1trial, collect(1:length(cellsOrd_fs)), color="gray", marker="o", ms=1.2, mec="None", linestyle="", alpha=1)
    # plot(rate_diff_fs[sortcells], collect(1:length(cellsOrd_fs)), color="magenta", lw=0.5, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # # ylim([0,length(cellsOrd)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # title("selectivity", fontsize=8)

    # savefig(dirfig * "fs_RL_raster.png", dpi=600)
    # savefig(dirfig * "fs_RL_raster.pdf", dpi=600)





    # fig = figure(figsize=(2.25,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=2, figure=fig, wspace=0.2, hspace=0.5, left=0.18, bottom=0.3, right=0.9, width_ratios=[1,2])

    # ax = fig.add_subplot(gs[1])
    # plot(rate_diffnorm_fs[sortcells], collect(1:length(rate_diffnorm_fs)), color="black", lw=0.8, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-2.5,2.5])
    # # ylim([0,length(cellsOrd)])
    # xticks([-2,0,2], ["-2", "0","2"], fontsize=7)
    # yticks([0,200], fontsize=7)
    # ax.set_xlabel("selectivity", fontsize=7)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # # title("selectivity", fontsize=8)
        

    # ax = fig.add_subplot(gs[2])
    # imshow((transpose(rate_L_fs - rate_R_fs)./rate_mean_fs)[sortcells,:], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # # xlim([0,2])
    # # ylim([0,length(cellsOrd_fs)])
    # xticks([0,50,100], [-2,-1,0], fontsize=7)
    # yticks([])
    # # yticks([100, 200, 300], fontsize=8)
    # ax.set_xlabel("time (s)", fontsize=7)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # cbar = colorbar(ticks=[-1,0,1])
    # cbar.ax.set_yticklabels(["-1", "0", "1"])
    # cbar.ax.tick_params(labelsize=7)

    # savefig("fig_selectivity_fs.png", dpi=600)
    # # savefig(dirfig * "fs_RL_psth.pdf", dpi=600)



#-----------------------------------------------------

    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.18, bottom=0.3, right=0.9, width_ratios=[1,1,2])

    # ax = fig.add_subplot(gs[1])
    # rate_diff_pyr1trial = (ns_R_fs - ns_L_fs) / 2
    # plot(rate_diff_pyr1trial, collect(1:length(cellsOrd_fs)), color="gray", marker="o", ms=1.2, mec="None", linestyle="", alpha=1)
    # # plot(rate_diff_fs[sortcells], collect(1:length(cellsOrd_fs)), color="magenta", lw=0.5, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-30,30])
    # ylim([0,length(cellsOrd_fs)])
    # xticks([-25,0,25], ["-25", "0","25"], fontsize=8)
    # yticks([0,200], fontsize=8)
    # # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # ax.set_ylabel("FS cell", fontsize=8)
    # # title("selectivity", fontsize=8)


    # ax = fig.add_subplot(gs[2])
    # # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff_fs[sortcells], collect(1:length(cellsOrd_fs)), color="black", lw=0.8, linestyle="-", alpha=1)
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
    # imshow((rate_R_fs - rate_L_fs)[reverse(sortcells),:], vmin=-25, vmax=25, cmap="bwr_r", aspect="auto")
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

    # savefig(dirfig * "fs_RL_psth.png", dpi=600)
    # savefig(dirfig * "fs_RL_psth.pdf", dpi=600)


end
