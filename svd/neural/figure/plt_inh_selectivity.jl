function plt_inh_selectivity(dirfig, rate_diff, Inhrate_R, Inhrate_L)
        
    # inhibitory neurons only
    Ncells = 5000
    Ne = 2500
    # times_R = times_R[Ne+1:Ncells,:]
    # times_L = times_L[Ne+1:Ncells,:]
    # ns_R = ns_R[Ne+1:Ncells]
    # ns_L = ns_L[Ne+1:Ncells]

    # order inhibitory neurons by mean rate difference
    rate_mean = mean((Inhrate_R + Inhrate_L)/2, dims=1)[:]
    rate_diffnorm = rate_diff ./ rate_mean
    cellsOrd = sortperm(rate_diffnorm)


    figure(figsize=(2.0,1.5))
    startind = 3
    imshow((transpose(Inhrate_R - Inhrate_L)./rate_mean)[reverse(cellsOrd),startind:end], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # xlim([0,2])
    # ylim([0,length(cellsOrd_fs)])
    xticks([0,50,100-startind], [-2,-1,0], fontsize=7)
    yticks([0,2000], fontsize=7)
    # yticks([100, 200, 300], fontsize=8)
    xlabel("time (s)", fontsize=7)
    ylabel("Inh neuron", fontsize=7)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.set_yticklabels(["-1", "0", "1"])
    cbar.ax.tick_params(labelsize=7)
    tight_layout()

    # savefig("fig_selectivity_inh.png", dpi=600)
    savefig(dirfig * "fig_selectivity_inh.pdf", dpi=600)



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
    




    # fig = figure(figsize=(2.3,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=2, figure=fig, wspace=0.2, hspace=0.5, left=0.18, bottom=0.3, right=0.9, width_ratios=[1,2])

    # ax = fig.add_subplot(gs[1])
    # # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diffnorm[cellsOrd], collect(1:length(cellsOrd)), color="black", lw=0.8, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-2.5,2.5])
    # # ylim([0,length(cellsOrd)])
    # xticks([-2,0,2], ["-2", "0","2"], fontsize=7)
    # yticks([0,2000], fontsize=7)
    # ax.set_xlabel("selectivity", fontsize=7)
    # # ax.set_ylabel("Inh neuron", fontsize=8)
    # # title("selectivity", fontsize=8)
    # # title("averaged", fontsize=8)
        
    # ax = fig.add_subplot(gs[2])
    # startind = 3
    # imshow((transpose(Inhrate_L - Inhrate_R)./rate_mean)[cellsOrd,startind:end], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # # xlim([0,2])
    # # ylim([0,length(cellsOrd_fs)])
    # xticks([0,50,100-startind], [-2,-1,0], fontsize=7)
    # yticks([])
    # # yticks([100, 200, 300], fontsize=8)
    # ax.set_xlabel("time (s)", fontsize=7)
    # # ax.set_ylabel("FS cell", fontsize=8)
    # cbar = colorbar(ticks=[-1,0,1])
    # cbar.ax.set_yticklabels(["-1", "0", "1"])
    # cbar.ax.tick_params(labelsize=7)

    # savefig("fig_selectivity_inh.png", dpi=600)
    # # savefig(dirfig * "inh_RL_psth.pdf", dpi=600)



#-----------------------------------------------------    

    # fig = figure(figsize=(3.05,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.18, bottom=0.3, right=0.9, width_ratios=[1,1,2])

    # ax = fig.add_subplot(gs[1])
    # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff_1trial[cellsOrd], collect(1:length(cellsOrd)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # # plot(rate_diff[cellsOrd], collect(1:length(cellsOrd)), color="magenta", lw=0.8, linestyle="--", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-20,20])
    # # ylim([0,length(cellsOrd)])
    # xticks([-15,0,15], ["-15", "0","15"], fontsize=8)
    # yticks([0,2000], fontsize=8)
    # # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # ax.set_ylabel("Inh neuron", fontsize=8)
    # # title("single", fontsize=8)


    # ax = fig.add_subplot(gs[2])
    # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff[cellsOrd], collect(1:length(cellsOrd)), color="black", lw=0.8, linestyle="-", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-20,20])
    # # ylim([0,length(cellsOrd)])
    # xticks([-15,0,15], ["-15", "0","15"], fontsize=8)
    # yticks([])
    # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # # ax.set_ylabel("Inh neuron", fontsize=8)
    # # title("selectivity", fontsize=8)
    # # title("averaged", fontsize=8)
        
    # ax = fig.add_subplot(gs[3])
    # startind = 3
    # imshow(transpose(Inhrate_R - Inhrate_L)[reverse(cellsOrd),startind:end], vmin=-15, vmax=15, cmap="bwr_r", aspect="auto")
    # # xlim([0,2])
    # # ylim([0,length(cellsOrd_fs)])
    # xticks([0,50,100-startind], [-2,-1,0], fontsize=8)
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
