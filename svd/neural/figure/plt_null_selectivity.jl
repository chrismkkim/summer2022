function plt_null_selectivity(dirfig, rate_diff_null, rate_R_null, rate_L_null)
        
    rate_mean_null = mean((rate_R_null + rate_L_null)/2, dims=1)[:]
    rate_diffnorm_null = rate_diff_null ./ rate_mean_null
    sortcells = sortperm(rate_diffnorm_null)


    figure(figsize=(2.0,1.5))
    imshow((transpose(rate_R_null - rate_L_null)./rate_mean_null)[reverse(sortcells),:], vmin=-1, vmax=1, cmap="bwr_r", aspect="auto")
    # xlim([0,2])
    # ylim([0,length(cellsOrd_null)])
    xticks([0,50,100], [-2,-1,0], fontsize=7)
    yticks([0,1000, 2000], fontsize=7)
    # yticks([100, 200, 300], fontsize=8)
    xlabel("time (s)", fontsize=7)
    ylabel("Neuron", fontsize=7)
    # ax.set_ylabel("null cell", fontsize=8)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.set_yticklabels(["-1", "0", "1"])
    cbar.ax.tick_params(labelsize=7)
    tight_layout()

    # savefig("fig_selectivity_null.png", dpi=600)
    savefig(dirfig * "fig_selectivity_null.pdf", dpi=600)



    figure(figsize=(1.5,1.4))
    ax = subplot(111)
    hist(rate_diffnorm_null, bins=100, range=(-2,2), lw=0.8, color="gray", histtype="step", density=true)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="gray", lw=1.5, label="Null")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=7, handlelength=0.8)     
    xticks([-2,-1,0,1,2], ["-2","-1","0","1","2"],fontsize=7)
    yticks(fontsize=7)
    
    xlabel("selectivity", fontsize=7)
    ylabel("density", fontsize=7)
    tight_layout()

    # savefig("fig_selectivity_null_hist.png", dpi=600)
    savefig(dirfig * "fig_selectivity_null_hist.pdf", dpi=600)


    # # inhibitory neurons only
    # Ncells = 5000
    # Ne = 2500
    # ns_R = ns_R[Ne+1:Ncells]
    # ns_L = ns_L[Ne+1:Ncells]

    # # order inhibitory neurons by mean rate difference
    # cellsOrd = sortperm(rate_diff)

    # fig = figure(figsize=(3.0,1.5))
    # gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, wspace=0.2, hspace=0.5, left=0.15, bottom=0.3, right=0.9, width_ratios=[1,1,2])

    # ax = fig.add_subplot(gs[1])
    # rate_diff_1trial = (ns_R - ns_L) / 2
    # plot(rate_diff_1trial[cellsOrd], collect(1:length(cellsOrd)), color="gray", marker="o", ms=0.8, mec="None", linestyle="", alpha=1)
    # # plot(rate_diff[cellsOrd], collect(1:length(cellsOrd)), color="magenta", lw=0.8, linestyle="--", alpha=1)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)
    # xlim([-20,20])
    # # ylim([0,length(cellsOrd)])
    # xticks([-15,0,15], ["-15", "0","15"], fontsize=8)
    # yticks([])
    # # ax.set_xlabel("selectivity (Hz)", fontsize=8)
    # ax.set_ylabel("Null neuron", fontsize=8)
    # # title("single", fontsize=8)


    # ax = fig.add_subplot(gs[2])
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

    # savefig(dirfig * "null_RL_psth.png", dpi=600)
    # savefig(dirfig * "null_RL_psth.pdf", dpi=600)

end
