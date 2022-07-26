function plt_rate_dist_model(dirfig, mrate)

    Ne = 2500

    Erate = mrate[1:Ne]
    Irate = mrate[Ne+1:end]

    Erate_log = log10.(Erate)
    Irate_log = log10.(Irate)
    Erate_log = Erate_log[.!isinf.(Erate_log)]
    Irate_log = Irate_log[.!isinf.(Irate_log)]

    logbins = 10 .^ range(-2,2,length=101)
    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    cnt_E, bins_E = hist(Erate_log, bins=100, range=(-2,2), lw=0.8, color="limegreen", histtype="step", density=true)
    cnt_I, bins_I = hist(Irate_log, bins=100, range=(-2,2), lw=0.8, color="darkorange", histtype="step", density=true)
    # cnt_E, bins_E = hist(Erate, bins=logbins, lw=0.8, color="limegreen", histtype="step")
    # cnt_I, bins_I = hist(Irate, bins=logbins, lw=0.8, color="darkorange", histtype="step")
    # xscale("log")
    # ax.tick_params(axis="y", which="minor")
    pmax_E = maximum(cnt_E)
    pmax_I = maximum(cnt_I)
    Emean = mean(Erate)
    Imean = mean(Irate)
    # plot(Emean*ones(11), collect(0:pmax_E/10:pmax_E), color="limegreen", lw=0.8, linestyle="--")
    # plot(Imean*ones(11), collect(0:pmax_I/10:pmax_I), color="darkorange", lw=0.8, linestyle="--")

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="limegreen", lw=1.5, label="Exc"),
                        Line2D([0], [0], color="darkorange", lw=1.5, label="Inh")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     
    # xticks([0.01, 0.1, 1, 10, 100],fontsize=8)
    xticks([-2, -1, 0, 1, 2], ["-2","-1","0","1","2"], fontsize=8)
    yticks(fontsize=8)
    
    xlabel("log10 of spk/s", fontsize=8)
    ylabel("density", fontsize=8)
    tight_layout()

    # savefig(dirfig * "rate.png", dpi=600)
    savefig(dirfig * "rate_model.pdf", dpi=600)

end