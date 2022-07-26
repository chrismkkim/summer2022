function plt_rate_dist_data(dirfig, pyr_R, pyr_L, fs_R, fs_L)

    pyr_R = mean(pyr_R, dims=1)[:]
    pyr_L = mean(pyr_L, dims=1)[:]
    fs_R = mean(fs_R, dims=1)[:]
    fs_L = mean(fs_L, dims=1)[:]

    pyrR_log = log10.(pyr_R)
    pyrL_log = log10.(pyr_L)
    fsR_log = log10.(fs_R)
    fsL_log = log10.(fs_L)
    # pyrR_log = pyrR_log[.!isinf.(pyrR_log)]
    # pyrL_log = pyrL_log[.!isinf.(pyrL_log)]
    # fsR_log = fsR_log[.!isinf.(fsR_log)]
    # fsL_log = fsL_log[.!isinf.(fsL_log)]


    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    cnt_E, bins_E = hist(pyrR_log, bins=50, range=(-2,2), lw=0.8, color="limegreen", histtype="step", density=true)
    cnt_I, bins_I = hist(fsR_log, bins=50, range=(-2,2), lw=0.8, color="darkorange", histtype="step", density=true)
    PYRmean = log10(mean(pyr_R))
    FSmean = log10(mean(fs_L))
    # axvline(PYRmean, color="limegreen", lw=0.8, linestyle="--")
    # axvline(FSmean, color="darkorange", lw=0.8, linestyle="--")

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="limegreen", lw=1.5, label="Pyr"),
                        Line2D([0], [0], color="darkorange", lw=1.5, label="FS")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     
    # xticks([0.01, 0.1, 1, 10, 100],fontsize=8)
    xticks([-2, -1, 0, 1, 2], ["-2","-1","0","1","2"], fontsize=8)
    yticks(fontsize=8)
    
    xlabel("log10 of spk/s", fontsize=8)
    ylabel("density", fontsize=8)
    tight_layout()

    # savefig(dirfig * "rate.png", dpi=600)
    savefig(dirfig * "rate_data.pdf", dpi=600)

end