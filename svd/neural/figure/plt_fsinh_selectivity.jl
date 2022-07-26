function plt_fsinh_selectivity(dirfig, fs_diff, fs_R, fs_L, inh_diff, inh_R, inh_L)


    fs_mean = mean((fs_R + fs_L)/2, dims=1)[:]
    fs_diffnorm = fs_diff ./ fs_mean

    inh_mean = mean((inh_R + inh_L)/2, dims=1)[:]
    inh_diffnorm = inh_diff ./ inh_mean


    figure(figsize=(1.5,1.4))
    ax = subplot(111)
    hist(fs_diffnorm, bins=30, range=(-2,2), lw=0.8, color="black", histtype="step", density=true)
    hist(inh_diffnorm, bins=100, range=(-2,2), lw=0.8, color="darkorange", histtype="step", density=true)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="black", lw=1.5, label="FS"),
                        Line2D([0], [0], color="darkorange", lw=1.5, label="Inh")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=7, handlelength=0.8)     
    xticks([-2,-1,0,1,2], ["-2","-1","0","1","2"],fontsize=7)
    yticks(fontsize=7)
    
    xlabel("selectivity", fontsize=7)
    ylabel("density", fontsize=7)
    tight_layout()

    # savefig("fig_selectivity_fsinh.png", dpi=600)
    savefig(dirfig * "fig_selectivity_fsinh.pdf", dpi=600)


end


