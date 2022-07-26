function plt_pyrexc_selectivity(dirfig, pyr_diff, pyr_R, pyr_L, exc_diff, exc_R, exc_L, matchedCellsOrd)

    pyr_mean = mean((pyr_R + pyr_L)/2, dims=1)[:]
    pyr_diffnorm = pyr_diff ./ pyr_mean

    exc_mean = mean((exc_R + exc_L)[:,matchedCellsOrd]/2, dims=1)[:]
    exc_diffnorm = exc_diff[matchedCellsOrd] ./ exc_mean


    figure(figsize=(1.5,1.4))
    ax = subplot(111)
    hist(pyr_diffnorm, bins=100, range=(-2,2), lw=0.8, color="black", histtype="step", density=true, clip_on=false)
    hist(exc_diffnorm, bins=100, range=(-2,2), lw=0.8, color="limegreen", histtype="step", density=true, clip_on=false)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="black", lw=1.5, label="Pyr"),
                        Line2D([0], [0], color="limegreen", lw=1.5, label="Exc")]
    legend(handles=legend_elements, frameon=false, fontsize=7, handlelength=0.8, bbox_to_anchor=(0.00, 0.6),loc="lower left")     
    xticks([-2,-1,0,1,2], ["-2","-1","0","1","2"],fontsize=7)
    yticks(fontsize=7)
    
    xlabel("selectivity", fontsize=7)
    ylabel("density", fontsize=7)
    tight_layout()

    # savefig("fig_selectivity_pyrexc.png", dpi=600)
    savefig(dirfig * "fig_selectivity_pyrexc.pdf", dpi=600)


end


