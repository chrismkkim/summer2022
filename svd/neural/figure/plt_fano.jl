function plt_fano(dirfig, ff_R, ff_L)

    ff_R = ff_R[ff_R .> 0]
    ff_L = ff_L[ff_L .> 0]

    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    cnt_R, bins_R = hist(ff_R, bins=60, range=(0,3), histtype="step", color="blue", linewidth=0.8)
    cnt_L, bins_L = hist(ff_L, bins=60, range=(0,3), histtype="step", color="red", linewidth=0.8)
    maxcnt_R = maximum(cnt_R); ffmean_R = mean(ff_R)
    maxcnt_L = maximum(cnt_L); ffmean_L = mean(ff_L)
    plot(ffmean_R*ones(11), collect(0:50:500), color="blue", linestyle="--", linewidth=0.8)
    plot(ffmean_L*ones(11), collect(0:50:500), color="red", linestyle="--", linewidth=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)

    legend_elements = [ Line2D([0], [0], color="blue", lw=1.5, label="right"),
                        Line2D([0], [0], color="red", lw=1.5, label="left")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8, loc="upper right", bbox_to_anchor=(1.0, 1.0))     
    xticks([0,1,2,3],fontsize=8)
    yticks([0,250, 500],fontsize=8)
    ylim([0,600])

    # ticklabel_format(axis="x", useMathText=true, useOffset=true)
    # ticklabel_format(axis="y", style="scientific", scilimits=(1,2), useMathText=true, useOffset=true)
    # ax.yaxis.offsetText.set_fontsize(8)
    xlabel("Fano factor", fontsize=8)
    ylabel("count", fontsize=8)
    tight_layout()


    savefig(dirfig * "fano.png", dpi=600)
    savefig(dirfig * "fano.pdf", dpi=600)

end

