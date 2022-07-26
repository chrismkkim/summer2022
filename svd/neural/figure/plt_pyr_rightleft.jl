matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

clrGray = [cmap.Greys(x) for x in [0.4, 0.8]]
clrGreen = [cmap.Greens(x) for x in [0.4, 0.8]]

function plt_pyr_rightleft(dirfig, pyr_R, pyr_L)
    
    pyrDiff = mean(pyr_L - pyr_R, dims=1)[:]
    pyrsorted = sortperm(pyrDiff)

    timev = 20*collect(1:100)
    
    for ngrp = 1:18
        println(ngrp)
        figure(figsize=(16,16))
        for ii = 1:100
            ax = subplot(10,10,ii)
            nid = pyrsorted[(ngrp-1)*100 + ii]
            plot(timev, pyr_R[:,nid], color="blue")
            plot(timev, pyr_L[:,nid], color="red")
            ax.spines["top"].set_visible(false)
            ax.spines["right"].set_visible(false)            
            title("$(round(pyrDiff[nid], digits=1))")
        end
        tight_layout()
        savefig(dirfig * "pyr_RL$(ngrp).png", dpi=300)
        close()
    end

    

    # figure(figsize=(3.5,3))
    # # # ax = subplot(211)
    # plot(excDiff_1trial, color="blue", marker="x", ms=1, linestyle="", alpha=0.5)
    # plot(excDiff[sortperm(excDiff)], color="blue", marker="o", ms=1, linestyle="")
    # plot(collect(1:length(excDiff)), zeros(length(excDiff)), color="gray", linestyle="--")
    # # plot(fsDiff, color="darkorange", marker="o", ms=2, linestyle="")
    # xlabel("Exc neuron", fontsize=12)
    # ylabel(L"$\Delta rate (Hz)$", fontsize=12)
    # xticks(fontsize=12)
    # yticks(fontsize=12)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)

    # # ax = subplot(212)
    # # plot(pyrDiff, color="limegreen", marker="o", ms=2, linestyle="")
    # plot(inhDiff_1trial, color="red", marker="x", ms=1, linestyle="", alpha=0.5)
    # plot(inhDiff[sortperm(inhDiff)], color="red", marker="o", ms=1, linestyle="")
    # plot(collect(1:length(inhDiff)), zeros(length(inhDiff)), color="gray", linestyle="--")
    # xlabel("Inh neuron", fontsize=12)
    # ylabel(L"$\Delta rate (Hz)$", fontsize=12)
    # xticks(fontsize=12)
    # yticks(fontsize=12)
    # ax.spines["top"].set_visible(false)
    # ax.spines["right"].set_visible(false)

    # tight_layout()

    # savefig(dirfig * "selectivity_ei.png", dpi=300)



    # for ii = 1:length(navg_ind)
    #     ax = subplot(2,1,ii)
    #     cnt_R, bins_R = hist(pcor_R[:,ii], bins=100, range=(-1,1), color="blue", histtype="step")
    #     cnt_L, bins_L = hist(pcor_L[:,ii], bins=100, range=(-1,1), color="red", histtype="step")
    #     pmax_R = maximum(cnt_R)
    #     pmax_L = maximum(cnt_L)
    #     ptmpR = pcor_R[:,ii]
    #     ptmpL = pcor_L[:,ii]
    #     medR = median(ptmpR[.!isnan.(ptmpR)])
    #     medL = median(ptmpL[.!isnan.(ptmpL)])
    #     plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", linestyle="--")
    #     plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", linestyle="--")
    #     title(L"$N_{trial}=$" * "$(navg_list[navg_ind[ii]])", fontsize=12)

    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     if ii == 1
    #         legend_elements = [Line2D([0], [0], color="blue", lw=2, label="right"),
    #                            Line2D([0], [0], color="red", lw=2, label="left")]
    #         legend(loc=2, handles=legend_elements, frameon=false, fontsize=10, handlelength=1)     
    #         tick_params(labelbottom=false)
    #         yticks(fontsize=12)
    #     else
    #         xticks(fontsize=12)
    #         yticks(fontsize=12)
    #     end

    # end
    # xlabel("correlation", fontsize=12)
    # ylabel("count", fontsize=12)
    # tight_layout()

    # savefig(dirfig * "pyr_pcor.png", dpi=300)
     
end
    