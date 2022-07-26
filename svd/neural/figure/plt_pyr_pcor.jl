matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

clrGray = [cmap.Greys(x) for x in [0.4, 0.8]]
clrGreen = [cmap.Greens(x) for x in [0.4, 0.8]]

function plt_pyr_pcor(dirfig, dirsim, navg_list, pyr_R, pyr_L, almOrd, matchedCells)
    
    navg_ind = [1, 4]
    pcor_R = zeros(length(matchedCells), length(navg_ind))
    pcor_L = zeros(length(matchedCells), length(navg_ind))

    for (idx, ii) in enumerate(navg_ind)
        navg = navg_list[ii]

        # load usum_navg, rate_navg
        fname_rate_R = dirsim * "rate_R_navg$(navg).jld"
        fname_rate_L = dirsim * "rate_L_navg$(navg).jld"
        rate_R = load(fname_rate_R, "rate")
        rate_L = load(fname_rate_L, "rate")

        # performance: Pyr - Exc
        nrep = 20
        pcor_R[:,idx] = runperformance_nrep(pyr_R, rate_R, almOrd, matchedCells, nrep)
        pcor_L[:,idx] = runperformance_nrep(pyr_L, rate_L, almOrd, matchedCells, nrep)
    end

    figure(figsize=(1.7,1.5))
    ax = subplot(111)
    ii = 2 # navg = 400 (number of trial averages)
    cnt_R, bins_R = hist(pcor_R[:,ii], bins=50, range=(-1,1), lw=0.8, color="blue", histtype="step")
    cnt_L, bins_L = hist(pcor_L[:,ii], bins=50, range=(-1,1), lw=0.8, color="red", histtype="step")
    pmax_R = maximum(cnt_R)
    pmax_L = maximum(cnt_L)
    ptmpR = pcor_R[:,ii]
    ptmpL = pcor_L[:,ii]
    medR = median(ptmpR[.!isnan.(ptmpR)])
    medL = median(ptmpL[.!isnan.(ptmpL)])
    plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", lw=0.8, linestyle="--")
    plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", lw=0.8, linestyle="--")
    # title(L"$N_{trial}=$" * "$(navg_list[navg_ind[ii]])", fontsize=12)

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="right"),
                        Line2D([0], [0], color="red", lw=1.5, label="left")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=8, handlelength=0.8)     
    xticks([-1,-0.5,0,0.5,1], ["-1","-0.5","0","0.5","1"],fontsize=8)
    yticks(fontsize=8)
    
    xlabel("performance", fontsize=8)
    ylabel("count", fontsize=8)
    tight_layout()

    savefig(dirfig * "pyr_pcor.png", dpi=600)
    savefig(dirfig * "pyr_pcor.pdf", dpi=600)


    return pcor_R[:,ii], pcor_L[:,ii]


    # figure(figsize=(3.5,3))
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
    