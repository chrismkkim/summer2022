function plt_pcor_pcomp(dirfig, pcor_R, pcor_L, pcomp_fs_R, pcomp_fs_L)


    pc1R = log10.(abs.(pcomp_fs_R[1,:] .- 15))
    pc1L = log10.(abs.(pcomp_fs_L[1,:] .- 5))

    # remove NaN
    idxR = .!isnan.(pcor_R)
    idxL = .!isnan.(pcor_L)
    idx = idxR .* idxL
    pcor = [pcor_R[idx]; pcor_L[idx]]
    pc1 = [pc1R[idx]; pc1L[idx]]

    # linear regression: pc1R, pc1L ~ pcor_R, pcor_L
    data = DataFrame(X=pcor, Y=pc1)
    ols = lm(@formula(Y ~ X), data)
    rsquared = r2(ols)
    intercept, slope = coef(ols)
    xrange = collect(-0.2:0.01:1)
    yrange = intercept .+ slope * xrange



    # plot figure
    figure(figsize=(1.5, 1.4))
    ax = subplot(111)
    plot(pcor_R, pc1R, marker="o", ms=1, mfc="None", mew=0.2, c="blue", alpha=1, linestyle="")
    plot(pcor_L, pc1L, marker="o", ms=1, mfc="None", mew=0.2, c="red", alpha=1, linestyle="")
    plot(xrange, yrange, c="black", lw=0.8, linestyle="--", label=L"$r^2=$" * "$(round(rsquared,digits=2))")
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)

    xlim([-1.15,1.15])
    # ylim([-3,6])
    xticks(fontsize=7)
    yticks(fontsize=7)
    xlabel("correlation", fontsize=7)
    ylabel("log10 |PC1|", fontsize=7)
    legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$r^2=$" * "$(round(rsquared,digits=2))")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=7, handlelength=0, bbox_to_anchor=(-0.1, 1.1))     
    tight_layout()
    
    # savefig("fs_pcor_pcomp.png", dpi=600)

    savefig(dirfig * "fs_pcor_pcomp.pdf", dpi=600)
    


    # figure(figsize=(3,3))
    # for i = 1:4
    #     subplot(2,2,i)
    #     plot(pcor_R, pcomp_fs_R[i,:],marker="o", ms=1, mec="None", c="blue", linestyle="")
    #     title("PC$(i)", fontsize=8)
    #     xlim([-1,1])
    #     ylim([-60,60])
    #     xticks(fontsize=8)
    #     yticks(fontsize=8)
    #     xlabel("corr", fontsize=8)
    #     ylabel("pcomp", fontsize=8)
    # end
    
    # tight_layout()
    
    # savefig(dirfig * "pca_pcor_fs_R.png", dpi=600)
    # savefig(dirfig * "pca_pcor_fs_R.pdf", dpi=600)
    
    
    
    # figure(figsize=(3,3))
    # for i = 1:4
    #     subplot(2,2,i)
    #     plot(pcor_L, pcomp_fs_L[i,:],marker="o", ms=1, mec="None", c="red", linestyle="")
    #     title("PC$(i)")
    #     xlim([-1,1])
    #     ylim([-60,60])
    #     xticks(fontsize=8)
    #     yticks(fontsize=8)
    #     xlabel("corr", fontsize=8)
    #     ylabel("pcomp", fontsize=8)
    # end
    
    # tight_layout()
    
    # savefig(dirfig * "pca_pcor_fs_L.png", dpi=600)
    # savefig(dirfig * "pca_pcor_fs_L.pdf", dpi=600)
    
    


end

