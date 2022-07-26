function plt_syn_dist(dirfig, uplas_avg_R, ubal_avg_R, uebal_avg_R, uibal_avg_R,
    uplas_avg_L, ubal_avg_L, uebal_avg_L, uibal_avg_L)

    figure(figsize=(2.3,1.5))
    ax = subplot(111)
    uplas_avg_R_exc = uplas_avg_R[1:p.Ne]
    uplas_avg_R_exc = uplas_avg_R_exc[abs.(uplas_avg_R_exc) .> 0]
    hist(ubal_avg_R[1:p.Ne], bins=1000, range=(-1,5), color="black", histtype="step", log=true, lw=0.5, label=L"$u^{bal}$")
    hist(-uibal_avg_R[1:p.Ne], bins=1000, range=(-1,5), color="darkorange", histtype="step", log=true, lw=0.5, label=L"$u^{bal}_I$")
    hist(uebal_avg_R[1:p.Ne], bins=1000, range=(-1,5), color="limegreen", histtype="step", log=true, lw=0.8, label=L"$u^{bal}_E$")
    hist(uplas_avg_R_exc, bins=1000, range=(-1,5), color="cyan", histtype="step", log=true, lw=0.8, label=L"$u^{plas}$")
    # hist(uext_avg_R[1:p.Ne], bins=1000, range=(-1,5), color="gray", histtype="step", log=true, lw=0.8, label=L"$u^{ext}$")
    plot(ones(11), collect(0:30:300), color="magenta", linestyle="--", lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [ Line2D([0], [0], color="limegreen", lw=1.5, label=L"$u_{bal}^E$"),
                        Line2D([0], [0], color="darkorange", lw=1.5, label=L"$u_{bal}^I$"),
                        Line2D([0], [0], color="black", lw=1.5, label=L"$u_{bal}$"),
                        Line2D([0], [0], color="cyan", lw=1.5, label=L"$u_{plas}$")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
    xticks(fontsize=8)
    yticks([10, 100], fontsize=8)
    xlabel("synaptic inputs", fontsize=8)
    ylabel("count", fontsize=8)
    tight_layout()

    savefig(dirfig * "udist_R.png", dpi=600)
    savefig(dirfig * "udist_R.pdf", dpi=600)



    figure(figsize=(2.3,1.5))
    ax = subplot(111)
    uplas_avg_L_exc = uplas_avg_L[1:p.Ne]
    uplas_avg_L_exc = uplas_avg_L_exc[abs.(uplas_avg_L_exc) .> 0]
    hist(ubal_avg_L[1:p.Ne], bins=1000, range=(-1,5), color="black", histtype="step", log=true, lw=0.5, label=L"$u^{bal}$")
    hist(-uibal_avg_L[1:p.Ne], bins=1000, range=(-1,5), color="darkorange", histtype="step", log=true, lw=0.5, label=L"$u^{bal}_I$")
    hist(uebal_avg_L[1:p.Ne], bins=1000, range=(-1,5), color="limegreen", histtype="step", log=true, lw=0.8, label=L"$u^{bal}_E$")
    hist(uplas_avg_L_exc, bins=1000, range=(-1,5), color="cyan", histtype="step", log=true, lw=0.8, label=L"$u^{plas}$")
    plot(ones(11), collect(0:30:300), color="magenta", linestyle="--", lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [ Line2D([0], [0], color="limegreen", lw=1.5, label=L"$u_{bal}^E$"),
                        Line2D([0], [0], color="darkorange", lw=1.5, label=L"$u_{bal}^I$"),
                        Line2D([0], [0], color="black", lw=1.5, label=L"$u_{bal}$"),
                        Line2D([0], [0], color="cyan", lw=1.5, label=L"$u_{plas}$")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
    xticks(fontsize=8)
    yticks([10, 100], fontsize=8)
    xlabel("synaptic inputs", fontsize=8)
    ylabel("count", fontsize=8)
    tight_layout()

    savefig(dirfig * "udist_L.png", dpi=600)
    savefig(dirfig * "udist_L.pdf", dpi=600)



end