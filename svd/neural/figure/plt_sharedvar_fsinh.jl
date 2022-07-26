matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

function plt_sharedvar_fsinh(dirfig, fname, svar_cor, svarX, svarY, svarX_ksum, svarX_tsum, svarY_ksum, svarY_tsum)

    startind = 7
    timev = 20*collect(startind:100)
    
    fig=figure(figsize=(0.6,1.3))
    gs = gridspec.GridSpec(4, 1, figure=fig, wspace=0, hspace=0.01, left=0.3, bottom=0.0)
    ncomp = 4
    sgn = [-1,1,1,1]
    for ii = 1:ncomp
        ax = fig.add_subplot(gs[ii,1])
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)        
        xlim([0, 2000])
        ylim([-100, 70])

        plot(timev, svarX[:,ii], c="pink", lw=1.2, linestyle="-")
        plot(timev, svarY[:,ii], c="red", lw=0.8, linestyle="--")
        if ii == 1
            text(0.2, 1.0, L"$\rho =$" * "$(round(svar_cor[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        else
            text(0.5, 0.8, "$(round(svar_cor[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        end
        ylabel("SC$(ii)", fontsize=6)
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -90*ones(length(tbar)), c="black", lw=1.5)
        end        
    end
    
    # tight_layout()
    
    savefig(dirfig * "SVC_" * fname * ".png", dpi=600)
    savefig(dirfig * "SVC_" * fname * ".pdf", dpi=600)
    
    
    
    
    # fig=figure(figsize=(0.4,1.3))
    # gs = gridspec.GridSpec(4, 1, figure=fig, wspace=0.0, hspace=0.4, left=0.0, bottom=0.0)
    # ncomp = 4
    # for ii = 1:ncomp
    #     ax = fig.add_subplot(gs[ii,1])
    #     plot(svarX[:,ii], svarY[:,ii], c="magenta", marker="o", ms=1.0, mec="None", linestyle="", label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")
    #     # title("SVC $(ii)", fontsize=8)
    #     # xticks(fontsize=8)
    #     # yticks(fontsize=8)
    #     tick_params(labelleft=false, labelbottom=false)
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)

    #     legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")]
    #     legend(loc=2, handles=legend_elements, frameon=false, fontsize=6, handlelength=0)     
    # end
    # savefig(dirfig * "SVC_" * fname * "_cor.png", bbox_inches="tight", dpi=600)
    # savefig(dirfig * "SVC_" * fname * "_cor.pdf", bbox_inches="tight", dpi=600)
    
    
    
    fig = figure(figsize=(1.7,1.5))
    ax = fig.add_subplot(111)
    ncomp = 40
    plot(collect(1:ncomp), svar_cor[1:ncomp], c="red", marker="o", ms=2.0, lw=0.8, mec="None")
    ylim([-0.1,1.1])
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xlabel("SC", fontsize=8)
    ylabel("correlation " * L"($\rho$)", fontsize=8)
    xticks([0,20,40],fontsize=8)
    yticks(fontsize=8)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_corall.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_corall.pdf", dpi=600)
    
    
    
    fig = figure(figsize=(1.7,1.5))
    ax = fig.add_subplot(111)
    ncomp = 20
    plot(collect(1:ncomp), svarX_ksum[1:ncomp]/svarX_tsum, c="pink", marker="o", ms=2, lw=0.8)
    plot(collect(1:ncomp), svarY_ksum[1:ncomp]/svarY_tsum, c="red", marker="o", ms=2, lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ylim([5e-4,1])
    yscale("log")
    xlabel("SC", fontsize=8)
    ylabel("explained var", fontsize=8)
    legend_elements = [Line2D([0], [0], color="pink", lw=1.5, label="FS"),
                       Line2D([0], [0], color="red", lw=1.5, label="Inh")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=8, handlelength=1)     
    xticks(fontsize=8)
    yticks(fontsize=8)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_expvar.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_expvar.pdf", dpi=600)
    
    
    
    fig = figure(figsize=(1.7,1.5))
    ax = fig.add_subplot(111)
    ncomp = 20
    plot(collect(1:ncomp), cumsum(svarX_ksum[1:ncomp]/svarX_tsum), c="pink", marker="o", ms=2, mec="None", lw=0.8)
    plot(collect(1:ncomp), cumsum(svarY_ksum[1:ncomp]/svarY_tsum), c="red", marker="o", ms=2, mec="None", lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="pink", lw=2, label="FS"),
                       Line2D([0], [0], color="red", lw=2, label="Inh")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=8, handlelength=1)     
    
    ylim([0,1])
    # yscale("log")
    xlabel("SC", fontsize=8)
    ylabel("cumulative var", fontsize=8)
    xticks([0,5,10,15,20], fontsize=8)
    yticks(fontsize=8)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_cumvar.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_cumvar.pdf", dpi=600)
    
    

end






# startind = 7
# timev = 20*collect(startind:100)

# fig=figure(figsize=(2.5,4.5))
# ncomp = 4
# sgn = [-1,1,1,1]
# for ii = 1:ncomp
#     ax = fig.add_subplot(4,1,ii)
#     ax.set_frame_on(false)
#     ax.tick_params(left=false, bottom=false)
#     ax.tick_params(labelleft=false, labelbottom=false)
#     ylabel("SVC $(ii)", fontsize=12)
#     xlim([0, 2000])
#     ylim([-100, 70])

#     plot(timev, svarX[:,ii], c="darkorange", lw=2)
#     plot(timev, svarY[:,ii], c="blue", lw=2, linestyle="--")
#     if ii == 1
#         legend_elements = [Line2D([0], [0], color="darkorange", lw=2, label="FS"),
#                            Line2D([0], [0], color="blue", lw=2, label="Inh")]
#         legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
#     end
#     if ii == 4
#         tbar = 20*collect(75:100)
#         plot(tbar, -90*ones(length(tbar)), c="black", lw=2)
#         text(1500, -80, "500ms", fontsize=8)
#     end        
# end

# tight_layout()

# savefig(dirfig * "SVC" * fname * ".png", dpi=600)





    
# fig=figure(figsize=(0.4,1.3))
# ncomp = 4
# for ii = 1:ncomp
#     ax = fig.add_subplot(4,1,ii)
#     plot(svarX[:,ii], svarY[:,ii], c="magenta", marker="o", ms=0.5, mec="None", linestyle="", label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")
#     # title("SVC $(ii)", fontsize=8)
#     xticks(fontsize=8)
#     yticks(fontsize=8)
#     ax.spines["top"].set_visible(false)
#     ax.spines["right"].set_visible(false)

#     # if ii > 1
#         xticks(color="None")
#         yticks(color="None")
#     # end
#     # xlim([-2, 1.5])
#     # ylim([-4, 2.5])
#     legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")]
#     legend(loc=4, handles=legend_elements, frameon=false, fontsize=8, handlelength=0)     
# end
# # annotate("FS cells", xy=(0.45,0.03), xycoords="figure fraction", fontsize=12)
# # annotate("Inh model", xy=(0.01,0.5), xycoords="figure fraction", fontsize=8, rotation="vertical", annotation_clip=true)
# # ax = fig.add_subplot(111)
# # xlabel("FS cells", fontsize=12, labelpad=10)
# # ylabel("Inh model", fontsize=12, labelpad=10)
# # ax.set_frame_on(false)
# # ax.tick_params(left=false, bottom=false)
# # ax.tick_params(labelleft=false, labelbottom=false)


# tight_layout()