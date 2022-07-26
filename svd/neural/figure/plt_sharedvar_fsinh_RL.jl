matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

function plt_sharedvar_fsinh_RL(dirfig, fname, 
    svar_cor_R, svarX_R, svarY_R, svarX_ksum_R, svarX_tsum_R, svarY_ksum_R, svarY_tsum_R,
    svar_cor_L, svarX_L, svarY_L, svarX_ksum_L, svarX_tsum_L, svarY_ksum_L, svarY_tsum_L,
    bal_cor_R, balX_R, balY_R, balX_ksum_R, balX_tsum_R, balY_ksum_R, balY_tsum_R,
    bal_cor_L, balX_L, balY_L, balX_ksum_L, balX_tsum_L, balY_ksum_L, balY_tsum_L)

    startind = 7
    timev = 20*collect(startind:100)
    
    ########################################
    #----- shared variance components -----#
    #----- trained balanced network   -----#
    ########################################
    # lick right
    fig=figure(figsize=(0.6,1.3))
    gs = gridspec.GridSpec(4, 1, figure=fig, wspace=0, hspace=0.01, left=0.3, bottom=0.0)
    ncomp = 4
    sgn = [1,1,1,-1]
    for ii = 1:ncomp
        ax = fig.add_subplot(gs[ii,1])
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)        
        xlim([0, 2000])
        ylim([-100, 70])

        plot(timev, sgn[ii]*svarX_R[:,ii], c="lightskyblue", lw=1.2, linestyle="-")
        plot(timev, sgn[ii]*svarY_R[:,ii], c="blue", lw=0.8, linestyle="--")
        if ii == 1
            text(0.2, 1.0, L"$\rho =$" * "$(round(svar_cor_R[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        else
            text(0.5, 0.8, "$(round(svar_cor_R[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        end
        ylabel("SC$(ii)", fontsize=6)
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -90*ones(length(tbar)), c="black", lw=1.5)
        end        
    end
    
    savefig(dirfig * "SVC_" * fname * "_R.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_R.pdf", dpi=600)
    
    
    # lick left
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

        plot(timev, sgn[ii]*svarX_L[:,ii], c="pink", lw=1.2, linestyle="-")
        plot(timev, sgn[ii]svarY_L[:,ii], c="red", lw=0.8, linestyle="--")
        if ii == 1
            text(0.2, 1.0, L"$\rho =$" * "$(round(svar_cor_L[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        else
            text(0.5, 0.8, "$(round(svar_cor_L[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        end
        ylabel("SC$(ii)", fontsize=6)
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -90*ones(length(tbar)), c="black", lw=1.5)
        end        
    end
    
    savefig(dirfig * "SVC_" * fname * "_L.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_L.pdf", dpi=600)
        
    

    ########################################
    #----- shared variance components -----#
    #----- initial balanced network   -----#
    ########################################
    # lick right
    fig=figure(figsize=(0.6,1.3))
    gs = gridspec.GridSpec(4, 1, figure=fig, wspace=0, hspace=0.01, left=0.3, bottom=0.0)
    ncomp = 4
    sgn = [1,1,1,-1]
    for ii = 1:ncomp
        ax = fig.add_subplot(gs[ii,1])
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)        
        xlim([0, 2000])
        ylim([-100, 70])

        plot(timev, sgn[ii]*balX_R[:,ii], c="lightskyblue", lw=1.2, linestyle="-")
        plot(timev, sgn[ii]*balY_R[:,ii], c="blue", lw=0.8, linestyle="--")
        if ii == 1
            text(0.2, 1.0, L"$\rho =$" * "$(round(bal_cor_R[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        else
            text(0.5, 0.8, "$(round(bal_cor_R[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        end
        ylabel("SC$(ii)", fontsize=6)
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -90*ones(length(tbar)), c="black", lw=1.5)
        end        
    end
    
    savefig(dirfig * "SVC_" * fname * "_R_bal.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_R_bal.pdf", dpi=600)
    
    
    # lick left
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

        plot(timev, sgn[ii]*balX_L[:,ii], c="pink", lw=1.2, linestyle="-")
        plot(timev, sgn[ii]balY_L[:,ii], c="red", lw=0.8, linestyle="--")
        if ii == 1
            text(0.2, 1.0, L"$\rho =$" * "$(round(bal_cor_L[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        else
            text(0.5, 0.8, "$(round(bal_cor_L[ii], digits=2))", fontsize=6, color="black", transform=ax.transAxes, clip_on=false)
        end
        ylabel("SC$(ii)", fontsize=6)
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -90*ones(length(tbar)), c="black", lw=1.5)
        end        
    end
    
    savefig(dirfig * "SVC_" * fname * "_L_bal.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_L_bal.pdf", dpi=600)



    ################################
    #----- correlation of SCs -----#
    ################################
    fig = figure(figsize=(2.1,1.4))
    ax = fig.add_subplot(111)
    ncomp = 40
    plot(collect(1:ncomp), svar_cor_R[1:ncomp], c="blue", marker="o", ms=2.0, lw=0.8, mec="None")
    plot(collect(1:ncomp), svar_cor_L[1:ncomp], c="red", marker="o", ms=2.0, lw=0.8, mec="None")
    ylim([-0.1,1.1])
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xlabel("SC", fontsize=7)
    ylabel("correlation " * L"($\rho$)", fontsize=7)
    xticks([0,20,40],fontsize=7)
    yticks(fontsize=7)    
    legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="Right"),
                       Line2D([0], [0], color="red", lw=1.5, label="Left")]
    legend(handles=legend_elements, frameon=false, fontsize=7, handlelength=1, bbox_to_anchor=(1.0, 1.0), loc="upper left")     
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_corall.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_corall.pdf", dpi=600)
    

    fig = figure(figsize=(2.1,1.4))
    ax = fig.add_subplot(111)
    ncomp = 40
    plot(collect(1:ncomp), bal_cor_R[1:ncomp], c="blue", marker="o", ms=2.0, lw=0.8, mec="None")
    plot(collect(1:ncomp), bal_cor_L[1:ncomp], c="red", marker="o", ms=2.0, lw=0.8, mec="None")
    ylim([-0.1,1.1])
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xlabel("SC", fontsize=7)
    ylabel("correlation " * L"($\rho$)", fontsize=7)
    xticks([0,20,40],fontsize=7)
    yticks(fontsize=7)    
    legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="Right"),
                       Line2D([0], [0], color="red", lw=1.5, label="Left")]
    legend(handles=legend_elements, frameon=false, fontsize=7, handlelength=1, bbox_to_anchor=(1.0, 1.0), loc="upper left")     
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_corall_bal.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_corall_bal.pdf", dpi=600)
    
    

    ##################################
    #----- explained var by SCs -----#
    ##################################    
    # lick right
    fig = figure(figsize=(1.5,1.4))
    ax = fig.add_subplot(111)
    ncomp = 10
    plot(collect(1:ncomp), svarX_ksum_R[1:ncomp]/svarX_tsum_R, c="lightskyblue", marker="o", ms=2, lw=0.8)
    plot(collect(1:ncomp), svarY_ksum_R[1:ncomp]/svarY_tsum_R, c="blue", marker="o", ms=2, lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ylim([5e-4,1])
    yscale("log")
    xlabel("SC", fontsize=7)
    ylabel("explained var", fontsize=7)
    legend_elements = [Line2D([0], [0], color="lightskyblue", lw=1.0, label="FS"),
                       Line2D([0], [0], color="blue", lw=1.0, label="Inh")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
    xticks(fontsize=7)
    yticks(fontsize=7)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_expvar_R.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_expvar_R.pdf", dpi=600)
    
    # lick left
    fig = figure(figsize=(1.5,1.4))
    ax = fig.add_subplot(111)
    ncomp = 10
    plot(collect(1:ncomp), svarX_ksum_L[1:ncomp]/svarX_tsum_L, c="pink", marker="o", ms=2, lw=0.8)
    plot(collect(1:ncomp), svarY_ksum_L[1:ncomp]/svarY_tsum_L, c="red", marker="o", ms=2, lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ylim([5e-4,1])
    yscale("log")
    xlabel("SC", fontsize=7)
    ylabel("explained var", fontsize=7)
    legend_elements = [Line2D([0], [0], color="pink", lw=1.5, label="FS"),
                       Line2D([0], [0], color="red", lw=1.5, label="Inh")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
    xticks(fontsize=7)
    yticks(fontsize=7)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_expvar_L.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_expvar_L.pdf", dpi=600)
        
    
    ###################################
    #----- cumulative var by SCs -----#
    ###################################    
    # lick right
    fig = figure(figsize=(1.5,1.4))
    ax = fig.add_subplot(111)
    ncomp = 10
    plot(collect(1:ncomp), cumsum(balY_ksum_R[1:ncomp]/balY_tsum_R), c="gray", marker="o", ms=2, mec="None", lw=0.8)
    plot(collect(1:ncomp), cumsum(svarX_ksum_R[1:ncomp]/svarX_tsum_R), c="lightskyblue", marker="o", ms=2, mec="None", lw=0.8)
    plot(collect(1:ncomp), cumsum(svarY_ksum_R[1:ncomp]/svarY_tsum_R), c="blue", marker="o", ms=2, mec="None", lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="lightskyblue", lw=1.5, label="FS"),
                       Line2D([0], [0], color="blue", lw=1.5, label="Inh"),
                       Line2D([0], [0], color="gray", lw=1.5, label="Null")]
    legend(loc="upper left", bbox_to_anchor=[0,1], columnspacing=0.4, borderaxespad=0.0, handles=legend_elements, frameon=false, fontsize=7, handlelength=1, ncol=2)     
    
    ylim([0,1])
    # yscale("log")
    xlabel("SC", fontsize=7)
    ylabel("cumulative var", fontsize=7)
    xticks([1, 3, 5, 7, 9], fontsize=7)
    yticks(fontsize=7)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_cumvar_R.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_cumvar_R.pdf", dpi=600)
    

    # lick left
    fig = figure(figsize=(1.5,1.4))
    ax = fig.add_subplot(111)
    ncomp = 10
    plot(collect(1:ncomp), cumsum(balY_ksum_L[1:ncomp]/balY_tsum_L), c="gray", marker="o", ms=2, mec="None", lw=0.8)
    plot(collect(1:ncomp), cumsum(svarX_ksum_L[1:ncomp]/svarX_tsum_L), c="pink", marker="o", ms=2, mec="None", lw=0.8)
    plot(collect(1:ncomp), cumsum(svarY_ksum_L[1:ncomp]/svarY_tsum_L), c="red", marker="o", ms=2, mec="None", lw=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    legend_elements = [Line2D([0], [0], color="pink", lw=1.5, label="FS"),
                       Line2D([0], [0], color="red", lw=1.5, label="Inh"),
                       Line2D([0], [0], color="gray", lw=1.5, label="Null")]
    legend(loc="upper left", bbox_to_anchor=[0,1], columnspacing=0.4, borderaxespad=0.1, handles=legend_elements, frameon=false, fontsize=7, handlelength=1, ncol=2)     
    
    ylim([0,1])
    # yscale("log")
    xlabel("SC", fontsize=7)
    ylabel("cumulative var", fontsize=7)
    xticks([1, 3, 5, 7, 9], fontsize=7)
    yticks(fontsize=7)
    tight_layout()
    
    savefig(dirfig * "SVC_" * fname * "_cumvar_L.png", dpi=600)
    savefig(dirfig * "SVC_" * fname * "_cumvar_L.pdf", dpi=600)
    

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
#         text(1500, -80, "500ms", fontsize=7)
#     end        
# end

# tight_layout()

# savefig(dirfig * "SVC" * fname * ".png", dpi=600)





    
# fig=figure(figsize=(0.4,1.3))
# ncomp = 4
# for ii = 1:ncomp
#     ax = fig.add_subplot(4,1,ii)
#     plot(svarX[:,ii], svarY[:,ii], c="magenta", marker="o", ms=0.5, mec="None", linestyle="", label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")
#     # title("SVC $(ii)", fontsize=7)
#     xticks(fontsize=7)
#     yticks(fontsize=7)
#     ax.spines["top"].set_visible(false)
#     ax.spines["right"].set_visible(false)

#     # if ii > 1
#         xticks(color="None")
#         yticks(color="None")
#     # end
#     # xlim([-2, 1.5])
#     # ylim([-4, 2.5])
#     legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")]
#     legend(loc=4, handles=legend_elements, frameon=false, fontsize=7, handlelength=0)     
# end
# # annotate("FS cells", xy=(0.45,0.03), xycoords="figure fraction", fontsize=12)
# # annotate("Inh model", xy=(0.01,0.5), xycoords="figure fraction", fontsize=7, rotation="vertical", annotation_clip=true)
# # ax = fig.add_subplot(111)
# # xlabel("FS cells", fontsize=12, labelpad=10)
# # ylabel("Inh model", fontsize=12, labelpad=10)
# # ax.set_frame_on(false)
# # ax.tick_params(left=false, bottom=false)
# # ax.tick_params(labelleft=false, labelbottom=false)


# tight_layout()