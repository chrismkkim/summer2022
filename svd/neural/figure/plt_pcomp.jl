function plt_pcomp(dirfig, fname, startind, M_L, M_R)


timev = 20*collect(startind:100)

fig=figure(figsize=(3.5,3.0))
ncomp = 9
# sgn_pyr = [-1,-1,-1,-1]
if fname == "pyr"
    sgn_L = [-1, 1, -1, 1, 1, -1, -1, 1, 1]
    sgn_R = [-1, 1, -1, 1, -1, 1, -1, 1, -1]
elseif fname == "fs"
    sgn_L = [-1, 1, -1,-1, -1, 1, 1, 1, -1]
    sgn_R = [-1, 1, -1, 1, -1, -1, -1, 1, 1]
elseif fname == "exc"
    sgn_L = [-1, 1, 1, 1, 1, 1, -1, 1, -1]
    sgn_R = [-1, 1, -1, 1, -1, -1, -1, 1, 1]
elseif fname == "inh"
    sgn_L = [-1, 1, 1,-1, -1, 1, 1, 1, -1]
    sgn_R = [-1, 1, -1, 1, 1, -1, -1, 1, 1]
elseif fname == "bal"
    sgn_L = [-1, 1, 1,-1, -1, 1, 1, 1, -1]
    sgn_R = [-1, 1, -1, 1, 1, -1, -1, 1, 1]
end

for ii = 1:ncomp
    ax = fig.add_subplot(3,3,ii)
    ax.set_frame_on(false)
    ax.tick_params(left=false, bottom=false)
    ax.tick_params(labelleft=false, labelbottom=false)
    title("PC $(ii)", fontsize=8)
    xticks(color="None", fontsize=0)
    yticks(color="None", fontsize=0)
    xlim([0, 2000])
    ylim([-0.25, 0.25])

    plot(timev, sgn_L[ii]*M_L.proj[:,ii], c="red")
    plot(timev, sgn_R[ii]*M_R.proj[:,ii], c="blue")
    
end

tight_layout()

# savefig(dirfig * "pca_pc_" * fname * "_smooth.png", dpi=600)
savefig(dirfig * "pca_pc_" * fname * "_smooth.pdf", dpi=600)





ncomp = 9
expvar_L = cumsum(principalvars(M_L))[1:ncomp] / tvar(M_L)
expvar_R = cumsum(principalvars(M_R))[1:ncomp] / tvar(M_R)

figure(figsize=(1.7,1.5))
plot(collect(1:ncomp), expvar_L, color="red", marker="o", lw=0.8, ms=2)
plot(collect(1:ncomp), expvar_R, color="blue", marker="o", lw=0.8, ms=2)
legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="Right"),
                   Line2D([0], [0], color="red", lw=1.5, label="Left")]
legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, loc=4)     

xlabel("PC", fontsize=8)
ylabel("cumulative var", fontsize=8)
xticks([1, 3, 5, 7, 9], fontsize=8)
yticks(fontsize=8)
ylim([0,1])
tight_layout()
# savefig(dirfig * fname * "_cumvar.png", dpi=600)
# savefig(dirfig * "pca_cumvar_" * fname * ".png", dpi=600)
savefig(dirfig * "pca_cumvar_" * fname * "_smooth.pdf", dpi=600)


end
