function plt_cormat(dirfig, fname, cor_fsinh, svd_fsinh)

svnorm = svd_fsinh.S / sum(svd_fsinh.S)

fig = figure(figsize=(1.7,1.5))
ax = fig.add_subplot(111)
lensvd = length(svnorm)
nsing = 20
plot(collect(1:nsing), svnorm[1:nsing], color="black", marker="o", ms=2, lw=0.8)
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
xlabel("rank", fontsize=8)
ylabel("singular value", fontsize=8)
xticks(fontsize=8)
yticks(fontsize=8)
# legend_elements = [Line2D([0], [0], color="red", lw=2, label="trained"),
#                    Line2D([0], [0], color="black", lw=2, label="not trained")]
# legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
tight_layout()

savefig(dirfig * "cormat_sv_" * fname * ".png", dpi=600)
savefig(dirfig * "cormat_sv_" * fname * ".pdf", dpi=600)



sortrow = sortperm(svd_fsinh.U[:,1])
fig = figure(figsize=(1.8,1.5))
ax = fig.add_subplot(111)
imshow(cor_fsinh[sortrow,:], cmap="bwr", vmin=-1,vmax=1, interpolation="None", aspect="auto")
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
xticks(fontsize=8)
yticks(fontsize=8)
# if fname == "fsinh"
xlabel("Inh neuron", fontsize=8)
ylabel("FS cell", fontsize=8)
# elseif fname == "pyrexc"
#     xlabel("Exc", fontsize=8)
#     ylabel("Pyr", fontsize=8)
# end
cbar = colorbar(ticks=[-1,0,1])
cbar.ax.tick_params(labelsize=8)
# cbar.ax.set_yticklabels(['< -1', '0', '> 1'])
tight_layout()

savefig(dirfig * "cormat_" * fname * ".png", dpi=600)
savefig(dirfig * "cormat_" * fname * ".pdf", dpi=600)



end

