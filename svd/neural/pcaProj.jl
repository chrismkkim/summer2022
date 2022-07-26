function pcaProj(usum, utargFS, explainedVar)

    # project to subspace
    M = fit(PCA, utargFS; pratio=explainedVar)
    pcomp = transform(M, usum)
    usum_proj = reconstruct(M,pcomp)

    return usum_proj

end