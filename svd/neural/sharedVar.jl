function sharedVar(svdCor, X, Y)

    nsvd = size(svdCor.U)[2]
    nstep = size(X)[1]
    svarX = zeros(nstep,nsvd)
    svarY = zeros(nstep,nsvd)
    for ii = 1:nsvd
        svarX[:,ii] = X * svdCor.U[:,ii]
        svarY[:,ii] = Y * svdCor.V[:,ii]
    end
    svarX_ksum = sum(svarX.^2, dims=1)[:]
    svarY_ksum = sum(svarY.^2, dims=1)[:]
    svarX_tsum = sum(X.^2)
    svarY_tsum = sum(Y.^2)
    
    svar_cor = zeros(nsvd)
    for ii = 1:nsvd
        svar_cor[ii] = cor(svarX[:,ii], svarY[:,ii])
    end
    
    return svar_cor, svarX, svarY, svarX_ksum, svarY_ksum, svarX_tsum, svarY_tsum

end

