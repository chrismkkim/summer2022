function funSample(p,X)

    stim_off = 1000.0
    Nsteps = Int(2000.0/p.dt)
    timev = collect(1:Nsteps)*p.dt
    idx = timev .>= stim_off + p.learn_every

    # timev = collect(1:p.Nsteps)*p.dt
    # idx = timev .>= p.stim_off + p.learn_every

    if ndims(X) == 1
        XpostStim = X[idx]
        Xsampled = XpostStim[1:p.learn_step:end]
    elseif ndims(X) == 2
        XpostStim = X[idx,:]
        Xsampled = XpostStim[1:p.learn_step:end,:]
    end

    return Xsampled

end