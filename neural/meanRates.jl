function meanRates(p, times)

    exc = collect(1:p.Ne)
    inh = collect(p.Ne+1:p.Ncells)

    times_exc = times[exc,:]
    times_inh = times[inh,:]

    idx1 = (times_exc .> p.stim_off) .* (times_exc .< p.train_time)
    idx2 = (times_inh .> p.stim_off) .* (times_inh .< p.train_time)

    rate_exc = sum(idx1) / p.train_duration * 1000 / length(exc)
    rate_inh = sum(idx2) / p.train_duration * 1000 / length(inh)

    return rate_exc, rate_inh
end