function poprate(p, times, ns, cellsTrained)

    twin = 10
    nstep = Int(p.train_time / twin) + 1
    timev = collect(1:nstep)
    rate = zeros(length(timev), p.Ncells)

    for ci = 1:p.Ncells
        for spk = 1:ns[ci]
            idx = Int(div(times[ci,spk], twin) + 1)
            rate[idx, ci] += 1 / (twin/1000)
        end
    end

    rate_exc = mean(rate[:, 1:p.Ne], dims=2)
    rate_inh = mean(rate[:, p.Ne+1:end], dims=2)

    return rate_exc, rate_inh

end
