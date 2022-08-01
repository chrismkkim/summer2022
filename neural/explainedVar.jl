function explainedVar(f, y)
    # mean centered
    f = f .- mean(f)
    y = y .- mean(y)

    # explained variance
    res = y - f
    stotal = var(y)
    sres = var(res)
    expvar = 1 - sres/stotal

    return expvar, stotal, sres

end