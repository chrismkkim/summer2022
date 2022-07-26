function explainedVar(y, f)
    # y : data
    # f : target    

    # # mean centered
    # f = f .- mean(f)
    # y = y .- mean(y)

    # # explained variance
    # res = y - f
    # stotal = var(f)
    # sres = var(res)
    # expvar = 1 - sres/stotal

    # mse
    res = y - f
    mse = 1 - sum(res.^2) / sum(f.^2)
    
    return mse

end