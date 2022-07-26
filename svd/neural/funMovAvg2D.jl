function funMovAvg2D(x,wid)

    if wid > 0
        Nsteps = size(x)[1]
        movavg = zeros(size(x))
        for i = 1:Nsteps
            Lind = maximum([i-wid, 1])
            Rind = minimum([i+wid, Nsteps])
            xslice = @view x[Lind:Rind,:]
            movavg[i,:] = mean(xslice, dims=1)
        end
            
        return movavg

    elseif wid == 0

        return x
        
    end

end