using DataFrames
using CSV
using PyPlot
using Statistics

dirpath = "/home/grabelmz/data/trainingData/"
fname1 = "target_fni16_exc_left.csv"
fname2 = "target_fni16_inh_left.csv"
fname3 = "target_fni16_exc_right.csv"
fname4 = "target_fni16_inh_right.csv"

target_16_exc_left = Matrix(DataFrame(CSV.File(dirpath * fname1)))
target_16_inh_left = Matrix(DataFrame(CSV.File(dirpath * fname2)))
target_16_exc_right = Matrix(DataFrame(CSV.File(dirpath * fname3)))
target_16_inh_right = Matrix(DataFrame(CSV.File(dirpath * fname4)))

target_16_exc_left = target_16_exc_left[:, 2:end]
target_16_inh_left = target_16_inh_left[:, 2:end]
target_16_exc_right = target_16_exc_right[:, 2:end]
target_16_inh_right = target_16_inh_right[:, 2:end]

target_avg1 = mean(target_16_exc_left, dims=1)[:]
target_avg2 = mean(target_16_inh_left, dims=1)[:]
target_avg3 = mean(target_16_exc_right, dims=1)[:]
target_avg4 = mean(target_16_inh_right, dims=1)[:]

nid1 = target_avg1 .> 0
nid2 = target_avg2 .> 0
nid3 = target_avg3 .> 0
nid4 = target_avg4 .> 0

target_avg1 = target_avg1[nid1]
target_avg2 = target_avg2[nid2]
target_avg3 = target_avg3[nid3]
target_avg4 = target_avg4[nid4]

logs1 = log10.(target_avg1) 
logs2 = log10.(target_avg2)
logs3 = log10.(target_avg3)
logs4 = log10.(target_avg4)

nid12 = logs1 .> -5
nid22 = logs2 .> -5
nid32 = logs3 .> -5
nid42 = logs4 .> -5

nl1 = logs1[nid12]
nl2 = logs2[nid22]
nl3 = logs3[nid32]
nl4 = logs4[nid42]


# mean(10^(log(f) + a)) = 5
# mean(log(f) + a) = log(5) # not correct
# a = log(5) - mean(log(f))

# mean(10^[(log(f) + a) / b]) = 5
# std(10^[(log(f) + a) / b]) = 5

#---- shifting the log normal distribution for the original traces to have mean 5
# a = log10(2) - mean(nl1)
# nl1_shift = nl1 .+ a
# traceShiftedAvg_16_exc_left = 10.0.^(nl1_shift) # shifted & time averaged
# mean(traceShiftedAvg_16_exc_left)
# std(traceShiftedAvg_16_exc_left)
# figure(); hist(log10.(traceShiftedAvg_16_exc_left), bins=100, histtype="step")

#---------- converting raw traces to spike rates -----------#
target_16_exc_left_trim = target_16_exc_left[:,nid1][:,nid12]

# if target_16_exc_left_trim[:,i] contains a negative element, then remove column i
# check how many i's were removed.
n1 = 0
in1 = []
for i in 1:size(target_16_exc_left_trim)[2]
    if any(target_16_exc_left_trim[:, i] .< 0)
        global n1
        global in1
        n1 = n1 + 1 # count number of negative neruons
        append!(in1, i) # add index of negative neuron
    end
end

target_16_exc_left_trim = target_16_exc_left_trim[:, Not(in1)] # remove negative neruons

targetShifted_16_exc_left = target_16_exc_left_trim * 10.0^3.0

target_16_inh_left_trim = target_16_inh_left[:,nid2][:,nid22]
n2 = 0
in2 = []
for i in 1:size(target_16_inh_left_trim)[2]
    if any(target_16_inh_left_trim[:, i] .< 0)
        global n2
        global in2
        n2 = n2 + 1 
        append!(in2, i)    
    end
end

target_16_inh_left_trim = target_16_inh_left_trim[:, Not(in2)] # remove negative neruons

targetShifted_16_inh_left = target_16_inh_left_trim * 10.0^3.0

target_16_exc_right_trim = target_16_exc_right[:,nid3][:,nid32]
n3 = 0
in3 = []
for i in 1:size(target_16_exc_right_trim)[2]
    if any(target_16_exc_right_trim[:, i] .< 0)
        global n3
        global in3
        n3 = n3 + 1
        append!(in3, i)
    end
end

target_16_exc_right_trim = target_16_exc_right_trim[:, Not(in3)] # remove negative neruons

targetShifted_16_exc_right = target_16_exc_right_trim * 10.0^3.0

target_16_inh_right_trim = target_16_inh_right[:,nid4][:,nid42]
n4 = 0
in4 = []
for i in 1:size(target_16_inh_right_trim)[2]
    if any(target_16_inh_right_trim[:, i] .< 0)
        global n4
        global in4
        n4 = n4 + 1
        append!(in4, i)
    end
end

target_16_inh_right_trim = target_16_inh_right_trim[:, Not(in4)] # remove negative neruons

targetShifted_16_inh_right = target_16_inh_right_trim * 10.0^3.0

#----- check the histograms of all four datasets
#targetShifted_16_exc_left = targetShifted_16_exc_left[.!isinf.(targetShifted_16_exc_left)]
#targetShifted_16_exc_left[targetShifted_16_exc_left .== 0] .= 10^-15

#targetShifted_16_inh_left = targetShifted_16_inh_left[.!isinf.(targetShifted_16_inh_left)]
#targetShifted_16_inh_left[targetShifted_16_inh_left .== 0] .= 10^-15

#targetShifted_16_exc_right = targetShifted_16_exc_right[.!isinf.(targetShifted_16_exc_right)]
#targetShifted_16_exc_right[targetShifted_16_exc_right .== 0] .= 10^-15

#targetShifted_16_inh_right = targetShifted_16_inh_right[.!isinf.(targetShifted_16_inh_right)]
#targetShifted_16_inh_right[targetShifted_16_inh_right .== 0] .= 10^-15


figure(); hist(log10.(vec(targetShifted_16_exc_left)), bins=300, histtype="step")
figure(); hist(log10.(vec(targetShifted_16_inh_left)), bins=300, histtype="step")
figure(); hist(log10.(vec(targetShifted_16_exc_right)), bins=300, histtype="step")
figure(); hist(log10.(vec(targetShifted_16_inh_right)), bins=300, histtype="step")




#figure(); hist(log10.(mean(targetShifted_16_exc_left,dims=1)[:]), bins=300, histtype="step")

#------- save targetShifted* as JLD2 files 
using JLD

@save "target16_el.jld" targetShifted_16_exc_left
@save "target16_il.jld" targetShifted_16_inh_left
@save "target16_er.jld" targetShifted_16_exc_right
@save "target16_ir.jld" targetShifted_16_inh_right

# K1 = 5 - mean(nl1)
# K2 = 5 - mean(nl2)
# K3 = 5 - mean(nl3)
# K4 = 5 - mean(nl4)



# target_16_exc_left[target_16_exc_left .< 0] .= 0
# target_16_exc_left_trim = target_16_exc_left[:,nid1][:,nid12] * 10^K1

# target_16_inh_left[target_16_inh_left .< 0] .= 0
# target_16_inh_left_trim = target_16_inh_left[:,nid2][:,nid22] * 10^K2

# target_16_exc_right[target_16_exc_right .< 0] .= 0
# target_16_exc_right_trim = target_16_exc_right[:,nid3][:,nid32] * 10^K3

# target_16_inh_right[target_16_inh_right .< 0] .= 0
# target_16_inh_right_trim = target_16_inh_right[:,nid4][:,nid42] * 10^K4

# trim1 = 5 .+ (target_16_exc_left_trim .- mean(target_16_exc_left_trim))
# trim2 = 5 .+ (target_16_inh_left_trim .- mean(target_16_inh_left_trim))
# trim3 = 5 .+ (target_16_exc_right_trim .- mean())

