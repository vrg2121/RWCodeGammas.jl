module DrawGammas
export StructParams, StructAllParams

import Random: rand, Random
import JLD2: @save
import Distributions: Normal

wd = pwd()
R = wd*"/Results" # path to results

D = wd*"/Data" # path to data

G = wd*"/Guesses" # path to guesses

# Set Up Parameters
include("./functions/ParamsFunctions.jl") 
include("./Params.jl")

import .ParamsFunctions: StructParams
import .Params: ModelParams, StructAllParams

println("Setting up parameters...")

P = ModelParams(D, G);

# overwrite the hardcoded gammas with draws from the random distribution
function gen_params_normal(P::StructAllParams, n::Int)
    gSolar = Vector{Float64}(undef, n)
    gWind = Vector{Float64}(undef, n)
    gBatteries = Vector{Float64}(undef, n)

    Random.seed!(123)

    gSolar .= rand(Normal(0.35, 0.0001), n) 
    gWind .= rand(Normal(0.2, 0.0001), n)
    gBatteries .= rand(Normal(0.35, 0.0001), n)

    wd = pwd()

    for i = 1:n
        P.params.gammaS = gSolar[i]
        P.params.gammaW = gWind[i]
        P.params.gammaB = gBatteries[i]

        @save "$wd/Data/Params/P_$i.jld2" P
    end

    println("All gammas have been drawn from the Normal and saved in ./Data/Params")
end

gen_params_normal(P, 3)

end # module DrawGammas
