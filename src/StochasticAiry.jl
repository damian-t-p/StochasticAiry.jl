module StochasticAiry

using LinearAlgebra
using Random, Distributions
using SpecialFunctions

import Distributions: ContinuousMatrixDistribution
import Random: rand

export rand

include("AiryU0.jl")
include("Kernels.jl")
include("SAi.jl")

end
