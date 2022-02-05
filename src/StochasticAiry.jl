module StochasticAiry

using LinearAlgebra
using Random, Distributions
using SpecialFunctions

import Distributions: ContinuousMatrixDistribution
import Random: rand

export SAi, rand

struct SAi <: ContinuousMatrixDistribution
    beta::Real
end

function rand(d::SAi, argtype::Symbol, basepoint::Number, argrange::AbstractRange{T}) where T <: Real

    if argtype == :time
        t_grid = collect(argrange)
        dt = step(argrange)

        U0 = rand(AiryU0(d.beta), argrange)

        return solve_sai(U0, t_grid, dt, basepoint)
    elseif argtype == :space

        return rand(d, :space, basepoint, collect(argrange))
    else
        throw("argtype must be :time or :space")
    end
end

function rand(d::SAi, argtype::Symbol, basepoint::Real, argrange::Array{T,dim};
              verbose=false, brownian_gridlen=3*10^3, brownian_gridmax=10) where {T <: Number, dim}

    if argtype != :space
        throw("Array argranges only defined for argtype == :space")
    end
    
    t = basepoint
    s = brownian_gridmax
    t_grid = range(s, t, length=brownian_gridlen);

    dt = (t - s)/brownian_gridlen

    U0 = rand(AiryU0(d.beta), t_grid)

    solve_phi(U0, collect(t_grid), dt, argrange; verbose=verbose)

end


include("AiryU0.jl")
include("Kernels.jl")
include("AirySolvers.jl")

end
