using StochasticAiry
using Test
using Random: seed!

@testset "StochasticAiry.jl" begin
    seed!(1)

    n_re_lambda = 1000
    n_im_lambda = 200

    re_grid = collect(range(-8, 15, length=n_re_lambda))
    im_grid = collect(range(-8, 8, length=n_im_lambda))

    lambda_grid = re_grid' .+ im*im_grid;

    ##
    
    t_grid_len = 3*10^3

    t = -20
    s = 10
    t_grid = range(s, t, length=t_grid_len);

    dt = (t - s)/t_grid_len;
    
    
    for beta in [1, Inf]
        phi, phi_dash= rand(SAi(beta), :space, 0, lambda_grid)

        @test size(phi) == size(lambda_grid)
        @test size(phi_dash) == size(lambda_grid)
        
        phi_re, phi_dash_re = rand(SAi(beta), :space, 0, -8:8, brownian_gridlen = 100)

        @test eltype(phi_re) <: Real
        @test eltype(phi_dash_re) <: Real

        sai, sai_dash = rand(SAi(beta), :time, 1 + im, t_grid)

        @test ndims(sai) == 1
        @test ndims(sai_dash) == 1

        sai_re, sai_dash_re = rand(SAi(beta), :time, 0, t_grid)

        @test eltype(sai_re) <: Real
        @test eltype(sai_dash_re) <: Real
    end
    
end
