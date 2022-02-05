export AiryU0

struct AiryU0 <: ContinuousMatrixDistribution
    beta::Real
end

function rand(d::AiryU0, trange::AbstractRange{T}) where T <: Real
    t_grid = collect(trange)
    dt = step(trange)

    B = rand_brownian(t_grid, dt)

    return make_U0(t_grid, dt, B, d.beta)
end

##

function rand_brownian(t_grid, dt)
    grid_len = length(t_grid)
    
    B_diff = randn(grid_len) * sqrt(abs(dt))
    cumsum(B_diff)
end

function make_U0(t_grid, dt, B, beta)
    U0_row = t_grid.^2/2 + B * sqrt(4/beta)
    U0_row .- U0_row'
end

