# StochasticAiry.jl

[![Build Status](https://github.com/damian-t-p/StochasticAiry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/damian-t-p/StochasticAiry.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/damian-t-p/StochasticAiry.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/damian-t-p/StochasticAiry.jl)

A Julia package for sampling the Stochatic Aity function as defined in [1].

# Usage

This package implements the type `SAi(beta)`, which denotes the distribution of the Stochastic Airy function with Dyson parameter `beta`.
When `beta = Inf`, this is the usual deterministic Airy function.

The Stochastic Airy function *SAi<sub>lambda</sub>(t)* is a function of a real-valued time parameter *t* and a complex-valued space parameter *lambda*.
This package allows for sampling of *SAi* and its time derivative *SAi'* with one of these parameters held fixed at a time using the `rand()` method.

## Function of the time parameter

For a fixed *lambda<sub>0</sub>*, we can sample over a range of real *t*-values.
The time grid must be a subtype of `AbstractRange`, as the t-values must be equally spaced.
For example, running
```
t_grid = range(10, -20, length=3*10^3)
sai, sai_dash = rand(SAi(beta), :time, lambda0, t_grid)
```
will produce 
![](/docs/img/SAi-path.png)

## Function of the space parameter

For a fixed *t<sub>0</sub>*, we can sample over an array of complex *lambda*-values.
In the following we run
```
n_re_lambda = 1000
n_im_lambda = 200

re_grid = collect(range(-8, 6, length=n_re_lambda))
im_grid = collect(range(-8, 8, length=n_im_lambda))

lambda_grid = re_grid' .+ im*im_grid;

sai, sai_dash = rand(SAi(beta), :space, y0, lambda_grid)
```

The following plot then depticts `log.(abs.(sai))`:
![](/docs/img/SAi-lambda-log-abs.png)

# References
[1] Lambert & Paquette, "Strong approximation of Gaussian beta-ensemble characteristic polynomials: the edge regime and the stochastic Airy function," arXiv:2009.05003 [math.PR] July 2021
