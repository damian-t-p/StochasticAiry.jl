function solve_sai(U0, t_grid, dt, lambda::Number)
    
    grid_len = length(t_grid)
    s = t_grid[1]
    
    t_outer = t_grid .- t_grid'
    
    K0 = make_volterra_kernel(U0)
    pert = make_volterra_kernel(t_outer)
    
    #c1 = airyai(s) * exp(-lambda * sqrt(s))
    #c2 = -sqrt(s) * c1
    
    c1 = airyai(s + lambda)
    c2 = airyaiprime(s + lambda)
    
    v1 = ones(grid_len)
    v2 = U0[:,1]
    v3 = t_outer[:,1]
    
    v = (c2 * v1) + (c1 * v2) + (c1 * lambda * v3)
    
    A = I - dt/2 * K0
    B = -dt/2 * pert
    
    sai_path = (A + lambda * B) \ v
    
    trapz_matrix = make_volterra_kernel(ones(size(U0))) * dt/2
    
    return (
        SAi = trapz_matrix * sai_path,
        DSAi = sai_path
    )
end


function solve_phi(U0, t_grid, dt, lambda_grid; verbose=false)
    grid_len = length(t_grid)
    
    s = t_grid[1]
    
    t_outer = t_grid .- t_grid'
    
    K0 = make_volterra_kernel(U0)
    pert = make_volterra_kernel(t_outer)
    
    u_deriv = [zeros(grid_len - 1); 1]
    u_full = dt/2 * ([0; ones(grid_len - 1)] + [ones(grid_len - 1); 0]);
    
    u = [u_deriv u_full]
    
    v1 = ones(grid_len)
    v2 = U0[:,1]
    v3 = t_outer[:,1]

    A = I - dt/2 * K0
    B = -dt/2 * pert

    BA_inv = B / A
    uA_inv = A' \ u

    N = BA_inv

    Nv_curr = [v1 v2 v3]
    
    # the 2*3 matrix u^t A^-1 v for u1-u2, v1-v3 for each power
    coefs_v = Array{Complex{Float64}, 3}(undef, (grid_len, 2, 3))
    
    for pow in 1:grid_len
        
        coefs_v[pow,:,:] = uA_inv' * Nv_curr * (-1)^(pow + 1)
        Nv_curr = N * Nv_curr
        
     end

    log_coefs_v = log.(coefs_v)
    
    phi_deriv_grid = undef_copy(lambda_grid)
    phi_full_grid = undef_copy(lambda_grid)
    
    for (i, lambda) in pairs(lambda_grid)

        c1 = airyai(s + lambda)
        c2 = airyaiprime(s + lambda)
        
        log_lambda = log(lambda)
        
        # Make the array [0, log(lambda), 2 log(lambda), ...]
        log_lambda_pows = accumulate(+, [0; log_lambda * ones(grid_len - 1)])
        
        # Make the array [C_0(lambda) lambda^0, C_1(lambda) lambda^1, ...]
        sums_deriv = sum(exp.(log_lambda_pows .+ log_coefs_v[:,1,:]), dims = 1)
        sums_full = sum(exp.(log_lambda_pows .+ log_coefs_v[:,2,:]), dims = 1)
        
        if any(isnan.(sums_deriv)) | any(isnan.(sums_full))
            ss = log_lambda_pows .+ log_coefs_v + maxval_log_cum
            display(ss)
            display(cumsum(ss, dims = 1))
            display(sums)
            throw("NAN sums")
        end
        
        phi_deriv_grid[i] = (c2 * sums_deriv[1,1] + c1 * sums_deriv[1,2] + c1 * lambda * sums_deriv[1,3]) + c1
        phi_full_grid[i] = (c2 * sums_full[1,1] + c1 * sums_full[1,2] + c1 * lambda * sums_full[1,3]) + c1
        
        if isnan(phi_deriv_grid[i]) | isnan(phi_full_grid[i])
            println(i)
            throw("NAN phi")
        end
    end
    
    if verbose
        return (phi_full_grid, phi_deriv_grid, coefs_v, log_coefs_v, maxval_logs)
    else
        return (
        SAi = phi_full_grid,
        DSAi = phi_deriv_grid)
    end
end;
