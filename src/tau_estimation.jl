function get_estimation_of_tau(A; target_ratio = 10.0)
    bandwidth = A.rowval[4]

    n_small = min(A.n, bandwidth*10)
    A_small = A[1:n_small, 1:n_small]
    nnz_A = nnz(A_small)

    #n_tiny = min(A.n, bandwidth*5)
    #A_tiny = A[1:n_tiny, 1:n_tiny]    
    
    function to_opt_temp(A_act)
        return tau -> begin
                if tau[1] < 1e-14
                    return 1000000
                else
                    LU = ilu(A_act, τ = tau[1])
                    nnz_ratio = nnz(LU)/nnz_A
                    return abs(nnz_ratio - target_ratio)
                end    
            end
    end

    #initial_tau = optimize(to_opt_temp(A_tiny), [1.0], 
    #            NelderMead(),
    #            Optim.Options(
    #                x_tol=1.0
    #            )
    #    ).minimizer[1]
    initial_tau = 1.0
    res = optimize(to_opt_temp(A_small), [initial_tau], 
                   NelderMead(),
                   Optim.Options(
                       x_tol=1.0
                   )
           ).minimizer[1]
    #@show res
    return res
end