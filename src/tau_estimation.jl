function get_estimation_of_tau(A; target_ratio = 10.0)
    n_small = Int(round(A.n/8))
    A_small = A[1:n_small, 1:n_small]
    nnz_A = nnz(A_small)

    function specialized_line_search(target_ratio; f_abs_tol)
        function f(tau)
            LU = ilu(A_small, τ = tau)
            return nnz(LU)/nnz_A
        end
        
        function get_tau_m(A, B, a, b)
            return B - (B-A)*(b/(b-a))
            #return (A + B)/2
        end

        #search for ascending trend
        fact = 1e-1
        
        tau = 1000
        f_A = f(tau)
        while f_A > target_ratio
            tau /= fact
            f_A = f(tau)            
        end

        tau *= fact
        f_B = f(tau)
        while (f_A > f_B) || f_B < target_ratio
            #@show tau, f_A, f_B
            tau *= fact
            f_A = f_B
            f_B = f(tau)
        end

        # bisection
        tau_A = tau/fact
        tau_B = tau

        tau_m = get_tau_m(tau_A, tau_B, (f_A-target_ratio), (f_B - target_ratio))
        f_m = f(tau_m)
        while abs(f_m - target_ratio) > f_abs_tol
            tau_m = get_tau_m(tau_A, tau_B, (f_A-target_ratio), (f_B - target_ratio))
            f_m = f(tau_m)
            #println(tau_A, " ", tau_m, " ", tau_B, "    ", f_A," ",  f_m," ",  f_B)
            
            if f_m < target_ratio
                tau_A = tau_m
                f_A = f_m
            else
                tau_B = tau_m
                f_B = f_m
            end
        end
        return tau_m
    end

    return specialized_line_search(target_ratio, f_abs_tol=0.5)
end