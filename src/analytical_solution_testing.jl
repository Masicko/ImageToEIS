function analytic_eval_1Dmatrix(matrix, prms)
    if size(matrix)[1] > 1
        println("ERROR: not a 1D domain !!!")
        return
    end
    aux_matrix = Int16.(zeros(length(matrix) + 3))
    aux_matrix .= 1
    aux_matrix[2:end-2] = matrix
    aux_matrix[end] = 6
    
    
    LSM_count = -1
    YSZ_count = 0
    YSZ_LSM_interfaces = 0
    for i in 1:(length(aux_matrix)-1)
        node_id = aux_matrix[i]
        next_id = aux_matrix[i+1]
        if node_id == 0
            YSZ_count += 1
            if next_id == 1
                YSZ_LSM_interfaces += 1
            end
        elseif node_id == 1
            LSM_count += 1
            if next_id == 0
                YSZ_LSM_interfaces += 1
            end
        end
        
    end
    
    R_YSZ = Dict(prms)["R_YSZ"]
    R_LSM = Dict(prms)["R_LSM"]
    C_p = Dict(prms)["C_pol_LSM"]
    R_p = Dict(prms)["R_pol_LSM"]
    
    
    R_HF = R_YSZ*YSZ_count + R_LSM*LSM_count
    R_LF = R_HF + YSZ_LSM_interfaces*Dict(prms)["R_pol_LSM"]
    
    n = YSZ_LSM_interfaces
    Z = w -> R_HF + (n*R_p)/(1 + im*w*R_p*C_p)
    
    
    return R_HF, R_LF, Z
end

function plot_analytic(matrix, prms, f_list)
    HF, LF, Z = analytic_eval_1Dmatrix(matrix, prms)
    ImageToEIS.nyquistPlot(Z.( f_list .* 2*pi), label= "anal")
     return Z.( f_list .* 2*pi)
 end

function compare_mat(matrix, prms, f_list=f_list)
    Z_anal = plot_analytic(matrix, prms, f_list)
    
    f, Z = image_to_EIS(matrix, #generate_matrix((6,6), 0.0, 0.0), 
           prms,
           #f_list=[0.0],
           f_list=f_list, 
           iterative_solver=false,
           #return_R_RC=true
           return_specific_impedance=false
           )
    return println("Error: ", sum(abs.(Z_anal .- Z) ./ abs.(Z)))
end



function test_1D()
    MM1D = [0 1 1 0 0 1 0 1 0 0 0 0 0 1 1 0 0 0 1 1 0 1 1 1 0 0 0 1 1 0]
    PRMS = ["R_YSZ" => 100.0, "R_LSM" => 0.1, "R_pol_LSM" => 3.0, "C_pol_LSM" => 0.005]
    f_list=[10^n for n in -4:0.5:4]

    compare_mat(MM1D, PRMS, f_list)
end













############ SymPy ########

function sym_get_A_b_from_matrix(matrix, prms, f)
    @syms w

    @syms R_Y, R_L, R_p, C_p



    RC(R,C,w) = R/(1+ im*R*C*w)

    my_RC = RC(R_p, C_p, w)
    
    
    function sym_get_Z_from_conscutive_nodes(n1, n2)
        if      n1 == 1
                if      n2 == 1
                  return (R_L)
                elseif  n2 == 0  
                  return (R_L/2 + R_Y/2 + RC(R_p, C_p, w))           
                end
        elseif  n1 == 0
                if      n2 == 1
                  return (R_L/2 + R_Y/2 + RC(R_p, C_p, w))            
                elseif  n2 == 0
                  return (R_Y)
                end
        else
            println("ERROR: get_Y_entry...")
        end        
    end
    
    function sym_make_Z_from_matrix(matrix)
        if size(matrix)[1] > 1
            println("ERROR: not a 1D domain !!!")
            return
        end
        aux_matrix = Int16.(zeros(length(matrix) + 3))
        aux_matrix .= 1
        aux_matrix[2:end-2] = matrix
        aux_matrix[end] = 6
    
        Zs = []
        for i in 1:length(aux_matrix)-2
            push!(Zs, 
                sym_get_Z_from_conscutive_nodes(aux_matrix[i], aux_matrix[i+1])
            )
        end
        return Zs
    end
    
    g_row(Zl, Zh, i, Us) = Us[i]*(-1/Zl - 1/Zh) + Us[i-1]/Zl + Us[i+1]/Zh
    
    function my_solve(Zs, sym_prms)
        num_nodes = length(Zs)+1
        Us = []
        for i in 1:num_nodes
            push!(Us, symbols("U$(i)"))
        end
        
        rows = []
        eqs =  []
        for i in 2:num_nodes-1
            push!(rows, 
                g_row(Zs[i-1], Zs[i], i, Us)
            )
            push!(eqs, 
                rows[i-1](sym_prms..., Us[1] => 1.0, Us[num_nodes] => 0.0)
            )
        end
        #@show Tuple(Us[2:num_nodes-1])
        #@show [eqs[i] for i in 1:length(eqs)]
        res = linsolve([eqs[i] for i in 1:length(eqs)], Tuple(Us[2:num_nodes-1]))
        #res = nothing
        return rows, eqs, res, Us
    end
    
    function sym_get_A_b(eqs, Us)
        num_unknowns = length(eqs)
        A = Matrix{Any}(undef, (num_unknowns, num_unknowns))
        b = Vector{Any}(undef, num_unknowns)
        for row in 1:num_unknowns
            for col in 1:num_unknowns
                A[row, col] = diff(eqs[row], Us[col+1])
            end
            b[row] = -eqs[row]([U .=> 0.0 for U in Us]...)
        end
        return A, b
    end
    
    function sym_get_sym_prms(prms)
        return [    
            R_Y => Dict(prms)["R_YSZ"],
            R_L => Dict(prms)["R_LSM"],
            R_p => Dict(prms)["R_pol_LSM"],
            C_p => Dict(prms)["C_pol_LSM"]
        ]
    end
    
    function sym_eval_A_b_at_w(A, b, input_w)
        @syms w
        B = deepcopy(A)
        for i in 1:length(A[:])
            B[i] = A[i](w => input_w)
        end
        c = deepcopy(b)
        for i in 1:length(c[:])
            c[i] = b[i](w => input_w)
        end
        return B, c
    end

    # function get_meas_I(res, Zs, sym_prms)
    #     (input_w) -> (1.0 - res[1](w => input_w))/Zs[1](sym_prms..., w=>input_w)
    # end

    sym_prms = sym_get_sym_prms(prms)
    Zs = sym_make_Z_from_matrix(matrix)
    rows, eqs, res, Us = my_solve(Zs, sym_prms)
    #return rows, eqs, res, Us
    A, b = sym_get_A_b(eqs, Us)
    #meas_I = get_meas_I(res, Zs, sym_prms)
    
    return Zs, rows, eqs, A, b, res, Us
    #return sym_eval_A_b_at_w(A, b, f*2*pi)

end

#ll = lambdify(res[1])

















#################### 2D testing
function get_2D_testing_matrix()
    MM2D = [       0 1 2 0 1 2 2 2 2 2 1 0; 
                    2 2 2 2 1 1 0 2 2 0 1 2; 
                    1 2 0 2 1 0 1 0 1 1 2 2; 
                    0 2 1 2 0 2 2 2 0 2 0 0; 
                    0 2 1 1 1 1 2 1 0 0 2 1; 
                    2 1 2 1 1 0 2 0 1 0 2 1; 
                    1 1 0 1 2 2 2 0 0 2 2 2; 
                    2 2 0 2 2 2 2 2 2 2 0 2; 
                    1 2 2 2 1 0 1 2 1 0 0 2; 
                    1 0 2 2 2 2 1 2 1 0 2 0; 
                    1 2 0 0 0 1 2 0 0 0 2 0; 
                    0 2 0 2 1 1 2 0 2 2 1 2]
end
