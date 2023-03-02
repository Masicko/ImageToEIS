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

#################### 2D testing
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

