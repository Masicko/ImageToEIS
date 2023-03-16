function get_prms_from_pairs(pairs)
  p = parameters()
  for pair in pairs
    setfield!(p, Symbol(pair[1]), float(pair[2]))
  end
  return p
end

function fill_static_entries!(list_eval, list)
  for i in eachindex(list)
    if !(typeof(list[i]) <: Function)
      list_eval[i] = list[i]      
    end
  end
end

function evaluate_nonstatic_entries_for_w!(eval, template, w)
  for i in eachindex(template)   
    if typeof(template[i]) <: Function
      eval[i] = template[i](w)
    end
  end
end

function set_initial_values(dims)
  x = Array{Float64}(undef, prod(dims))
  
  for k in 1:dims[3]
      for j in 1:dims[2]
          start = (k-1)*dims[1]*dims[2] + (j-1)*dims[1] + 1
          ending = start + dims[1] - 1
          x[start : ending] .= 1.0 - j/(dims[2]+1)
          #@show start, ending
      end
  end
  return x
end

function try_convert_system_to_real(sp_input, RHS)
  sp_input_values = []
  b_eval = []
  try
    sp_input_values = Float64.(sp_input[3])
    b_eval = Float64.(RHS)
  catch
    sp_input_values = sp_input[3]
    b_eval = RHS
  end
  return sp_input_values, b_eval
end

function material_matrix_to_impedance(
            material_matrix::Array = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #      
            f_list,
            #
            iterative_solver = "auto",
            verbose = false,
            return_only_linsys = false,
            tau = "auto",
            fill_in_ratio = 10
            )
  if iterative_solver == "auto"
    if prod(size(material_matrix)) < 15^3
      iterative_solver = false
    else
      iterative_solver = true
    end
  end

  if verbose
    @time params = get_prms_from_pairs(prms_pairs)          

    Z_list = []

    for f in f_list
      verbose && @show f        
      w = 2*pi*f

      print("-- build lin sys:      "); 
      @time (header, sp_input, RHS), current_measurement = material_matrix_to_lin_sys(material_matrix, params, w)
      sp_input_values, b_eval = try_convert_system_to_real(sp_input, RHS)     
      A_eval = sparse(sp_input[1], sp_input[2], sp_input_values)


      if return_only_linsys
        return A_eval, b_eval, current_measurement
      end
      p = LinearProblem(A_eval,b_eval)
      if iterative_solver
        print("-- estiamtion of tau:   "); 
        @time tau == "auto" && (tau = get_estimation_of_tau(A_eval; target_ratio = fill_in_ratio))

        print("-- solver --- ilu:      "); @time LU = ilu(A_eval, τ = tau)
        print("      |          \\--------------- fill_in_ratio: ", nnz(LU)/nnz(A_eval), "\n")
        print("      \\------ bicgstab: "); @time x = bicgstabl(A_eval, b_eval, 2, Pl = LU
          ,max_mv_products = 2000
          )
      else                
        print("-- solver ---- LU:      "); @time x = LinearSolve.solve(p)
      end
      push!(Z_list, 1.0/current_measurement(x, w))
    end    
  else
    params = get_prms_from_pairs(prms_pairs)          

    Z_list = []

    for f in f_list
      verbose && @show f        
      w = 2*pi*f
      
      (header, sp_input, RHS), current_measurement = material_matrix_to_lin_sys(material_matrix, params, w)
      sp_input_values, b_eval = try_convert_system_to_real(sp_input, RHS)     
      A_eval = sparse(sp_input[1], sp_input[2], sp_input_values)

      if return_only_linsys
        return A_eval, b_eval, current_measurement
      end

      p = LinearProblem(A_eval,b_eval)
      if iterative_solver
        tau == "auto" && (tau = get_estimation_of_tau(A_eval; target_ratio = fill_in_ratio))
        LU = ilu(A_eval, τ = tau)
        x = bicgstabl(A_eval, b_eval, 2, Pl = LU
          ,max_mv_products = 2000
          )       
      else                
        x = LinearSolve.solve(p)
      end
      push!(Z_list, 1.0/current_measurement(x, w))
    end     
  end

  return f_list, Z_list
end