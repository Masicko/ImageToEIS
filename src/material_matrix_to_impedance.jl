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

function material_matrix_to_impedance(
            material_matrix::Array = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #      
            f_list,            
            #
            complex_type=ComplexF64,
            iterative_solver = "auto",
            verbose = false,
            return_only_linsys = false
            )
  if iterative_solver == "auto"
    if prod(size(material_matrix)) < 1e3
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
      
      @time (header, sp_input, b_eval), current_measurement = material_matrix_to_lin_sys(material_matrix, params, w)
      A_eval = sparse(sp_input...)

      if return_only_linsys
        return A_eval, b_eval, current_measurement
      end
      p = LinearProblem(A_eval,b_eval)
      @time if iterative_solver
        x = LinearSolve.solve(p, KrylovJL_BICGSTAB(); Pl = ILUZero.ilu0(A_eval))
      else                
        x = LinearSolve.solve(p)
      end
      push!(Z_list, 1.0/current_measurement(x, w))
    end    
  else
    params = get_prms_from_pairs(prms_pairs)          

    Z_list = []

    for f in f_list
      verbose && @show f        
      w = 2*pi*f
      
      (header, sp_input, b_eval), current_measurement = material_matrix_to_lin_sys(material_matrix, params, w)
      A_eval = sparse(sp_input...)

      if return_only_linsys
        return A_eval, b_eval, current_measurement
      end

      p = LinearProblem(A_eval,b_eval)
      if iterative_solver       
        x = LinearSolve.solve(p, KrylovJL_BICGSTAB(); Pl = ILUZero.ilu0(A_eval))
      else                
        x = LinearSolve.solve(p)
      end
      push!(Z_list, 1.0/current_measurement(x, w))
    end     
  end

  return f_list, Z_list
end