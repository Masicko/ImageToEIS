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
            iterative_solver = false,
            verbose = false,
            return_only_linsys = false
            )
           
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
      @time if iterative_solver       
        x = bicgstabl(A_eval, b_eval)
      else                
        x = A_eval \ b_eval
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
      if iterative_solver       
        x = bicgstabl(A_eval, b_eval)
      else                
        x = A_eval \ b_eval
      end
      push!(Z_list, 1.0/current_measurement(x, w))
    end     
  end

  return f_list, Z_list
end