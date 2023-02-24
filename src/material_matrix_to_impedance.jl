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
    @time (header, sp_input, RHS), current_measurement = material_matrix_to_lin_sys(material_matrix, params)          

    sp_input_vals_eval = Array{complex_type}(undef, length(sp_input[3]))          
    b_eval = Array{complex_type}(undef, length(RHS))
    
    Z_list = []
    
    @time begin
      fill_static_entries!(sp_input_vals_eval, sp_input[3])
      fill_static_entries!(b_eval, RHS)
    end 

    for f in f_list
      verbose && @show f        
      w = 2*pi*f
      @time begin
        evaluate_nonstatic_entries_for_w!(sp_input_vals_eval, sp_input[3], w)
        evaluate_nonstatic_entries_for_w!(b_eval, RHS, w)
      end
      
      A_eval = sparse(sp_input[1], sp_input[2], sp_input_vals_eval)
      
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
    (header, sp_input, RHS), current_measurement = material_matrix_to_lin_sys(material_matrix, params)          

    sp_input_vals_eval = Array{complex_type}(undef, length(sp_input[3]))          
    b_eval = Array{complex_type}(undef, length(RHS))
    
    Z_list = []
    
    fill_static_entries!(sp_input_vals_eval, sp_input[3])
    fill_static_entries!(b_eval, RHS)
    
    for f in f_list      
      w = 2*pi*f
      
      evaluate_nonstatic_entries_for_w!(sp_input_vals_eval, sp_input[3], w)
      evaluate_nonstatic_entries_for_w!(b_eval, RHS, w)
   
      A_eval = sparse(sp_input[1], sp_input[2], sp_input_vals_eval)
      
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