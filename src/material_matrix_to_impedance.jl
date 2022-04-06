function evaluate_matrix_for_w!(A_eval, A, w)
  m,n = size(A)
  for i in 1:m
    for j in 1:n
      if typeof(A[i,j]) <: Function
        A_eval[i,j] = A[i,j](w)
      else
        #A_eval[i,j] = A[i,j]
      end
    end
  end 
end


function get_prms_from_pairs(pairs)
  p = parameters()
  for pair in pairs
    setfield!(p, Symbol(pair[1]), float(pair[2]))
  end
  return p
end



function material_matrix_to_impedance(
            material_matrix::Array = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #      
            f_list,            
            #
            complex_type=ComplexF64,
            iterative_solver = false,
            verbose = false
            )
           
  if verbose
    @time params = get_prms_from_pairs(prms_pairs)
    @time header, A, b = material_matrix_to_lin_sys(material_matrix, params)          
  
    Z_list = []
    A_eval = Matrix{complex_type}(undef, size(A)...)
    
    @time for i in 1:length(A[:])
      if !(typeof(A[i]) <: Function)
        
        A_eval[i] = A[i]      
      end  
    end
    
    
    b = convert.(complex_type, b)
    for f in f_list
      verbose && @show f        
      @time evaluate_matrix_for_w!(A_eval, A, 2*pi*f) 
      
      #return A_eval, b
      @time if iterative_solver       
        x = bicgstabl(A_eval, b)
      else                
        x = A_eval \ b
      end      
      push!(Z_list, 1/x[1])
    end    
  else
    params = get_prms_from_pairs(prms_pairs)         
    header, A, b = material_matrix_to_lin_sys(material_matrix, params)          
  
    Z_list = []
    A_eval = Matrix{complex_type}(undef, size(A)...)
        
    for i in 1:length(A[:])
      if !(typeof(A[i]) <: Function)
        
        A_eval[i] = A[i]      
      end  
    end
            
    b = convert.(complex_type, b)
    for f in f_list
      verbose && @show f        
      evaluate_matrix_for_w!(A_eval, A, 2*pi*f) 
      
      #return A_eval, b
      if iterative_solver       
        x = bicgstabl(A_eval, b)
      else                
        x = A_eval \ b
      end      
      push!(Z_list, 1/x[1])
    end
  end
    
  

  
  
  return f_list, Z_list      
end
