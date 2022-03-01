function subs_parameters_to_matrix!(A, p)
  m,n = size(A)
  for i in 1:m
    for j in 1:n
      if typeof(A[i,j]) == String
        for name in fieldnames(parameters)    
          #@show "p."*string(name), A[i,j]
          A[i,j] = replace(A[i,j], "p."*string(name) => getfield(p, name))
        end
      end
      #A[i,j] = simplify!(Sym(A[i,j]))
    end
  end
end


function evaluate_matrix_for_w!(A_eval, A, w)
  m,n = size(A)
  for i in 1:m
    for j in 1:n
      if typeof(A[i,j]) == String       
        A_eval[i,j] = eval(Meta.parse(
            replace(A[i,j], "w" => w)
          )
        )        
      else
        A_eval[i,j] = A[i,j]
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



function image_to_EIS(
            matrix_input::Union{Array, String} = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #            
            f_range="TPE",
            #f_range=[2.0^i for i in -3:10], 
            #            
            complex_type=ComplexF64,
            iterative_solver = false,
            verbose = false,
            pyplot=true)
           
  if typeof(matrix_input) == String
    material_matrix = file_to_matrix(matrix_input)
  else
    material_matrix = matrix_input
  end
  
  if size(material_matrix)[1] < 1 || size(material_matrix)[2] < 1
    println("ERROR: wrong input matrix")
    return
  end
  
  two_point_extrapolation = false
  if typeof(f_range)==String && (f_range == "two_point_extrapolation" || f_range == "TPE")
    two_point_extrapolation = true
    f_range = [0.1, 10000]
  end
  
  params = get_prms_from_pairs(prms_pairs)         
  
  header, A, b = material_matrix_to_equations(material_matrix, params)        
    
  begin
    subs_parameters_to_matrix!(A, params)
    
    Z_range = []
    A_eval = Matrix{complex_type}(undef, size(A)...)
    b = convert.(complex_type, b)
    for f in f_range
      verbose && @show f
      # critical line !!!
      evaluate_matrix_for_w!(A_eval, A, 2*pi*f) 
      
      #return A_eval, b
      if iterative_solver       
        x = bicgstabl(A_eval, b)
      else                
        x = A_eval \ b
      end      
      push!(Z_range, 1/x[1])
    end
  end
  
  
  
  
  if two_point_extrapolation
    omegas = 2*pi .* [f_range[1], f_range[2]]
    Zs = [Z_range[1], Z_range[2]]
    R_ohm, R, C = get_R_RC_prms_from_2Z(omegas, Zs)
    
    f_range = [2.0^i for i in -5: 0.8 : 17]
    Z_range = get_Z_from_R_RC(f_range, R_ohm, R, C)    
  end
    
  pyplot && nyquistPlot(Z_range)
  
  return f_range, Z_range  
end
