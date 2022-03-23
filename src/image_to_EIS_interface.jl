
function get_material_matrix(matrix_input)
  input_path = ""
  if typeof(matrix_input) == String
    input_path = matrix_input
    material_matrix = file_to_matrix(matrix_input)
  else
    material_matrix = matrix_input
  end
  if size(material_matrix)[1] < 1 || size(material_matrix)[2] < 1
    println("ERROR: wrong input matrix")
    return
  end
  return input_path, material_matrix
end

function change_extension_to(path, new_ext)
  last_dot_position = findall(x -> x == '.', path)[end]
  if last_dot_position < 2
    println("ERROR: input_path: last_dot_position = $(last_dot_position) < 2")
  end
  return path[1:last_dot_position-1]*"."*new_ext
end

function image_to_EIS(
            matrix_input::Union{Array, String} = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #            
            f_list="TPE",
            TPE_f_list = [2.0^n for n in -5 : 0.5 : 15],
            TPE_warning = true,
            #f_list=[2.0^i for i in -3:10], 
            #            
            complex_type=ComplexF64,
            iterative_solver = false,
            verbose = false,
            pyplot=true,
            return_R_RC=false,
            export_z_file="",
            #export_z_file="!use_file_name",
            save_also_image="",
            #save_also_image="!input",
            
            
            keep_constant_physical_height = true,
            store_R_RC=""
            )
  
  input_path, material_matrix = get_material_matrix(matrix_input)
  
  two_point_extrapolation = false
  if typeof(f_list)==String && (f_list == "two_point_extrapolation" || f_list == "TPE")
    two_point_extrapolation = true
    f_list = [0.1, 10000]
  end
  
  
  
  
  
  extract_R_RC = false  
  if return_R_RC || two_point_extrapolation || store_R_RC != ""
    extract_R_RC = true
  end
                    
  f_list, Z_list =  matrix_to_impedance(
    material_matrix,         
    prms_pairs,
    #      
    f_list = f_list,            
    #
    complex_type=complex_type,
    iterative_solver = iterative_solver,
    verbose = verbose
  )
  
  if keep_constant_physical_height
    #normalize output
    Z_list .*= size(material_matrix)[1] / size(material_matrix)[2]
  end
  
  
  if extract_R_RC 
    R_ohm, R, C = (1.0, 1.0, 1.0)
    if length(f_list) < 2
      println("ERROR: extract R_ohm, R, C: length(f_list) < 2, $(length(f_list)) < 2")
      return
    else
      omegas = 2*pi .* [f_list[1], f_list[end]]
      Zs = [Z_list[1], Z_list[end]]
      try
        R_ohm, R, C = get_R_RC_prms_from_2Z(omegas, Zs)      
      catch        
        println("ERROR: R_RC circuit evaluation failed for (Z1, Z2) = $(Zs)  ... returning (f_list, Z_list)\n")
        extract_R_RC = false
      end      
    end
  end
  
  if two_point_extrapolation && extract_R_RC    
    f_list = TPE_f_list
    Z_list = get_Z_from_R_RC(f_list, R_ohm, R, C)
    TPE_warning && println("WARNING: TPE (two_point_extrapolation) has been used!")
  end

  pyplot && nyquistPlot(Z_list)
  
  if store_R_RC != "" && extract_R_RC
    open(store_R_RC, "a") do f
      write(f, @sprintf("%s\t%s\t%E\t%E\t%E\n", 
        Dates.now(), 
        (input_path == "" ? "<matrix_input>" : input_path ), 
        R_ohm, R, C)
      )
    end
  end
  
  if export_z_file == "!use_file_name"
    if input_path != ""      
      export_EIS_to_Z_file(change_extension_to(input_path, "z"), f_list, Z_list)
    else
      println("ERROR: cannot export_z_file: input_path == \"\"!")
    end
  elseif export_z_file != ""
    
    export_EIS_to_Z_file(export_z_file, f_list, Z_list)
        
    if save_also_image == "!input" && input_path != ""
      cp(
        input_path, 
        change_extension_to(export_z_file, input_path[end-2 : end]),
        force=true
      )
    elseif save_also_image != ""      
      matrix_to_file(save_also_image, material_matrix)
    end
  end
  
  if extract_R_RC && return_R_RC
    return R_ohm, R, C
  else
    return f_list, Z_list  
  end
end
