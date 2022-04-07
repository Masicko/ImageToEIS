
function get_I_name(x,y)
    return "I_$(x)_$(y)"
end

function get_Z_name(x,y)
    return "Z_$(x)_$(y)"
end



# auxilary_A has 2 on left, m*n+3 on the right, -1 on top and bottom
function dict_to_lin_sys(Z_dict, auxilary_A)
  matrix_header = String[]
  
  equation_rows = []
  
  sparse_input_row_idxs = Int64[]
  sparse_input_col_idxs = Int64[]
  sparse_input_vals = Union{Float64, Function}[]
  sparse_input = (sparse_input_row_idxs, sparse_input_col_idxs, sparse_input_vals)
  sys_row_idx = 0
  
  RHS = Float64[]
  
  
  m,n = size(auxilary_A) .-(2,2)
  
  function find_idx(token, array)
      res = findall(x->x==token, array)
      if length(res)>0
          return res[1]
      else
          return -1
      end
  end
  
  function add_value_to_matrix_row(matrix_header, equation_row, I_name, value)
      idx = find_idx(I_name, matrix_header)
      if idx == -1
        idx = length(matrix_header) + 1
        push!(matrix_header, I_name)
      end
      push!(equation_row, (idx, value))
  end
  
  # TODO optimize
  function previous(i, j)    
    if j == 0
      return [1]
    else
      res = [auxilary_A[i+1, j-1+1], auxilary_A[i-1+1, j+1]]
      deleteat!(res, findall(x->x==-1,res))
      return res
    end
  end
  
  function next(i,j)
    if j == 0      
      return auxilary_A[2:end-1, 2]
    else
      res = [auxilary_A[i+1+1, j+1], auxilary_A[i+1, j+1+1]]
      deleteat!(res, findall(x->x==-1,res))
      return res
    end
  end
  
  function add_row_to_sparse_input!(sparse_input,
                                    new_row, 
                                    row_idx)
    
    
    for item in new_row
      push!(sparse_input[1], row_idx)
      push!(sparse_input[2], item[1])
      push!(sparse_input[3], typeof(item[2]) <: Number ? float(item[2]) : item[2])
    end
  end
  
  # I vertices
  function add_equation_I_row(i, j)
    new_row = []
    act_id = auxilary_A[i+1, j+1]    
    for node_id in previous(i,j)        
        I_name = get_I_name(node_id, act_id)
        add_value_to_matrix_row(matrix_header, new_row, I_name, -1)
    end
    for node_id in next(i,j)
        I_name = get_I_name(act_id, node_id)
        add_value_to_matrix_row(matrix_header, new_row, I_name, 1)
    end

    sort!(new_row, by = first)
    
    sys_row_idx += 1
    add_row_to_sparse_input!(sparse_input, new_row, sys_row_idx)
    
    push!(RHS, 0.0) 
    return
  end

  # initial vertex "2"
  add_equation_I_row(1,0)
  
  # other vertices
  for i in 1:m, j in 1:n
      add_equation_I_row(i,j)
  end
      
    
  # paths
  paths_storage = []

  function get_sufficient_paths(auxilary_A)
    A = auxilary_A[2:end-1, 2:end-1]
    for row in 1:m-1      
      for fall_point in n+1:-1:1        
        new_path = [2]            
        column = 1
        act_row = row
        while column < n + 1
          push!(new_path, A[act_row, column])
          if fall_point == column && act_row == row
            act_row += 1
          else
            column += 1
          end          
        end
        push!(new_path, m*n + 3)
        push!(paths_storage, deepcopy(new_path))
      end
    end
    push!(paths_storage, deepcopy(auxilary_A[m+1, :]))
    return paths_storage
  end
    
  # U equations
  for path in get_sufficient_paths(auxilary_A)
      new_row = []
      
      for (i, idx1) in enumerate(path[1:end-1])
          idx2 = path[i+1]
          add_value_to_matrix_row(
              matrix_header, 
              new_row, 
              get_I_name(idx1, idx2),
              Z_dict[get_Z_name(idx1, idx2)]
          )
      end
    
      sort!(new_row, by = first)
    
      sys_row_idx += 1
      add_row_to_sparse_input!(sparse_input, new_row, sys_row_idx)
      
      push!(RHS, 1.0)
  end
  
  return matrix_header, sparse_input, RHS
  
#   sys_A = Matrix(undef, length(equation_rows), length(matrix_header)-1)
#   sys_A .= 0
#   sys_b = Vector(undef, length(equation_rows))
#   sys_b .= 0
#   for i in 1:size(sys_A)[1]
#       sys_A[i,1:length(equation_rows[i])-1] = equation_rows[i][2:end]
#       sys_b[i] = equation_rows[i][1]
#   end
#   return matrix_header, sys_A, sys_b
end




function get_auxilary_graph_matrix(m, n)
    a = Vector(undef, m*n)
    for i in 1:length(a)
        a[i] = i+2
    end
    pre_A = reshape(a, m, n)
    A = Matrix(undef, m+2, n+2)
    A[2:end-1, 2:end-1] = pre_A
    A[:, 1] .= 2
    A[:, end] .= m*n+3
    A[1, :] .= -1
    A[end, :] .= -1
    return A
end


function get_large_material_matrix(pre_A)
    (m,n) = size(pre_A)
    A = Matrix(undef, m+2, n+2)
    A[2:end-1, 2:end-1] = pre_A
    A[:, 1] .= i_LSM
    A[:, end] .= i_LSM
    A[1, :] .= -1
    A[end, :] .= -1    
    return A
end



function material_matrix_to_lin_sys(
            material_A::Matrix=[0 0 0; 1 1 1], 
            p::parameters=parameters()
    )
    (m,n) = size(material_A)
    aux_A = get_auxilary_graph_matrix(m,n)
    large_material_A = get_large_material_matrix(material_A)
    
    function add_Z_to_dict(i,j, k, l)
        n1 = large_material_A[i,j]
        n2 = large_material_A[i+k, j+l]
        if n2 == -1
            return
        end
        
        idx1 = aux_A[i, j]
        idx2 = aux_A[i+k, j+l]
        if idx2 < idx1
            aux = idx1
            idx1 = idx2
            idx2 = aux
        end
        current_key = get_Z_name(idx1, idx2)
        if haskey(Z_dict, current_key)            
            f = deepcopy(Z_dict[current_key])            
            Z_dict[current_key] = w -> f(w) +
                get_Z_entry_from_material_matrix_codes(n1, n2, p)(w)
        else            
            Z_dict[current_key] = 
                get_Z_entry_from_material_matrix_codes(n1, n2, p)
        end        
    end
    
    Z_dict = Dict()
    for i in 2:m+1
        for j in 2:n+1
            add_Z_to_dict(i,j, +1, 0) 
            add_Z_to_dict(i,j, -1, 0)
            add_Z_to_dict(i,j, 0, +1)
            add_Z_to_dict(i,j, 0, -1)
        end
    end
    for i in 2:m+1
        add_Z_to_dict(i,1, 0, +1)
        add_Z_to_dict(i,n+2, 0, -1)
    end    
    return dict_to_lin_sys(Z_dict, aux_A)    
end
