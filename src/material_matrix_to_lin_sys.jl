# idx in A_aux matrix
# assuming that x, y can be in [1, 2, ... m*n+3]
function aux_matrix_connectivity_entries_to_lin_idx(x,y, m,n)
    if (x, y) == (1, 2)
      return 1
    elseif x == 2
      # y - 2 + 1
      return y - 1
    elseif y == m*n + 3
      #          y                - 2 + 1 ... for special case of right edge 
      # (m*n + 3) + mod(x - 3, m) - 2 + 1
      return m*n + 2 + mod(x - 3, m)
    # general case within the matrix
    # ... horizontal ... indexed first
    elseif mod(x - y, m) == 0
      return y - 1
    # ... vertical 
    else
      # max index for horiz. case  +  idx of x  -  col_number (because there are 1 less vertical connections than #of rows)
      # (m*n + 3) + m - 2          +  x - 2     -  div(x-3, n)
      return m*(n + 1) + x - 1 - div(x-3, m)
    end
end

# 
function aux_matrix_connectivity_entries_to_lin_idx(x,y, m,n,s)
    if (x, y) == (1, 2)
      return 1
    elseif x == 2
      # left end of electrode -> horizontal
      # idx of y          + transpose to layer
      # y - 2 + 1         + m*(n+1)*div(y-3, m*n)
      return y - 1 + m*div(y-3, m*n)
    elseif y == m*n*s + 3
      # right end of electrode -> horizontal
      # 1 + #of standard horizontal of the first layer    + #row                + transpose to layer
      # 1 + m*n                                           + mod(x-3, m) + 1     + m*(n + 1)*div(x-3, m*n)
      return m*n + 2 + mod(x-3,m) + m*(n+1)*div(x-3, m*n)
    elseif div(x-3, m*n) == div(y-3, m*n) 
      if mod(x - y, m) == 0
        # horizontal
        # idx of y          + transpose to layer
        # y - 2 + 1         + m*(n+1)*div(x-3, m*n)
        return y - 1 + m*div(y-3, m*n)
      else
        # vertical
        # 1 + #horizontal    + idx of x    - #col   (because there is one connection less than row entries in each column)
        # 1 + m*(n+1)*s     + x - 2       - div(x-3, m)
        return m*(n+1)*s - 1 + x - div(x-3, m)
      end          
    else
      # interlayer
      # 1 + # horizontal  + #vertical     + idx of x
      # 1 + m*(n+1)*s     + (m-1)*n*s     + x - 2
      return (m*(2*n+1) - n)*s + x - 1     
    end
end



function ans_322(in_f, x, y)
  f(x, y) = in_f(x, y, 3,2,2)
  
  if (x,y) == (1,2); return f(x,y) == 1
  elseif (x,y) == (2,3); return f(x,y) == 2
  elseif (x,y) == (2,4); return f(x,y) == 3
  elseif (x,y) == (2,5); return f(x,y) == 4
  #
  elseif (x,y) == (3,6); return f(x,y) == 5
  elseif (x,y) == (4,7); return f(x,y) == 6
  elseif (x,y) == (5,8); return f(x,y) == 7
  #
  elseif (x,y) == (6,15); return f(x,y) == 8
  elseif (x,y) == (7,15); return f(x,y) == 9
  elseif (x,y) == (8,15); return f(x,y) == 10
  #
  #
  elseif (x,y) == (2,9); return f(x,y) == 11
  elseif (x,y) == (2,10); return f(x,y) == 12
  elseif (x,y) == (2,11); return f(x,y) == 13
  #
  elseif (x,y) == (9,12); return f(x,y) == 14
  elseif (x,y) == (10,13); return f(x,y) == 15
  elseif (x,y) == (11,14); return f(x,y) == 16
  #
  elseif (x,y) == (12,15); return f(x,y) == 17
  elseif (x,y) == (13,15); return f(x,y) == 18
  elseif (x,y) == (14,15); return f(x,y) == 19
  else
    return "ach jo"
  end
end

max_lin_idx(m,n) = aux_matrix_connectivity_entries_to_lin_idx(m*n + 1, m*n + 2, m, n)
max_lin_idx(m,n,s) = aux_matrix_connectivity_entries_to_lin_idx(m*n*(s-1) + 2, m*n*s + 2, m, n, s)








function add_value_to_matrix_row(matrix_header, equation_row, previous_id, next_id, value, m, n)
    lin_idx = aux_matrix_connectivity_entries_to_lin_idx(previous_id, next_id, m, n)     
    push!(equation_row, (lin_idx, value))
end














get_previous_list() = [(0, -1), (-1, 0)]
get_next_list() = [(1, 0), (0, 1)]

function previous(i, j, auxilary_A)    
  if j == 0
    return [1]
  else
    res = Array{Int64}(undef, length(get_previous_list()))
    for (k, d) in enumerate(get_previous_list())        
      res[k] = auxilary_A[(i+1, j+1) .+ d...]                      
    end      
    deleteat!(res, findall(x->x==-1,res))      
    return res
  end
end

function next(i,j, auxilary_A)
  if j == 0      
    return auxilary_A[2:end-1, 2]
  else
    res = Array{Int64}(undef, length(get_next_list()))
    for (k, d) in enumerate(get_next_list())
      res[k] = auxilary_A[(i+1, j+1) .+ d...]        
    end      
    deleteat!(res, findall(x->x==-1,res))
    return res
  end
end

function add_row_to_sparse_input!(sparse_input,
                                  new_row, 
                                  row_idx)        
  #@show row_idx
  for item in new_row
    push!(sparse_input[1], row_idx)
    push!(sparse_input[2], item[1])
    push!(sparse_input[3], typeof(item[2]) <: Number ? float(item[2]) : item[2])
  end
end

function add_equation_I_row(i, j, auxilary_A, m, n, matrix_header, sys_row_idx, sparse_input, RHS)
  new_row = []
  act_id = auxilary_A[i+1, j+1]               
  #
  for node_id in previous(i,j, auxilary_A)                
      add_value_to_matrix_row(matrix_header, new_row, node_id, act_id, -1, m, n)
  end
  for node_id in next(i,j, auxilary_A)        
      add_value_to_matrix_row(matrix_header, new_row, act_id, node_id, 1, m, n)
  end
  
  sort!(new_row, by = first)
  
  sys_row_idx[1] += 1
  add_row_to_sparse_input!(sparse_input, new_row, sys_row_idx[1])    
  push!(RHS, 0.0) 
  return
end






















function add_U_eq_for_each_relevant_path(Z_vector, auxilary_A, m, n, matrix_header, sys_row_idx, sparse_input, RHS)
  A = auxilary_A[2:end-1, 2:end-1]
  
  for row in 1:m
    for fall_point in n+1:-1: (row != m ? 1 : n+1)  
      new_row = []
      # first index of each path is 2
      idx1 = 2
      column = 1
      act_row = row
      while column < n + 1        
        idx2 = A[act_row, column]
        add_value_to_matrix_row(
            matrix_header, 
            new_row, 
            idx1, idx2,
            Z_vector[aux_matrix_connectivity_entries_to_lin_idx(idx1, idx2, m, n)],
            m, n
        )
        idx1 = idx2
        #
        if fall_point == column && act_row == row
          act_row += 1
        else
          column += 1
        end          
      end
      # last index of each path is m*n + 3
      idx2 = m*n + 3
      add_value_to_matrix_row(
          matrix_header, 
          new_row, 
          idx1, idx2,
          Z_vector[aux_matrix_connectivity_entries_to_lin_idx(idx1, idx2, m, n)],
          m, n
      )
      #
      sort!(new_row, by = first)
      
      sys_row_idx[1] += 1
      add_row_to_sparse_input!(sparse_input, new_row, sys_row_idx[1])   
      push!(RHS, 1.0)      
    end
  end
  
  return
end


# auxilary_A has 2 on left, m*n+3 on the right, -1 on top and bottom
function vector_to_lin_sys(Z_vector, auxilary_A)  
  # sparse_input =  (     sparse_input_row_idxs, 
  #                       sparse_input_col_idxs,
  #                       sparse_input_vals )
  matrix_header = String[]
  
  sparse_input = (Int64[], Int64[], Union{Float64, Function}[])
  sys_row_idx = [0]
  
  RHS = Float64[]

  m,n = size(auxilary_A) .-(2,2)
  
  # I verticesnext(i,j)
  # initial vertex "2"
  add_equation_I_row(1,0, auxilary_A, m, n, matrix_header, sys_row_idx, sparse_input, RHS)
  
  # other vertices
  for i in 1:m, j in 1:n
      add_equation_I_row(i,j, auxilary_A, m, n, matrix_header, sys_row_idx, sparse_input, RHS)
  end
    
  # U equations
  add_U_eq_for_each_relevant_path(Z_vector, auxilary_A, m, n, matrix_header, sys_row_idx, sparse_input, RHS)
  
  return matrix_header, sparse_input, RHS  
end






































function get_auxilary_graph_matrix(m, n)
    A = Matrix{Int64}(undef, m+2, n+2)
    #
    aux_vec = Vector{Int64}(undef, m*n)
    for i in 1:length(aux_vec)
        aux_vec[i] = i+2
    end
    pre_A = reshape(aux_vec, m, n)
    #
    A[2:end-1, 2:end-1] = pre_A
    A[:, 1] .= 2
    A[:, end] .= m*n+3
    A[1, :] .= -1
    A[end, :] .= -1
    return A
end


function get_large_material_matrix(pre_A)
    (m,n) = size(pre_A)
    A = Matrix{Int64}(undef, m+2, n+2)
    A[2:end-1, 2:end-1] = pre_A
    A[:, 1] .= i_LSM
    A[:, end] .= i_LSM
    A[1, :] .= -1
    A[end, :] .= -1    
    return A
end

function add_Z_to_vector!(Z_vector, i,j, k, l, aux_A, large_material_A, p, m, n)
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
    
    current_key = aux_matrix_connectivity_entries_to_lin_idx(idx1, idx2, m, n)
    if Z_vector[current_key] == 0
        Z_vector[current_key] = 
            get_Z_entry_from_material_matrix_codes(n1, n2, p)            
    else            
        f = deepcopy(Z_vector[current_key])            
        Z_vector[current_key] = w -> f(w) +
            get_Z_entry_from_material_matrix_codes(n1, n2, p)(w)
    end
end


function material_matrix_to_lin_sys(
            material_A::Matrix{<:Integer}=[0 0 0; 1 1 1],
            p::parameters=parameters()
    )
    (m,n) = size(material_A)
    aux_A = get_auxilary_graph_matrix(m,n)
    large_material_A = get_large_material_matrix(material_A)
    
    Z_vector = Vector{Union{Function}}(undef, max_lin_idx(m,n))
    Z_vector .= w -> 0
    for i in 2:m+1
      for j in 2:n+1
        for d in vcat(get_previous_list(), get_next_list())          
          add_Z_to_vector!(Z_vector, i,j, d..., aux_A, large_material_A, p, m, n)
        end
      end
    end
    for i in 2:m+1
      add_Z_to_vector!(Z_vector, i,1, 0, +1, aux_A, large_material_A, p, m, n)
      add_Z_to_vector!(Z_vector, i,n+2, 0, -1, aux_A, large_material_A, p, m, n)
    end    
    #return Z_vector
    return vector_to_lin_sys(Z_vector, aux_A)    
end
