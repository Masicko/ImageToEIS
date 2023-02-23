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

max_lin_idx(m,n) =  if m==1
                      aux_matrix_connectivity_entries_to_lin_idx(m*n + 2, m*n + 3, m, n)
                    else
                      aux_matrix_connectivity_entries_to_lin_idx(m*n + 1, m*n + 2, m, n)                      
                    end
max_lin_idx(m,n,s) =  if s==1
                        max_lin_idx(m,n)
                      else
                        aux_matrix_connectivity_entries_to_lin_idx(m*n*(s-1) + 2, m*n*s + 2, m, n, s)
                      end








function add_value_to_matrix_row(matrix_header, equation_row, unknown_id, value)
    push!(matrix_header, (unknown_id, (unknown_id)))
    
    push!(equation_row, (unknown_id, value))
end














get_previous_list(m, n) = [(0, -1), (-1, 0)]
get_next_list(m, n) = [(1, 0), (0, 1)]

get_previous_list(m, n, s) = [(0, -1, 0), (-1, 0, 0), (0, 0, -1)]
get_next_list(m, n, s) = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]



function previous(i, j, k, auxilary_A)
  if j == 0
    return [1]
  else
    res = Array{Int64}(undef, length(get_previous_list(0,0,0)))
    for (l, d) in enumerate(get_previous_list(0,0,0))        
      res[l] = auxilary_A[(i+1, j+1, k+1) .+ d...]                      
    end      
    deleteat!(res, findall(x->x==-1,res))      
    return res
  end
end

function previous(i, j, auxilary_A)    
  if j == 0
    return [1]
  else
    res = Array{Int64}(undef, length(get_previous_list(0,0)))
    for (k, d) in enumerate(get_previous_list(0,0))        
      res[k] = auxilary_A[(i+1, j+1) .+ d...]                      
    end      
    deleteat!(res, findall(x->x==-1,res))      
    return res
  end
end

function next(i,j,k, auxilary_A)
  if j == 0          
    return auxilary_A[2:end-1, 2, 2:end-1][:]
  else
    res = Array{Int64}(undef, length(get_next_list(0,0,0)))        
    for (l, d) in enumerate(get_next_list(0,0,0))      
      res[l] = auxilary_A[(i+1, j+1, k+1) .+ d...]        
    end          
    deleteat!(res, findall(x->x==-1,res))
    return res
  end
end

function next(i,j, auxilary_A)
  if j == 0      
    return auxilary_A[2:end-1, 2]
  else
    res = Array{Int64}(undef, length(get_next_list(0,0)))
    for (k, d) in enumerate(get_next_list(0,0))
      res[k] = auxilary_A[(i+1, j+1) .+ d...]        
    end      
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



function add_equation_I_row(pos, Z_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
  new_row = []
  act_id = auxilary_A[pos .+ 1 ...]               
  #
  Z_sum = 0 + 0*im
  for node_id in previous(pos..., auxilary_A)
    minus_one_over_Z_ji = w -> -1/(Z_vector[aux_matrix_connectivity_entries_to_lin_idx(node_id, act_id, dims...)](w))
    f(w) = deepcopy(Z_sum)(w)
    Z_sum = w -> f(w) + minus_one_over_Z_ji(w)
    add_value_to_matrix_row(matrix_header, new_row, 
                            node_id, 
                            minus_one_over_Z_ji)
  end
  for node_id in next(pos..., auxilary_A)    
    one_over_Z_ij = w -> 1/(Z_vector[aux_matrix_connectivity_entries_to_lin_idx(act_id, node_id, dims...)](w))
    f(w) = deepcopy(Z_sum)(w)
    Z_sum = w -> f(w) + one_over_Z_ij(w)
    add_value_to_matrix_row(matrix_header, new_row, 
                            node_id, 
                            one_over_Z_ij)
  end
  add_value_to_matrix_row(matrix_header, new_row, 
    act_id, 
    Z_sum
  )  
  
  sort!(new_row, by = first) 
  
  sys_row_idx[1] += 1
  add_row_to_sparse_input!(sparse_input, new_row, sys_row_idx[1])    
  push!(RHS, 0.0) 
  return
end







function get_readable_matrix_header(matrix_header)
  sort!(matrix_header, by = first)
  res = [matrix_header[1]]  
  for item in matrix_header
    if item[1] != res[end][1]
      push!(res, item)
    end
  end
  return [item[2] for item in res]
end


# auxilary_A has 2 on left, m*n+3 on the right, -1 on top and bottom
function vector_to_lin_sys(Z_vector, auxilary_A)  
  # sparse_input =  (     sparse_input_row_idxs, 
  #                       sparse_input_col_idxs,
  #                       sparse_input_vals )
  matrix_header = []
  
  sparse_input = (Int64[], Int64[], Union{Float64, Function}[])
  sys_row_idx = [0]
  
  RHS = Float64[]

  dims = size(auxilary_A) .- 2
  
  mode_3D = length(dims) == 3

  # other vertices
  if mode_3D
    for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]
      add_equation_I_row((i,j, k), Z_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
    end
  else
    for i in 1:dims[1], j in 1:dims[2]
        add_equation_I_row((i,j), Z_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
    end
  end
  
  matrix_header_output = get_readable_matrix_header(matrix_header)
  
  return matrix_header_output, sparse_input, RHS
end



































function get_auxilary_graph_matrix(m, n)
    A_aux = Array{Int64}(undef, m+2, n+2)
    #
    aux_vec = Vector{Int64}(undef, m*n)
    for i in 1:length(aux_vec)
        aux_vec[i] = i+2
    end
    pre_A = reshape(aux_vec, m, n)
    #
    A_aux[2:end-1, 2:end-1] = pre_A
    A_aux[:, 1] .= 2
    A_aux[:, end] .= m*n+3
    A_aux[1, :] .= -1
    A_aux[end, :] .= -1    
    return A_aux
end


function get_auxilary_graph_matrix(m, n, s)
    A_aux = Array{Int64}(undef, m+2, n+2, s+2)
    #
    aux_vec = Vector{Int64}(undef, m*n*s)
    for i in 1:length(aux_vec)
        aux_vec[i] = i+2
    end
    pre_A = reshape(aux_vec, m, n, s)
    #
    A_aux[2:end-1, 2:end-1, 2:end-1] = pre_A
    A_aux[:, 1, :] .= 2
    A_aux[:, end, :] .= m*n*s+3
    A_aux[1, :,  :] .= -1
    A_aux[end, :, :] .= -1
    A_aux[:, :, 1] .= -1
    A_aux[:, :, end] .= -1
    return A_aux
end

function get_large_material_matrix(pre_A::Array{<:Integer, 2})
  (m,n) = size(pre_A)
  A = Array{Int64}(undef, m+2, n+2)
  A[2:end-1, 2:end-1] = pre_A
  A[:, 1] .= i_LSM
  A[:, end] .= i_LSM
  A[1, :] .= -1
  A[end, :] .= -1              
  return A
end

function get_large_material_matrix(pre_A::Array{<:Integer, 3})            
    (m,n,s) = size(pre_A)
    A = Array{Int64}(undef, m+2, n+2, s+2)
    A[2:end-1, 2:end-1, 2:end-1] = pre_A
    A[:, 1, :] .= i_LSM
    A[:, end, :] .= i_LSM
    A[1, :, :] .= -1
    A[end, :, :] .= -1        
    A[:, :, 1] .= -1
    A[:, :, end] .= -1
    return A
end

function add_Z_to_vector!(Z_vector, pos, dir, aux_A, large_material_A, p, dims)
    n1 = large_material_A[pos...]
    n2 = large_material_A[pos .+ dir...]
    if n2 == -1
        return  
    end
    
    idx1 = aux_A[pos...]
    idx2 = aux_A[pos .+ dir...]
    if idx2 < idx1
        aux = idx1
        idx1 = idx2
        idx2 = aux
    end
    
    current_key = aux_matrix_connectivity_entries_to_lin_idx(idx1, idx2, dims...)
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
            material_A::Array{<:Integer}=[0 0 0; 1 1 1],
            p::parameters=parameters()
    )
    
    dims = size(material_A)
    mode_3D = length(dims) == 3
    aux_A = get_auxilary_graph_matrix(dims...)
    large_material_A = get_large_material_matrix(material_A)
    
    Z_vector = Vector{Union{Function}}(undef, max_lin_idx(dims...))
    Z_vector .= w -> 0
    for d in vcat(get_previous_list(dims...), get_next_list(dims...))  
      for i in 2:dims[1]+1
        for j in 2:dims[2]+1                
          if mode_3D
            for k in 2:dims[3]+1
              add_Z_to_vector!(Z_vector, (i, j, k), d, aux_A, large_material_A, p, dims)
            end
          else
            add_Z_to_vector!(Z_vector, (i, j), d, aux_A, large_material_A, p, dims)
          end
        end
      end
    end
    for i in 2:dims[1]+1
      if mode_3D
        for k in 2:dims[3]+1
          add_Z_to_vector!(Z_vector, (i,1, k), (0, +1, 0), aux_A, large_material_A, p, dims)
          add_Z_to_vector!(Z_vector, (i,dims[2]+2, k), (0, -1, 0), aux_A, large_material_A, p, dims)
        end
      else
        add_Z_to_vector!(Z_vector, (i,1), (0, +1), aux_A, large_material_A, p, dims)
        add_Z_to_vector!(Z_vector, (i,dims[2]+2), (0, -1), aux_A, large_material_A, p, dims)
      end
    end    
    #return Z_vector, aux_A, large_material_A
    return vector_to_lin_sys(Z_vector, aux_A)    
end
