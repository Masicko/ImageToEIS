
function generate_matrix(dimensions, porosity::Float64, LSM_ratio::Float64, pore_prob::Union{Float64, Nothing}=nothing; recursion_depth=1e12, check_connectivity=true)
  if typeof(pore_prob) == Nothing
    the_domain = Array{Int16}(undef, dimensions)
    repetition_number = 1
    while repetition_number < recursion_depth
      for i in 1:length(the_domain[:])
        if rand() <= porosity
          the_domain[i] = i_hole
        else
          if rand() <= LSM_ratio
            the_domain[i] = i_LSM
          else
            the_domain[i] = i_YSZ
          end
        end
      end      
      if !check_connectivity || check_material_connection(the_domain)
        return the_domain
      else
        repetition_number += 1
      end
    end
    
    println("ERROR: recursion_depth = 0 ... no trials left ... material_connectivity is not ensured")
    return -1
  else
    return shoot_pores(dimensions, porosity, LSM_ratio, pore_prob)
  end
end


function generate_submatrix_to_matrix(matrix, left_upper::Union{Tuple, Array}, right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
  submatrix = generate_matrix(right_lower .- left_upper .+ (1,1), porosity, LSM_ratio)
  matrix[left_upper[1] : right_lower[1], left_upper[2] : right_lower[2]] = deepcopy(submatrix)
  return 
end

function generate_submatrix_to_matrix(matrix, right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
  generate_submatrix_to_matrix(matrix, (1,1), right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
end

function get_default_submatrix_list()
  # 1 item in the resulting list is
  # 
  # left upper, right lower corner, porosity, LSM_ratio

  first_block_column = 1
  second_block_column = 21
  third_block_column = 41
  width = 42
  height = 80
  return [
    [(1,1), (height, first_block_column), 1.0, 0.5],
    [(1,first_block_column + 1), (height, second_block_column), 0.2, 0.5],
    [(1,second_block_column + 1), (height, third_block_column), 0.2, 1.0],
    [(1,third_block_column + 1), (height, width), 0.2, 1.0],
    #
    [(10, 1), (15, first_block_column), 0.0, 1.0],
    [(40, 1), (45, first_block_column), 0.0, 1.0],
    [(60, 1), (65, first_block_column), 0.0, 1.0]
  ]
end

function generate_matrix(submatrix_list::Array=get_default_submatrix_list())
  max_height, max_width = -1, -1
  for submatrix in submatrix_list
    s = submatrix[2]
    if s[1] > max_height
      max_height = s[1]
    end
    if s[2] > max_width
      max_width = s[2]
    end    
  end  
  
  output_matrix = Matrix{Int16}(undef, max_height, max_width)
  output_matrix .= -1

  for submatrix in submatrix_list
    generate_submatrix_to_matrix(output_matrix, submatrix...)    
  end 
  invalid_submatrix_list = false
  foreach(x -> x==-1 ?  invalid_submatrix_list = true : false, output_matrix)
  if invalid_submatrix_list
    println("ERROR: invalid_submatrix_list input")
    return throw(Exception)
  else
    return output_matrix
  end
end



function three_column_domain_template(LSM_ratio1, LSM_ratio2, LSM_ratio3; 
                              #
                              porosity1=0.5, porosity2=0.5, porosity3=0.5,
                              #                              
                              positions_of_contacts=[15, 50], height_of_contacts=5, 
                              #
                              column_width = 5,
                              #
                              height = 70
                              )
  # 1 item in the resulting list is
  # 
  # left upper, right lower corner, porosity, LSM_ratio
  first_block_column = 1
  second_block_column = column_width + 1
  third_block_column = 2*column_width + 1
  fourth_block_column = 3*column_width + 1
  width = fourth_block_column + 1
  
  output = [
    [(1,1), (height, first_block_column), 1.0, 0.5],
    [(1,first_block_column + 1), (height, second_block_column), porosity1, LSM_ratio1],
    [(1,second_block_column + 1), (height, third_block_column), porosity2, LSM_ratio2],
    [(1,third_block_column + 1), (height, fourth_block_column), porosity3, LSM_ratio3],    
    [(1, fourth_block_column + 1), (height, width), 0.0, 1.0]
  ]
  
  # contacts
  for p_of_contact in positions_of_contacts
    push!(output, 
        [(p_of_contact, 1), (p_of_contact + height_of_contacts - 1, first_block_column), 0.0, 1.0]
    )
  end
  
  return output
end

function chess_matrix(m, n)
  A = Matrix{Int64}(undef, m,n)
  A .= 0
  for i in 1:m, j in 1:n
      if mod(i+j, 2) == 0
        A[i,j] = 1
      end
  end
  return A
end

function aux_domain(domain::Array{T} where T<: Integer; inner_number=nothing, boundary_number=0)
  if length(size(domain)) == 2
    aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
    if typeof(inner_number) == Nothing
      aux_matrix[2:end-1, 2:end-1] = domain
    else
      aux_matrix[:, :] .= inner_number
    end
    aux_matrix[:, 1] .= boundary_number
    aux_matrix[:, end] .= boundary_number
    #
    aux_matrix[1, :] .= boundary_number
    aux_matrix[end, :] .= boundary_number
    #
    return aux_matrix
  elseif length(size(domain)) == 3
    aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
    if typeof(inner_number) == Nothing
      aux_matrix[2:end-1, 2:end-1, 2:end-1] = domain
    else
      aux_matrix[:, :, :] .= inner_number
    end

    aux_matrix[:, end, :] .= boundary_number
    aux_matrix[:, 1, :]   .= boundary_number
    #
    aux_matrix[1, :, :]   .= boundary_number
    aux_matrix[end, :, :] .= boundary_number
    #
    aux_matrix[:, :, 1]   .= boundary_number
    aux_matrix[:, :, end] .= boundary_number
    #
    return aux_matrix
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end  
end

function search_dirs(domain::Array{T} where T<: Integer)
  if length(size(domain)) == 2
    return search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
  elseif length(size(domain)) == 3
    return search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end
end




function check_material_connection(domain::Array{T} where T <: Integer; return_alone_percentage=false)
  material_sign = 8
  end_sign = 6
  
  if length(size(domain)) == 2
    dims = 2
  else
    dims = 3
  end  
  
  if dims == 2
    aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
    aux_matrix .= 1 
    aux_matrix[:, end] .= end_sign
    aux_matrix[:, 1] .= 0
    #
    aux_matrix[1, :] .= 0
    aux_matrix[end, :] .= 0
    
    search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
  end
  
  if dims == 3
    aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
  
    aux_matrix .= 1 
    aux_matrix[:, end, :] .= end_sign
    
    aux_matrix[:, 1, :] .= 0
    #
    aux_matrix[1, :, :] .= 0
    aux_matrix[end, :, :] .= 0
    #
    aux_matrix[:, :, 1] .= 0
    aux_matrix[:, :, end] .= 0
    
    search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
  end
  
  is_connected=false  
  
  
  if dims==2    
    for row in 1:size(domain)[1]
      list_to_process = [(row +1, 2)]
      
      while length(list_to_process) > 0
        x, y = list_to_process[end]
        deleteat!(list_to_process, length(list_to_process))
        
        if aux_matrix[x,y] == end_sign        
          is_connected=true
          return true
        else
          if (aux_matrix[x,y] == 1)
            #@show x,y
            if (domain[x-1, y-1] in i_material_list)          
              aux_matrix[x,y] = material_sign         
              for dir in search_dirs
                push!(list_to_process, (x,y) .+ dir)
                #search_around_this(((x,y) .+ dir)...)
              end
            else
              aux_matrix[x,y] = 0
            end
          end
        end
      end            
    end
  
  
  else

    
    for row in 1:size(domain)[1], layer in 1:size(domain)[3]
      list_to_process = [(row +1, 2, layer+1)]
      
      while length(list_to_process) > 0
        x, y, z = list_to_process[end]
        deleteat!(list_to_process, length(list_to_process))
        
        if aux_matrix[x,y,z] == end_sign        
          is_connected=true
          return true
        else
          if (aux_matrix[x,y,z] == 1)           
            if (domain[x-1, y-1, z-1] in i_material_list)                      
              aux_matrix[x,y,z] = material_sign         
              for dir in search_dirs        
                push!(list_to_process, (x,y,z) .+ dir)                
              end
            else
              aux_matrix[x,y, z] = 0
            end        
          end
        end
      end
      
    end    
  end
  
  
  if return_alone_percentage == true
    if dims == 2
      sub_aux_matrix = aux_matrix[2:end-1, 2:end-1]
    elseif dims == 3
      sub_aux_matrix = aux_matrix[2:end-1, 2:end-1, 2:end-1]
    end
    #
    material_pixels = 0
    alone_pixels = 0
    for i in 1:length(domain[:])      
      if domain[i] in i_material_list
        material_pixels += 1
        if sub_aux_matrix[i] != material_sign
          alone_pixels += 1    
        end
      end
    end    
    return alone_pixels/material_pixels
  else
    return is_connected
  end
end



































function typical_use_shoot_pores()
  mm = ImageToEIS.shoot_pores((20, 20, 20), 0.7, 0.4, 0.0)
  ImageToEIS.matrix_to_file("jojojjoooojoooo.png", mm)
  println(ImageToEIS.check_material_connection( 
      mm
      )
  )

end




function test_shoot_pores()
  return shoot_pores((10, 10), 0.5, 0.4, 0.1)
end


function ext_correction(domain)
  if length(size(domain)) == 2
    return (1, 1)
  else
    return (1,1,1)
  end  
end

function check_pore_is_boundary(ext_domain, pos)
  is_boundary = false
  ext_corr = ext_correction(ext_domain)
  if ext_domain[pos .+ ext_corr...] == -1
    return false
  end
  for dir in search_dirs(ext_domain)
    if (ext_domain[pos .+ dir .+ ext_corr...] != i_hole) && (ext_domain[pos .+ dir .+ ext_corr...] != -1)
      is_boundary = true
    end
  end
  return is_boundary
end

function get_indicies_of_random_neighbour(ext_domain, pos)
  neigh_list = []
  for dir in search_dirs(ext_domain)
    test_pos = pos .+ dir .+ ext_correction(ext_domain)
    if (ext_domain[test_pos...] != i_hole) && (ext_domain[test_pos...] != -1)
      push!(neigh_list, test_pos)
    end
  end
  return neigh_list[Int32(rand(1:length(neigh_list)))] .- ext_correction(ext_domain)  
end

function get_standard_domain(ext_domain)
  if length(size(ext_domain)) == 2
    return ext_domain[2:end-1, 2:end-1]
  else
    return ext_domain[2:end-1, 2:end-1, 2:end-1]
  end
end

function get_body_list(domain)
  dims = size(domain)
  if length(dims) == 2
    return [(x,y) for x in 1:dims[1], y in 1:dims[2]][:]
  elseif length(size(domain)) == 3
    return [(x,y,z) for x in 1:dims[1], y in 1:dims[2], z in 1:dims[3]][:]
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end  
end

function shoot_pores(dims, porosity, LSM_ratio, pore_prob; recursion_depth=10)
  domain = Array{Integer}(undef, dims)
  boundary_pore_list = []
  body_list = get_body_list(domain)
  #
  pix_tot = prod(dims)
  pix_por = Int32(round(porosity*pix_tot))
  
  
  extended_domain = aux_domain(domain, inner_number=i_YSZ, boundary_number=-1)
  
  
  for i in 1:pix_por        
    mother_item_idx = nothing
    if (rand() <= pore_prob) && (length(boundary_pore_list) > 0)        
        mother_item_idx = rand(1:length(boundary_pore_list))        
        swapping_item_indices = get_indicies_of_random_neighbour(extended_domain, boundary_pore_list[mother_item_idx])        
        swapping_item_idx = findall(x -> x == swapping_item_indices, body_list)[1]
    else        
        swapping_item_idx = rand(1:length(body_list)) 
        swapping_item_indices = body_list[swapping_item_idx]
    end
    
    
    extended_domain[swapping_item_indices .+ ext_correction(domain)...] = i_hole
    
    if (typeof(mother_item_idx) != Nothing) && !check_pore_is_boundary(extended_domain, boundary_pore_list[mother_item_idx])
      deleteat!(boundary_pore_list, mother_item_idx)      
    end
    
    for dir in search_dirs(domain)
      
        if !check_pore_is_boundary(extended_domain, swapping_item_indices .+ dir)          
          search_result = findall(x -> x == swapping_item_indices .+ dir, boundary_pore_list)
      
          if length(search_result) > 0
            deleteat!(boundary_pore_list, search_result[1])
          end
        end
    end
    
    if check_pore_is_boundary(extended_domain, swapping_item_indices)
      if !(swapping_item_indices in boundary_pore_list)
        push!(boundary_pore_list, swapping_item_indices)
      end
    end    
    deleteat!(body_list, swapping_item_idx)
  end
  
  for i in 1:length(extended_domain[:])
    if (extended_domain[i] != -1) && (extended_domain[i] != i_hole)
      if rand() < LSM_ratio
        extended_domain[i] = i_LSM
      else
        extended_domain[i] = i_YSZ
      end
    end
  end
  
  the_domain = get_standard_domain(extended_domain)
  if check_material_connection(the_domain)
    return the_domain
  else
    if recursion_depth > 0
      return shoot_pores(dims, porosity, LSM_ratio, pore_prob, recursion_depth=recursion_depth-1)
    else
      println("ERROR: recursion_depth = 0 ... no trials left ... material_connectivity is not ensured")
      return -1
    end
  end
end



















# function make_domain_connected(domain :: Array{T} where T <: Integer)
#   if length(size(domain)) == 2
#     dims = 2
#   else
#     dims = 3
#   end  
#   
#   if dims == 2
#     aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
#     aux_matrix = -1
#     aux_matrix[2:end-1, 2:end-1] = domain
#     
#     search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
#     
#     list_to_process = [(x,y) for x in 1:size(domain[1]), y in 1:size(domain[2])]
#     while length(list_to_process) > 0
#       act_idx = rand(1:length(list_to_process))
#       act_indices = list_to_process[act_idx]
#       
#       if !check_point_connection(aux_matrix, aux_matrix[act_indices .+ 1])
#         pore_to_swap = find_a_place_for_the_point(aux_matrix)
#         
#       end
#     end
#   end
#   
#   if dims == 3
#     aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
#     aux_matrix .= -1
#     aux_matrix[2:end-1, 2:end-1, 2:end-1] = domain
#     
#     search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
#   end
#   
# end
