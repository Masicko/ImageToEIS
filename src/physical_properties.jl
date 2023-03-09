const i_YSZ = 0
const i_LSM = 1
const i_hole = 2

const i_material_list = [0,1]


Base.@kwdef mutable struct parameters        
    R_YSZ::Float64 = 1/0.045 # S/cm
    R_pol_YSZ::Float64 = 0
    C_pol_YSZ::Float64 = 0.001
    #
    R_LSM::Float64 = 1/290 # S/cm
    R_pol_LSM::Float64 = 40
    C_pol_LSM::Float64 = 0.005
    #
    R_hole::Float64 = 1e12
end

function RC_el(R, C, w)
  return (R/(1 + C*R*w*im))
end

function get_Y_entry_from_material_matrix_codes(n1, n2, p::parameters)
    if      n1 == i_LSM
            if      n2 == i_LSM
              return w -> 1/(p.R_LSM/2)
            elseif  n2 == i_YSZ
              if p.R_pol_LSM == 0.0
                return w -> 1/(p.R_LSM/2)
              else
                return w -> 1/(p.R_LSM/2 + RC_el(p.R_pol_LSM, p.C_pol_LSM, w))
              end
            elseif  n2 == i_hole
              return w ->  1/(p.R_LSM/2)
            end
    elseif  n1 == i_YSZ
            if      n2 == i_LSM
              if p.R_pol_YSZ == 0.0
                return w -> 1/(p.R_YSZ/2)
              else          
                return w -> 1/(p.R_YSZ/2 + RC_el(p.R_pol_YSZ, p.C_pol_YSZ, w))
              end             
            elseif  n2 == i_YSZ
              return w -> 1/(p.R_YSZ/2)
            elseif  n2 == i_hole
              return w -> 1/(p.R_YSZ/2)
            end
    elseif  n1 == i_hole            
            return w -> 1/(p.R_hole/2)
    else
        println("ERROR: get_Y_entry...")
    end        
end



function file_to_matrix(path="src/geometry.png")
  RGB_m = load(path)
  m = Matrix{Int}(undef, size(RGB_m)...)
  for i in 1:size(m)[1]
    for j in 1:size(m)[2]
      (r,g,b) = (RGB_m[i,j].r, RGB_m[i,j].g, RGB_m[i,j].b)
      if (r > 0.5) && (g > 0.5) && (b > 0.5)
        m[i,j] = i_hole
      elseif (r > 0.5) && (g > 0.5) 
        m[i,j] = i_YSZ
      else
        m[i,j] = i_LSM
      end
    end
  end
  return m
end

function matrix_to_file(path, matrix, colors=((1,1,0), (0,0,0)))
  save(path, map(x -> if     x == i_hole
                   RGB(1, 1, 1)
                elseif x == i_YSZ
                   RGB(colors[1]...)
                elseif x == i_LSM
                   RGB(colors[2]...)
                end, matrix)
  )
  return
end
