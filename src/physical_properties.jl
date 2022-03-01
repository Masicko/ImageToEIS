const i_YSZ = 0
const i_LSM = 1
const i_dira = 2


Base.@kwdef mutable struct parameters
    R_YSZ::Float64 = 100
    R_LSM::Float64 = 1
    R_pol::Float64 = 40
    C_pol::Float64 = 0.001
    R_dira::Float64 = 1000000
end

function get_Z_entry_from_material_matrix_codes(n1, n2, p::parameters)
    
    if      n1 == i_LSM
            if      n2 == i_LSM
              return "p.R_LSM/2"
            elseif  n2 == i_YSZ
              return "p.R_LSM/2 + (p.R_pol/(1 + p.C_pol*p.R_pol*w*im))"
            elseif  n2 == i_dira
              return "p.R_LSM/2"
            end
    elseif  n1 == i_YSZ
            return "p.R_YSZ/2"
    elseif  n1 == i_dira
            return "p.R_dira/2"
    else
        println("ERROR: get_Z_entry...")
    end    
end


function file_to_matrix(path="src/geometry.png")
  RGB_m = load(path)
  m = Matrix(undef, size(RGB_m)...)
  for i in 1:size(m)[1]
    for j in 1:size(m)[2]
      (r,g,b) = (RGB_m[i,j].r, RGB_m[i,j].g, RGB_m[i,j].b)
      if (r > 0.5) && (g > 0.5) && (b > 0.5)
        m[i,j] = i_dira
      elseif (r > 0.5) && (g > 0.5) 
        m[i,j] = i_YSZ
      else
        m[i,j] = i_LSM
      end
    end
  end
  return m
end
