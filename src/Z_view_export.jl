function get_frequency_descendent!(f_range, Z_range)
  if f_range[1] <= f_range[end]
    des_f_range = reverse(f_range)
    des_Z_range = reverse(Z_range)
  else
    des_f_range = f_range
    des_Z_range = Z_range
  end
  return des_f_range, des_Z_range
end

function export_EIS_to_Z_file(path, f_range, Z_range)
  if length(f_range) != length(Z_range)
    println("ERROR: length(f_range) != length(Z_range) >> $(length(f_range)) != $(length(Z_range))")
    return throw(Exception)
  end  
  
  des_f_range, des_Z_range = get_frequency_descendent!(f_range, Z_range)
  
  src_dir = @__DIR__  
  mkpath(dirname(path))
  
  open(path, "w") do f
    open(src_dir*"/Z_view_template.txt", "r") do template            
      for line in eachline(template)
        write(f, line*"\n")
      end
    end
    
    for i in 1:length(des_f_range)
      write(f, @sprintf("%E\t%E\t%E\t%E\t%E\t%E\t%E\t%i\t%i\t", 
                des_f_range[i], 0, 0, 0, real(des_Z_range[i]), imag(des_Z_range[i]), 0, 0, 0
                )*"\n"
      )
    end
  end
end
