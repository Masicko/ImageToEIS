function export_EIS_to_Z_file(path, f_range, Z_range)
  if length(f_range) != length(Z_range)
    println("ERROR: length(f_range) != length(Z_range) >> $(length(f_range)) != $(length(Z_range))")
    return throw(Exception)
  end  
  src_dir = @__DIR__  
  open(path, "w") do f
    open(src_dir*"/Z_view_template.txt", "r") do template            
      for line in eachline(template)
        write(f, line*"\n")
      end
    end
    for i in 1:length(f_range)
      write(f, @sprintf("%E\t%E\t%E\t%E\t%E\t%E\t%E\t%i\t%i\t", 
                f_range[i], 0, 0, 0, real(Z_range[i]), imag(Z_range[i]), 0, 0, 0
                )*"\n"
      )
    end
  end
end
