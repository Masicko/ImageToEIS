function plot_results(n_list, timings)
  figure(1)
  title("Speed of ImageToEIS")  
  xlabel("n (nxn image)")
  ylabel("t [s]")
  plot(Int.(n_list), timings, "-x")
end


function elapsed_test()
  timings = []
  
  n_list = collect(5 : 5 : 50)
  for n in n_list
      println(" --- n = $(n) ---")
      push!(timings, 
          @elapsed image_to_EIS_2D(zeros(n,n), f_range=[0.01, 100000], pyplot=false)
      )
  end
  plot_results(n_list, timings)
  return n_list, timings
end






function par_study(;
                    LSM_ratio_list = collect( 0 : 0.2 : 1.0),
                    hole_ratio = 0.2,
                    dimensions = (20,20),
                    repeted_trials = 10,
                    parameters = []
  )
  function generate_random_specification(LSM_ratio, hole_ratio)
    v = Array{Int16}(undef, 100)    
    
      
    hole_length = Int32(round(hole_ratio*100))
    LSM_length = Int32(round((1 - hole_ratio)*LSM_ratio*100))        
    
    v .= i_YSZ
    v[1 : hole_length] .= i_hole
    v[hole_length + 1 : hole_length + LSM_length] .= i_LSM
    return v
  end
  
  R_list = []
  R_pol_list = []
  C_pol_list = []
  
  for LSM_ratio in LSM_ratio_list
    @show LSM_ratio
    
    res = []    
    for i in 1:repeted_trials
      push!(res,
        image_to_EIS( rand(
                          generate_random_specification(LSM_ratio, hole_ratio), 
                          dimensions...
                      ),
                      return_R_RC=true, 
                      TPE_warning=false,                      
                      pyplot=false 
        )
      )
    end
    
    push!(R_list, [a[1] for a in res]/repeted_trials)
    push!(R_pol_list, [a[2] for a in res]/repeted_trials)
    push!(C_pol_list, [a[3] for a in res]/repeted_trials)     
  end
  return LSM_ratio_list, R_list, R_pol_list, C_pol_list
end


function run_par_study()
  #todo
end


function evaluate_slurm_results(dir="src/")
  
  function get_the_strs_from_file(path)
    output_line = ""
    hole_ratio_string = ""
    open(path) do f
      for line in eachline(f)                 
        if length(line) >= 3 && line[1:3] == "Ima"              
          output_line = deepcopy(line)         
        elseif length(line) >= 10 && line[1:10] == "hole_ratio"
          hole_ratio_string = deepcopy(line)
        end
      end
      end
    return hole_ratio_string, output_line
  end
  
  dict_for_hole_ratios = Dict()
  
  hole_ratio_identifier = ""
  for file in readdir(dir)
    if length(file) >= 6 && file[1:6] == "slurm-"
      hole_ratio_string, my_str = get_the_strs_from_file(dir*file)
      if hole_ratio_string == "" || my_str == ""
        println("file $(file) skipped")
      else        
        last_is_equal = findall(x -> x == '=', my_str)[end]
        
        #head = my_str[1:last_is_equal-1]
        #first_bracket = findall(x -> x == '(', head)[1]
        
        ##############################
        # has fields (x, R, Rp, Cp)
        # -> x = LSM_ratio_list
        # -> R = [R(x) for x in LSM_ratio_list]
        output_tuple = eval(Meta.parse(my_str[last_is_equal + 1 : end]))
        
        
        hole_ratio_identifier = string(
            eval(Meta.parse(
            split(hole_ratio_string, '=')[end]
          ))
        )
        
              
        #@show dict_for_hole_ratios
        if haskey(dict_for_hole_ratios, hole_ratio_identifier)
          recent_tuple = dict_for_hole_ratios[hole_ratio_identifier]
          if recent_tuple[1] != output_tuple[1]
            println("ERROR: LSM_ratio_list mismatch: $(recent_tuple[1]) != $(output_tuple[1])")
          end
          for (i, ratio) in enumerate(output_tuple[1])
            for prm_identifier in 2:4
              append!(recent_tuple[prm_identifier][i], output_tuple[prm_identifier][i])
            end
          end
        else                
          dict_for_hole_ratios[hole_ratio_identifier] = output_tuple        
        end
      end
    end
  end
  
  return dict_for_hole_ratios
end



function plot_par_study_results(x, R, Rp, Cp, label)
  figure(5)
  suptitle("YSZ-LSM-hole 50x50 random distribution test - 10 repetitions for each conditions")
  subplot(221)
  plot(x, [sum(R[n])/length(R[n]) for n in 1:length(R)], label=label, "-x")
  xlabel("LSM ratio")
  ylabel("R_ohm")
  legend()
  
  subplot(223)    
  plot(x, [sum(Rp[n]/length(Rp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel("LSM ratio")
  ylabel("R_pol")
  legend()
  
  subplot(122)  
  plot(x, [sum(Cp[n]/length(Cp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel("LSM ratio")
  ylabel("C_pol")
  legend()
end







