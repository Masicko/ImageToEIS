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
                      parameters,
                      #
                      return_R_RC=true, 
                      TPE_warning=false,                      
                      pyplot=false,
                      
        )
      )
    end
    
    push!(R_list, [a[1] for a in res])
    push!(R_pol_list, [a[2] for a in res])
    push!(C_pol_list, [a[3] for a in res])     
  end
  return LSM_ratio_list, R_list, R_pol_list, C_pol_list
end


function run_par_study()
  #todo
end




function plot_par_study_results(x, R, Rp, Cp, label="")
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
  
  subplot(222)  
  plot(x, [sum(Cp[n]/length(Cp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel("LSM ratio")
  ylabel("C_pol")
  legend()
  
  subplot(224)  
  plot(x, [sum(R[n])/length(R[n]) for n in 1:length(R)] .+ [sum(Rp[n]/length(Rp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel("LSM ratio")
  ylabel("R_tot")
  legend()
end


function evaluate_slurm_results(dir="src/", changing_prm_name="hole_ratio"; plot=false, plot_range = nothing )
  
  function get_the_strs_from_file(path)
    output_line = ""
    changing_prm_direct_string = ""
    open(path) do f
      for line in eachline(f)         
        if length(line) >= 3 && line[1:3] == "Ima"              
          output_line = deepcopy(line)         
        elseif (length(line) >= length(changing_prm_name) &&
                line[1:length(changing_prm_name)] == changing_prm_name)
          changing_prm_direct_string = deepcopy(line)
        elseif (length(line) >= length("INPUT") &&
                line[1:length("INPUT")] == "INPUT")
          eval(Meta.parse(line))          
        end
      end
    end
    return changing_prm_direct_string, output_line
  end
  
  output_dict = Dict()
    
  for file in readdir(dir)
    if length(file) >= 6 && file[1:6] == "slurm-"
      changing_prm_direct_string, my_str = get_the_strs_from_file(dir*file)
      if my_str == ""
        println("file $(file) skipped")
      else        
        last_is_equal = findall(x -> x == '=', my_str)[end]        
        
        if changing_prm_direct_string != ""          
          prm_value = eval(
              Meta.parse(
                split(changing_prm_direct_string, '=')[end]
              )
          )                        
        else
          head = my_str[1:last_is_equal-1]
          first_bracket = findall(x -> x == '(', head)[1]
          #
          NT_prms = eval(Meta.parse(head[first_bracket : end]))
          physical_prms = NT_prms[:parameters]
          for pair in physical_prms
            if pair[1] == changing_prm_name
              prm_value = pair[2]
            end
          end                                        
        end
        
        ##############################
        # has fields (x, R, Rp, Cp)
        # -> x = LSM_ratio_list
        # -> R = [R(x) for x in LSM_ratio_list]        
        output_tuple = eval(Meta.parse(my_str[last_is_equal + 1 : end]))
        prm_value_identifier = prm_value
        
        #@show output_dict        
        if haskey(output_dict, prm_value_identifier)
          recent_tuple = output_dict[prm_value_identifier]
          if recent_tuple[1] != convert.(Float32, output_tuple[1])
            println("ERROR: LSM_ratio_list mismatch: $(recent_tuple[1]) != $(output_tuple[1])")
          end
          for (i, ratio) in enumerate(output_tuple[1])
            for prm_identifier in 2:4                            
              append!(recent_tuple[prm_identifier][i], output_tuple[prm_identifier][i])
            end
          end
        else                          
          output_dict[prm_value_identifier] = map(x ->  map(
                                                          y -> convert.(Float32, y)
                                                          ,
                                                          x
                                                        )                                                        
                                                  , 
                                                  output_tuple
                                              )
        end
      end
    end
  end

  
  if plot 
    for plotted_prm in sort!([deepcopy(keys(output_dict))...])
      key = plotted_prm 
      ImageToEIS.plot_par_study_results(output_dict[key]..., "$(changing_prm_name) = $(key)")
    end
  end
  
  return output_dict
end








