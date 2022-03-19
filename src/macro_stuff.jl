# time testing
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




function homogenous_matrix(par_study_prms::Dict)
  return rand(
            generate_random_specification(
              par_study_prms["LSM_ratio"], 
              par_study_prms["hole_ratio"]
            ), 
            par_study_prms["dimensions"]...
  )
end

function three_column_domain_matrix(p::Dict)
  return generate_matrix(
    three_column_domain_template(p["LSM_ratio1"], p["LSM_ratio2"], p["LSM_ratio3"], 
                              #
                              hole_ratio1=p["hole_ratio1"], hole_ratio2=p["hole_ratio2"], hole_ratio3=p["hole_ratio3"],
                              #                              
                              positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
                              #
                              column_width = p["column_width"],
                              #
                              height = p["height"]
    )
  )
end

function three_column_domain_LSM_ratios(p::Dict)
  return generate_matrix(
    three_column_domain_template(p["LSM_ratios"]..., 
                              #
                              hole_ratio1=p["hole_ratio1"], hole_ratio2=p["hole_ratio2"], hole_ratio3=p["hole_ratio3"],
                              #                              
                              positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
                              #
                              column_width = p["column_width"],
                              #
                              height = p["height"]
    )
  )
end

#
#
function for_each_prms_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else      
      if typeof(prms_lists[active_idx]) <: Array
        list_for_iteration = prms_lists[active_idx]
      else  
        list_for_iteration = [prms_lists[active_idx]]
      end
      for parameter in list_for_iteration
        recursive_call(push!(deepcopy(output_set),parameter), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end




function par_study(
  # if value is array, the parametric study will iterate over it  
  input_prms_dict = Dict();
  save_to_file = ""
  )
  
  function extract_physical_parameters(par_study_prms)
    parameters = []
    for physical_prms_key in filter( 
          key -> length(key) > 2 && key[1:2] == "p.",
          collect(keys(par_study_prms))
        )
      push!(parameters, physical_prms_key[3:end] => par_study_prms[physical_prms_key])
    end   
    return parameters
  end

  if typeof(input_prms_dict) <: Vector
    input_prms_dict = Dict(input_prms_dict)
  end
  
  
  par_study_prms = Dict{String, Any}(
    "trials_count" => 2.0,
  )
  
  if !haskey(input_prms_dict, "matrix_template")
    merge!(par_study_prms,
      Dict(
        "hole_ratio" => 0.2,
        "LSM_ratio" => 0.2,    
        "dimensions" => (5, 5),
        "matrix_template" => homogenous_matrix
      )
    )
  end
  
  merge!(par_study_prms, input_prms_dict)

  output_data_frame = DataFrame()
  
  function execute_for_specific_values(specific_values)
    local_par_study_prms = Dict(collect(keys(par_study_prms)) .=> specific_values)     
    local_parameters = extract_physical_parameters(local_par_study_prms)
      
    for i in 1:par_study_prms["trials_count"]
      R, R_pol, C_pol = convert.(Float64,
        image_to_EIS(
                      local_par_study_prms["matrix_template"](local_par_study_prms),
                      local_parameters,
                      #
                      return_R_RC=true, 
                      TPE_warning=false,                      
                      pyplot=false,
        )
      )
      
      append!(output_data_frame, 
        merge(
          local_par_study_prms,
          Dict("trial_number" => i, "R" => R, "R_pol" => R_pol, "C_pol" => C_pol)
        )
      )
    end
  end
  
  for_each_prms_in_prms_lists(
    collect(values(par_study_prms)), 
    execute_for_specific_values
  )
  
  @show output_data_frame
  if save_to_file != ""
    CSV.write(save_to_file, output_data_frame)
  end
  return output_data_frame
end  

# par_study
function OLD_par_study(;
                    LSM_ratio_list = collect( 0 : 0.2 : 1.0),
                    hole_ratio = 0.2,
                    dimensions = (20,20),
                    matrix_template = default_matrix_template,                    
                    repeted_trials = 10,
                    parameters = []
  )
  
  R_list = []
  R_pol_list = []
  C_pol_list = []
  
  for LSM_ratio in LSM_ratio_list
    @show LSM_ratio
    
    res = []    
    for i in 1:repeted_trials
      push!(res,
        image_to_EIS( matrix_template(LSM_ratio, hole_ratio, dimensions),
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








