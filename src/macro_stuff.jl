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
              par_study_prms["porosity"]
            ), 
            par_study_prms["dimensions"]...
  )
end

function three_column_domain_matrix(p::Dict)
  return generate_matrix(
    three_column_domain_template(p["LSM_ratio1"], p["LSM_ratio2"], p["LSM_ratio3"], 
                              #
                              porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
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
                              porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
                              #                              
                              positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
                              #
                              column_width = p["column_width"],
                              #
                              height = p["height"]
    )
  )
end

# function three_column_domain_LSM_ratios(p::Dict)  
#   return generate_matrix(
#     three_column_domain_template(p["LSM_ratios"]..., 
#                               #
#                               porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
#                               #                              
#                               positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
#                               #
#                               column_width = p["column_width"],
#                               #
#                               height = p["height"]
#     )
#   )
# end





















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
    "repetition_idx" => 1.0,
  )
  
  if !haskey(input_prms_dict, "matrix_template")
    merge!(par_study_prms,
      Dict(
        "porosity" => 0.2,
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
        Dict("R" => R, "R_pol" => R_pol, "C_pol" => C_pol)
      )
    )  
  end
  
  for_each_prms_in_prms_lists(
    collect(values(par_study_prms)), 
    execute_for_specific_values
  )
  
  #@show output_data_frame
  if save_to_file != ""
    slash_idxs = findall(x -> x in ['/'], save_to_file)
    if length(slash_idxs) > 0
      mkpath(save_to_file[1:slash_idxs[end]])
    end      
    CSV.write(save_to_file, output_data_frame)
  end
  return output_data_frame
end











































function template_par_study_three_domain()
  run_par_study(
    par_study_prms_dict = Dict(
                    "matrix_template" => ImageToEIS.three_column_domain_LSM_ratios,
                #
                    "repetition_idx" => collect(1:1),
                    #
                    #
                    "LSM_ratios" => [                                     
                                    INPUT1 .* INPUT2 for INPUT1 in collect(0.0 : 0.5 : 0.0), INPUT2 in [(1, 1, 1), (0, 1, 1), (0, 0, 1)]
                                    ],
                    #
                    "porosity" => collect(0.0 : 0.05 : 0.1),
                    #
                    "positions_of_contacts" => (2, 7),
                    "height_of_contacts" => 1,
                    #
                    "column_width" => 1,
                    "height" => 10
    ),
    scripted_prms_names = [
                    "LSM_ratios", 
                    "porosity" => ["porosity1", "porosity2", "porosity3"]
    ],
    save_to_file_prefix = "3_domain_",
    direct = false,
    shell_command = "echo"
  )
end


function template_par_study_homogenous_matrix()
  ImageToEIS.run_par_study(
    par_study_prms_dict = Dict(
                            "matrix_template" => ImageToEIS.homogenous_matrix,
                            #
                            "repetition_idx" => collect(1:1),
                            #
                            "LSM_ratio" => collect(0.0 : 0.5 : 1.0),                            
                            "porosity" => collect(0.0 : 0.5 : 0.5),
                            #
                            "dimensions" => (5,5),                                                        
                        ), 
    scripted_prms_names = ["LSM_ratio", "porosity"],
    save_to_file_prefix = "homog_",
    direct = false,
    shell_command = "echo"
  )
end

function template_par_study_temperature()
  for T in [600, 650, 700, 750, 800]
         ImageToEIS.run_par_study(
           par_study_prms_dict = Dict(
                                   "matrix_template" => ImageToEIS.homogenous_matrix,
                                   #
                                   "repetition_idx" => collect(1:1),
                                   #
                                   "LSM_ratio" => collect(0.44 : 0.5 : 0.44),                            
                                   "porosity" => collect(0.5 : 0.5 : 0.5),
                                   #
                                   "dimensions" => (4,4),
                                   #
                                   #
                                   "T" => T,
                                   "R_YSZ" => TI("R_YSZ", T),
                                   "R_LSM" => TI("R_LSM", T)
                               ), 
           scripted_prms_names = ["repetition_idx", "LSM_ratio", "porosity"],
           save_to_file_prefix = "por_study_$(T)_",
           direct = false,
           shell_command = "echo"
         )
       end
end

function construct_explicit_par_study_dict(par_study_prms_dict, parameter_dependency)
  
end

function run_par_study(;shell_command="echo",
                        script_file="run_EX3_ysz_fitting_ImageToEIS.jl",
                        mode="go!", direct=false,
                        par_study_prms_dict::Dict,
                        paramter_dependency::Expr=:(),
                        scripted_prms_names::Array,
                        save_to_file_prefix = "default_"
                      )
  
  
  
  scripted_prms_lists = [                          
                          typeof(prm_name) <: Pair ?                           
                            par_study_prms_dict[prm_name[1]] :
                            par_study_prms_dict[prm_name]                          
                          for prm_name in scripted_prms_names
                        ]
                        
  function run_par_study_for_one_dict(prms_tuple)                
    DICT = deepcopy(par_study_prms_dict)
    update_pairs = []
    for i in 1:length(scripted_prms_names)      
      if typeof(scripted_prms_names[i]) <: Pair        
        push!(update_pairs, [identifier => prms_tuple[i] for identifier in scripted_prms_names[i][2]]...)
        delete!(DICT, scripted_prms_names[i][1])
      else
        push!(update_pairs, scripted_prms_names[i] => prms_tuple[i])
      end
    end      
    DICT = merge(DICT, Dict(update_pairs))
    DICT_str = string(DICT)
    save_to_file = save_to_file_prefix*"$(prms_tuple).csv"
    #    
    if mode == "go!" || mode == "try_one_prms"
      if direct
        par_study(DICT, save_to_file=save_to_file)
      else      
        run(`$(shell_command) $(script_file) $(DICT_str) $(save_to_file)`)
      end
    elseif mode == "only_save_image"        
        @show DICT
        matrix_to_file("save_first_image.png", DICT["matrix_template"](DICT))
    end
  end
  
  if mode == "only_save_image" || mode == "try_one_prms"
    @show scripted_prms_lists
    first_prms_array = [list[Int(floor(end/2)+1)] for list in scripted_prms_lists]
    @show first_prms_array
    run_par_study_for_one_dict(first_prms_array)
  else
    for_each_prms_in_prms_lists(scripted_prms_lists, run_par_study_for_one_dict)
  end
  println("Done :) ")
  return
end


















































#### Evaluation

function test_get_processed_df()
  x_axis = "T"
  other_parameters = [
    "porosity" => [0.1, 0.2, 0.3],
    "LSM_ratio"
  ]
  return show_plots(x_axis, other_parameters, "snehurka/par_studies/por_study_50x50/")  
end


function collect_df_files_in_folder(dir)
  collected_df = DataFrame()
  for file_name in readdir(dir)
    actual_DF = DataFrame(CSV.File(dir*"/$(file_name)"))
    collected_df = vcat(collected_df, actual_DF)
  end
  return collected_df
end

function test_if_all_changing_prms_were_assigned(df, keys)
  for group in groupby(df, [keys..., "repetition_idx"])
    if size(group)[1] != 1
      return false
    end
  end
  return true
end

function get_grouped_processed_df(collected_df, x_axis, other_parameters; specific_symbol="", throw_exception=true)
  processed_df = deepcopy(collected_df)    
  
  keys_to_group = []
  for other_parameter in other_parameters  
    if typeof(other_parameter) <: Pair
      filter!(row -> row[other_parameter[1]] in other_parameter[2], processed_df)
      push!(keys_to_group, other_parameter[1])
    else
      push!(keys_to_group, other_parameter)
    end
  end
  
  gdf = groupby(processed_df, [keys_to_group..., x_axis]) 
  #
  if !test_if_all_changing_prms_were_assigned(processed_df, [keys_to_group..., x_axis])
    println("\n <<< ! ! ! ERROR! not all changing parameters were specified! >>> \n")
    if throw_exception
      return throw(Exception)
    end
  end
  #
  df_to_plot = combine(gdf, ["R", "R_pol", "C_pol"] .=> arr -> sum(arr)/length(arr), renamecols = false)
  return groupby(df_to_plot, keys_to_group)
end

function show_plots(x_axis, other_parameters, dir="snehurka/par_study/"; 
        specific_symbol="", 
        apply_func= x -> x, 
        throw_exception=true)
  collected_df = collect_df_files_in_folder(dir)
  #@show collected_df
  grouped_df = get_grouped_processed_df(collected_df, x_axis, other_parameters, specific_symbol=specific_symbol, throw_exception=throw_exception)
  
  for (key, sub_df) in pairs(grouped_df)
    legend_entry = "$(key)"[12:end-1]
    sort!(sub_df, x_axis)
    @show sub_df
    plot_par_study_results(
      sub_df[!, x_axis], 
      apply_func.(sub_df[!, :R]), 
      apply_func.(sub_df[!, :R_pol]),  
      apply_func.(sub_df[!, :C_pol]),
      label=legend_entry,
      x_axis_label=x_axis
    )    
  end
  return grouped_df
end



function plot_par_study_results(x, R, Rp, Cp; label="", x_axis_label)
  figure(5)
  #suptitle("YSZ-LSM-hole 50x50 random distribution test - repetitions for each conditions")
  
  subplot(221)
  @show R
  plot(x, [sum(R[n])/length(R[n]) for n in 1:length(R)], label=label, "-x")
  xlabel(x_axis_label)
  ylabel("R_ohm")
  legend()
  
  subplot(223)    
  plot(x, [sum(Rp[n]/length(Rp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel(x_axis_label)
  ylabel("R_pol")
  legend()
  
  subplot(222)  
  plot(x, [sum(Cp[n]/length(Cp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel(x_axis_label)
  ylabel("C_pol")
  legend()
  
  subplot(224)  
  plot(x, [sum(R[n])/length(R[n]) for n in 1:length(R)] .+ [sum(Rp[n]/length(Rp[n])) for n in 1:length(R)], label=label, "-x")
  xlabel(x_axis_label)
  ylabel("R_tot")
  legend()
end















