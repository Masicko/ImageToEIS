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


function template_par_study_three_domain()
  run_par_study(
    par_study_prms_dict = Dict(
                    "matrix_template" => ImageToEIS.three_column_domain_LSM_ratios,
                #
                    "trials_count" => 1,
                    #
                    #
                    "LSM_ratios" => [                                     
                                    INPUT1 .* INPUT2 for INPUT1 in collect(0.0 : 0.5 : 0.0), INPUT2 in [(1, 1, 1), (0, 1, 1), (0, 0, 1)]
                                    ],
                    #
                    "hole_ratio" => collect(0.0 : 0.05 : 0.1),
                    #
                    "positions_of_contacts" => (2, 7),
                    "height_of_contacts" => 1,
                    #
                    "column_width" => 1,
                    "height" => 10
    ),
    scripted_prms_names = [
                    "LSM_ratios", 
                    "hole_ratio" => ["hole_ratio1", "hole_ratio2", "hole_ratio3"]
    ],
    save_to_file_prefix = "3_domain_",
    direct = false,
    shell_command = "echo"
  )
end


function template_par_study_homogenous_matrix()
  run_par_study(
    par_study_prms_dict = Dict(
                            "matrix_template" => ImageToEIS.homogenous_matrix,
                            #
                            "trials_count" => 2,
                            #
                            "LSM_ratio" => collect(0.0 : 0.5 : 1.0),                            
                            "hole_ratio" => collect(0.0 : 0.5 : 0.5),
                            #
                            "dimensions" => (5,5),                                                        
                        ), 
    scripted_prms_names = ["LSM_ratio", "hole_ratio"],
    save_to_file_prefix = "homog_",
    direct = false
  )
end

function run_par_study(;shell_command="echo",
                        script_file="run_EX3_ysz_fitting_ImageToEIS.jl",
                        mode="go!", direct=false,
                        par_study_prms_dict::Dict,
                        scripted_prms_names::Array,
                        save_to_file_prefix = "default_"
                      )
  
  dict_list = []
  prm_name = "hole_ratio"
  
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
    save_to_file = save_to_file_prefix*"$(prms_tuple)"
    #
    if direct
      par_study(DICT, save_to_file=save_to_file)
    else
      if mode == "go!"
        run(`$(shell_command) $(script_file) $(DICT) $(save_to_file)`)
      end
    end            
  end
  
  for_each_prms_in_prms_lists(scripted_prms_lists, run_par_study_for_one_dict)
  println("Done :) ")
  return
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

function evaluate_DFs(dir, y_axis_labels)  
  total_DF = DataFrame()
  for file_name in readdir(dir)
    actual_DF = DataFrame(CSV.File(dir*"/$(file_name)"))
    total_DF = vcat(total_DF, actual_DF)
  end
  
  res_DF = DataFrame()
  if total_DF[!, :matrix_template][1] == "three_column_domain_LSM_ratios"        
    for subDF in groupby(total_DF, :LSM_ratios)
      LSM_tuple = eval(Meta.parse(subDF[!, :LSM_ratios][1]))
      LSM_ratio = LSM_tuple[3]
      if LSM_tuple[2] == 0.0 && LSM_tuple[3] != 0.0
          cell_type = 3
      elseif LSM_tuple[1] == 0.0 && LSM_tuple[2] != 0.0
          cell_type = 2
      else
          cell_type = 1
      end
      for primitive_DF in groupby(subDF, :hole_ratio1)                
        append!(
          res_DF, 
          Dict(
            :LSM_ratio => LSM_ratio,
            :cell_type => cell_type,
            :hole_ratio => primitive_DF[!, :hole_ratio1][1],
            [y_axis => sum(primitive_DF[!, y_axis])/length(primitive_DF[!, y_axis]) for y_axis in y_axis_labels]...
          )
        )
      end            
    end    
  elseif total_DF[!, :matrix_template][1] == "homogenous_matrix"    
    for subDF in groupby(total_DF, :LSM_ratio)      
      for primitive_DF in groupby(subDF, :hole_ratio)
        append!(
          res_DF, 
          Dict(
            :LSM_ratio => subDF[!, :LSM_ratio][1],                    
            :hole_ratio => primitive_DF[!, :hole_ratio][1],
            [y_axis => sum(primitive_DF[!, y_axis])/length(primitive_DF[!, y_axis]) for y_axis in y_axis_labels]...
          )
        )
      end
    end
  else 
    println("ERROR: unknown \"matrix_template\": $(total_DF[!, :matrix_template][1])")
    throw(Exception)
  end
  return res_DF
end 

function show_plots(x_axis, prms_choice, dir="snehurka/par_study/")  
  processed_df = evaluate_DFs(dir, [:R, :R_pol, :C_pol])
  #return processed_df
  
  primitive_DF = subset(processed_df, 
          [Symbol(pair[1]) => ByRow(==(pair[2])) for pair in prms_choice]...  
  )
  sort!(primitive_DF, x_axis)
  plot_par_study_results(
    primitive_DF[!, x_axis], 
    primitive_DF[!, :R], 
    primitive_DF[!, :R_pol],
    primitive_DF[!, :C_pol],
    ""
  )
  return
end
















