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
                      ["R_pol_YSZ"=> 0.0],
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



function plot_par_study_results(x, R, Rp, Cp)
  figure(5)
  title("R")
  plot(x, [sum(R[n]) for n in 1:length(R)])
  
  figure(6)
  title("R_pol")
  plot(x, [sum(R[n]) for n in 1:length(R)])
  
  figure(7)
  title("C_pol")
  plot(x, [sum(R[n]) for n in 1:length(R)])
end







