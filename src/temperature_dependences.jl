#using Interpolations


function TI(label)
  if label == "R_LSM"
    return CubicSplineInterpolation(      
        600 : 50 : 800
      , 1 ./reverse([
      286.138977466499
      205.521051620173
      150.497871331355
      124.493586820521
      107.754096313306
      ])
    )
  elseif label == "R_YSZ"
    return CubicSplineInterpolation(      
        600 : 50 : 800
      , 1 ./reverse([  
        0.044905404632136
        0.03005779931599
        0.01725995101764
        0.010115351717062        
        0.005578168458813
      ])
    )
 end
end

TI(label, T) = TI(label)(T)
