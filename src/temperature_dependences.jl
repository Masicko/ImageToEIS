function fixed_R_YSZ()
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

# temperature data are in degrees of Celsia
# resistanace data are in Ohm/cm
# ... a bit different interpretation of porosity
function TI_clank(label)
  if label == "R_LSM"
    
    return CubicSplineInterpolation(      
        600 : 50 : 800
      , 
      # fitted to 3D with "clankova" porozita      
      [0.004640194824206025, 0.0040162711411057365, 0.0033223061268364204, 0.0024328408017493925, 0.0017474026238125484]
      # fitted to 2D
      #[0.005764, 0.005022, 0.00413, 0.003021, 0.002185]      
    )
  elseif label == "R_YSZ"
    # ok data from pure YSZ (porosity = 0.0)
    return fixed_R_YSZ()
 end  
end

function TI(label)
  if label == "R_LSM"
    
    return CubicSplineInterpolation(      
        600 : 50 : 800
      , 
      # fitted to 3D
      [0.006418221532022811, 0.005556938763855908, 0.004598968863545539, 0.00337112453237228, 0.0024249253299838836]
      # fitted to 2D
      #[0.005764, 0.005022, 0.00413, 0.003021, 0.002185]      
    )
  elseif label == "R_YSZ"
    # ok data from pure YSZ (porosity = 0.0)
    return fixed_R_YSZ()
 end
end

function TI_2D(label)
  if label == "R_LSM"
    
    return CubicSplineInterpolation(      
        600 : 50 : 800
      , 
      # fitted to 3D
      #[0.006418221532022811, 0.005556938763855908, 0.004598968863545539, 0.00337112453237228, 0.0024249253299838836]
      # fitted to 2D
      [0.005764, 0.005022, 0.00413, 0.003021, 0.002185]      
    )
  elseif label == "R_YSZ"
    # ok data from pure YSZ (porosity = 0.0)
    return fixed_R_YSZ()
 end
end

function TI_por_LSM(T)
  #0.12 porosity LSM 
  return CubicSplineInterpolation(      
      600 : 50 : 800
    , 1 ./reverse([
    286.138977466499
    205.521051620173
    150.497871331355
    124.493586820521
    107.754096313306
    ])
  )(T)
end

TI(label, T) = TI(label)(T)
TI_2D(label, T) = TI_2D(label)(T)
TI_clank(label, T) = TI_clank(label)(T)
