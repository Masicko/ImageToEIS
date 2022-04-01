function nyquistPlot(Z; fignum=1, label="")
    s = subplot(111)
    title("Nyquist plot")
    xlabel("Re\$(Z) \\ [\\Omega]\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    PyPlot.plot(real(Z), -imag(Z), label=label, "-x")
    if label!=""
      legend()
    end
    grid(true)
    s.set_aspect(1.0)  
end

function get_RC_prms_from_Z(omega, Z)
  a = real(Z)
  b = imag(Z)
  return (
    (a*a + b*b)/a,
    - b/(omega*(a*a + b*b))
    )
end

function get_R_RC_prms_from_2Z(omegas, Zs)
  a1, b1 = real(Zs[1]), imag(Zs[1])
  a2, b2 = real(Zs[2]), imag(Zs[2])
  
  tol = 1.0e-4
  if abs(a1 - a2) < tol*(a1+a2)/2
    R_ohm = a1
    R = 0
    C = 0
  else
        
    A = a1-a2
    B = ((a1^2 + b1^2) - (a2^2 + b2^2))
    C = a2*(a1^2 + b1^2) - a1*(a2^2 + b2^2)
        
    D = B^2 - 4*A*C
    
    s1 = (-B + sqrt(D))/   
            (2*A)
    
    R_ohm = -s1  
    R, C = get_RC_prms_from_Z(omegas[1], Zs[1] - R_ohm)
    #get_RC_prms_from_Z(omegas[2], Zs[2] - R_ohm)
  end
  
  return R_ohm, R, C
end

function get_Z_from_R_RC(f_range, R_ohm, R, C)
  [R_ohm + R/(1 + im*2*pi*f*R*C) for f in f_range]
end
