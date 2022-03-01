# ImageToEIS.jl

## Introduction

The goal is to simulate an impedance measurement of an electrochemical cell of defined structure consisting of YSZ, LSM and pores using electrical elements. 

### Geometry
The structure is defined by a material matrix which contains integers

- *0* : for a YSZ particle (yellow),
- *1* : for a LSM particle (black),
- *2* : for a hole (white).

Material matrix can be also defined by a raster image using colors specified above in brackes.

!["Specified cell geometry"](images/geometry.png?raw=true )

In addition, left and right sides of the matrix is grasped by 1 column of LSZ layer as current collectors. 

!["Cell geometry with electordes"](images/geometry_with_electrodes.png?raw=true )

### Electrochemistry
Physical simulation is done using standard electrical elements - R elements for ohmic resistance and RC element for a double-layer behavior on interfaces of LSM | YSZ. We have the following parameters:

- *R_YSZ* : resistance of YSZ
- *R_LSM* : resistance of LSM
- *R_hole* : resistance of hole (is set to very high number as default)
- *C_pol* : polarization capacitance
- *R_pol* : polarization resistance

Geometry specific behavior is simulated via interactions between pixels by the manner described in the following picture.

!["Particle impedance scheme"](images/scheme.png?raw=true )


The each impedance *Z* is spedified using the information about starting material *M1* and the ending material *M2* (see the picture) in the following manner:

!["Interaction scheme"](images/scheme_interaction.png?raw=true )

- *M1 = YSZ* => *Z1 = R_YSZ/2*    (factor 1/2 is there so a total resistance through the whole YSZ particle is *R_YSZ*)
- *M1 = hole* => *Z1 = R_hole/2*
- *M1 = LSM* 
  - *M2 = YSZ*  => *Z1 = R_LSM/2 + Z_RC(R_pol, C_pol)* (including the double-layer)
  - *M2 = LSM*  => *Z1 = R_LSM/2*
  - *M2 = hole* => *Z1 = R_LSM/2*


## Installation
The package can be then installed using 
```julialang
] add https://github.com/Masicko/ImageToEIS.jl
```


## Usage
Supposing we have either `material_matrix`(with values in {0,1,2}) of bitmap image *my_image* (with colors {yellow, black, white}). In addition, we can specify parameters via pairs. The following are the default parameters:

```julialang
physical_parameters = ["R_YSZ" => 100, "R_LSM" => 1, "R_pol" => 40, "C_pol" => 0.001, "R_dira" => 1000000]
```

If less parameters are specified, the others are supposed to be default, i.e.

```julialang
physical_parameters = ["R_YSZ" => 73]
```

The core function is
```julialang
f_range, Z_range = image_to_EIS(material_matrix, physical_parameters)
```
or

```julialang
f_range, Z_range = image_to_EIS(physical_parameters, path_to_file="images/geometry.png")
```

which returns frequencies `f_range` for which impedances `Z_range` are computed.





















