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

- *R_YSZ* : ohmic resistance of YSZ
- *C_pol_YSZ* : polarization capacitance on YSZ side
- *R_pol_YSZ* : polarization resistance on YSZ side
- *R_LSM* : ohmic resistance of LSM
- *C_pol_LSM* : polarization capacitance on LSM side
- *R_pol_LSM* : polarization resistance on LSM side
- *R_hole* : resistance of hole (is set to very high number as default)


Geometry specific behavior is simulated via interactions between pixels by the manner described in the following picture.

!["Particle impedance scheme"](images/scheme.png?raw=true )


The each impedance *Z* is spedified using the information about starting material *M1* and the ending material *M2* (see the picture) in the following manner:

!["Interaction scheme"](images/scheme_interaction.png?raw=true )

- *M1 = YSZ*
  - *M2 = YSZ*  => *Z = R_YSZ/2* 
  - *M2 = LSM*  => *Z = R_YSZ/2 + Z_RC(R_pol_YSZ, C_pol_YSZ)* (including the double-layer on YSZ side)
  - *M2 = hole* => *Z = R_YSZ/2*
- *M1 = LSM*
  - *M2 = YSZ*  => *Z = R_LSM/2 + Z_RC(R_pol_LSM, C_pol_LSM)* (including the double-layer on LSM side)
  - *M2 = LSM*  => *Z = R_LSM/2*
  - *M2 = hole* => *Z = R_LSM/2*
- *M1 = hole* => *Z = R_hole/2*

The factor 1/2 is there so a total resistance through the whole (e. g.) YSZ particle is *R_YSZ*.

## Installation
The package can be then installed via 
```julialang
] add https://github.com/Masicko/ImageToEIS.jl
```


## Usage

Before using the package, you have to execute

```julialang
using ImageToEIS
```

### Basics

Supposing we have either 

- `material_matrix = [1 1 1; 0 1 2]`(with values in {0,1,2}) or 
- bitmap image *my_image.png* (with colors {yellow, black, white}). 

In addition, we can specify parameters as a set of pairs. The following are the default parameters:

```julialang
physical_parameters = [ "R_YSZ" => 100, 
                        "R_pol_YSZ" => 0, 
                        "C_pol_YSZ" => 0.001, 
                        #
                        "R_LSM" => 1, 
                        "R_pol_LSM" => 40, 
                        "C_pol_LSM" => 0.005, 
                        #
                        "R_hole" => 1000000]
```

If less parameters are specified, the others are supposed to be default, i.e.

```julialang
physical_parameters = ["R_YSZ" => 73]
```

The core function is

```julialang
f_list, Z_list = image_to_EIS(material_matrix, physical_parameters)
```

or using path to image file

```julialang
f_list, Z_list = image_to_EIS("images/geometry.png", physical_parameters)
```

or specifying parameters in a function call 

```julialang
f_list, Z_list = image_to_EIS("images/geometry.png", ["R_YSZ" => 73])
```

or without specifying parameters


```julialang
f_list, Z_list = image_to_EIS("images/geometry.png")
```

which returns (by default) frequencies `f_list` for which impedances `Z_list` are computed.

### Additinal options

Practically useful keyword parameters are

- `f_list = [1, 10, 100]` : specification of array of frequencies for which EIS simulation will run. Good format is `= [2.0^n for n in (-5 : 0.5 : 15)]`.
  - `= "two_point_extrapolation"` : the simulation is run only for `f_list = [0.001, 1000]` yielding two impedances, 
      R-RC circuit is fitted to the two computed impedances. The output Z_list is computed using this R-RC circuit for 
      frequencies in TPE_f_list 
  - default value is `= "TPE"` : which is a shortcut for "two_point_extrapolation" with the same meaning
- `TPE_f_list = [2.0^n for n in (-5 : 0.5 : 15)]` 
- `pyplot = true` : if *false*, no Nyquist plot is plotted
- `return_R_RC = false` :
  - if `= true` : the output of function `image_to_EIS` is a tripple (R_ohm, R_pol, C_pol) from R-RC circuit
  - if `= false` : the output is a tuple `(f_list, Z_list)`
- `export_z_file = ""` : decides whether a standard file for z_view is exported
  - default value is `= ""`, which means *do nothing*
  - if `= "some_file.z"` : exports to this file
  - if `= "!use_file_name"` : this option is valid only when the function `image_to_EIS` was **called with a path of image**, e. g. "images/geometry.png"
  and it means that z_file will have a form "images/geometry.z", i. e. changes only the extension to ".z"
- `save_also_image = false`
  - if `= true` : if `export_z_file != ""` and `image_to_EIS` was **called with a path of image**, than the input image is copied with a name of `export_z_file` but with the extension of *input_image*

Advanced keyword parameters are 

- `complex_type = ComplexF64` : changes the data type in which the impedance calculation is performed
- `iterative_solver = false` : 
  - if `= false` : the system of equations is solved by a direct solver
  - if `= true` : the system is solved by iterative solver using Biconjugate gradient stabilized method


### Example

```julialang
image_to_EIS(   [1 1 1; 0 1 2], 
                ["R_YSZ" => 73],
                #
                export_z_file="test.z", 
                return_R_RC=true,
                save_also_image=true
                )
```

















