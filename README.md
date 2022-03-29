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

Note that `"R_pol_YSZ" => 0` which means no second arc should appear in Nyquist diagram and so *two point extrapolation* (see below in section Additional options) is viable as resulting to a dramatic speedup. If less parameters are specified, the others are supposed to be default, i.e.

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
      frequencies in TPE_f_list. When using this option, a warning is printed to output.
  - default value is `= "TPE"` : which is a shortcut for "two_point_extrapolation" with the same meaning
- `TPE_f_list = [2.0^n for n in (-5 : 0.5 : 15)]` 
- `TPE_warning = true` : if true and `"TPE" = true`, a warning about extrapolation is printed to standard output. 
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
- `store_R_RC` : turns on the evaluation of R_RC element from two points of computed impedance and append the output to a specified file
  - if `= ""` : means *do nothing* (which is default)
  - if `= "storage.txt"` : append a line formated as <dateTtime> tab <input_file_name> tab <R_ohm> tab <R_pol> tab <C_pol> . If the function was called with a matrix (not an path to file), `"<matrix_input>"` is written instead of <input_file_name>.
  - if `= "storage.csv"` : csv extension works too
- `return_specific_impedancde' = true` : returns specific impedance in *ohms x cm*. Also values *R_RC* will be returned in specific units.

Advanced keyword parameters are 

- `complex_type = ComplexF64` : changes the data type in which the impedance calculation is performed
- `iterative_solver = false` : 
  - if `= false` : the system of equations is solved by a direct solver
  - if `= true` : the system is solved by iterative solver using Biconjugate gradient stabilized method


### Real Example

```julialang
image_to_EIS(   [1 1 1; 0 1 2], 
                ["R_YSZ" => 73],
                #
                export_z_file="test.z", 
                return_R_RC=true,
                save_also_image=true,
                store_R_RC = ""
                )
```


## Generate random domain with structure

Automated generating of random structure is essential for statistical testing of system behavior. There are a few helping features using two main parameters. 

- `hole_ratio` in [0, 1] : the ratio of holes over the total points (material points + holes) in the picture. 
- `LSM_ratio` in [0, 1] : probability that the material point will be LSM.

### Simple homogenous matrix
The simplest example is a homogenous domain of `dimensions = (m, n)`, where *m* is a number of rows and *n* a number of columns. Matrix of this type can be constructed via

```julialang
homogenous_matrix = generate_matrix(dimensions, hole_ratio, LSM_ratio)
```

### Structure using multiple submatrices
For more complicated domains composed of several different homogenous subdomains, there is a possibility to construct appropriate matrix. Suppose we want to construct *m* x *n* matrix consisting of 2 different submatrices. First, a list of submatrices must be created such that one submatrix is represented by its *location* (left upper corner and right lower corner) in the resulting matrix and *hole_ratio* and *LSM_ratio*. 

```julialang
    # the structure for each submatrix in the list is >>
    # 
    # [left upper coord, right lower coord,  hole_ratio,   LSM_ratio]
    
submatrix_list = [        
    [(1,1),             (10, 5),            0.2,          0.5],   # this is the first submatrix
    [(1,6),             (10, 10),           0.0,          1.0]    # this is the second submatrix
]

two_subdomain_matrix = generate_matrix(submatrix_list)
```

There can be more subdomains in `submatrix_list`. Dimensions of the resulting `two_subdomain_matrix` is computed as en rectangular envelope of all *locations* in `submatrix_list`. The submatrices can overlap (in this case, the latter has priority in evaluation of matrix), but **every pixel must be covered** by the submatrices.

#### Three column domain template
There is a template using the upper structure of defying submatrices, which generates a domain of 3 columns with defined material specification (*hole_ratio* and *LSM_ratio*). In addition, there are contacts (of width 1 and optional height) on the left side providing an interesting distribution of electrical current through the system. The right side consists of a continuous one layer of LSM as a connection to conductive electrolyte. The obligate input parameters are 

- `LSM_ratio1`, `LSM_ratio2`, `LSM_ratio3`

Optional parameters (with default values) are 

- `hole_ratio1=0.5, hole_ratio2=0.5, hole_ratio3=0.5`
- `positions_of_contacts=[15, 50]` : starting row for each LSM contact 
- `height_of_contacts=5`
- `column_width=5`
- `height=70`

in this default case, the LSM contacts will be between [15, 20] and [50, 55] pixels of the first column while the whole matrix has height of 70 pixels. The submatrix list can be then obtained by function `three_column_domain_template` and then inserted to `generate_matrix`.

```julialang
LSM_ratio1 = 0.2
LSM_ratio2 = 1.0
LSM_ratio3 = 0.5

template_submatrix_list = three_column_domain_template(LSM_ratio1, LSM_ratio1, LSM_ratio1,
                                     #                              
                                     column_width = 10,
                                     hole_ratio1 = 0.0, hole_ratio3 = 0.2,
                                     positions_of_contacts=[20, 45], height_of_contacts=4
                                 )
                        
three_column_matrix = generate_matrix(template_submatrix_list)
```

### Material matrix visualization

A material matrix can be saved to a file with specified `path` using `matrix_to_file` function. For example

```julialang
matrix_to_file("images/three_column_domain.png", three_column_matrix)
```

!["Three column domain"](images/three_column_domain.png?raw=true )

## Parametric studies

















