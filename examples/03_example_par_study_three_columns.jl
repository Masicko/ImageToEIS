module 03_example_par_study_three_columns



function LSM_ratios_as_tuple()
  return ImageToEIS.par_study(
                 Dict(
                     "matrix_template" => ImageToEIS.three_column_domain_LSM_ratios,
                     #              
                     "trials_count" => 20,
                     #
                     #
                     "LSM_ratios" => [
                                       (0.5, 0.5, 0.5), 
                                       (0.0, 0.5, 0.5),
                                       (0.0, 0.0, 0.5)
                                       ],
                     #
                     "porosity1" => 0.5,
                     "porosity2" => 0.5,
                     "porosity3" => 0.5,
                     #
                     "positions_of_contacts" => (10, 30),
                     "height_of_contacts" => 5,
                     #
                     "column_width" => 5,
                     "height" => 60
                 ),
                 save_to_file = "jojo"
       )
end


           
end # of module
