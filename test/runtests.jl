using ImageToEIS
using Test

function test_matrix_3x3x3_0_0()
	PRMS = ["R_YSZ" => 10.0, "R_LSM" => 1.0, "R_pol_LSM" => 66.0, "C_pol_LSM" => 0.005]
	res3D = image_to_EIS(generate_matrix((3,3,3), 0.0, 0.0), PRMS, return_R_RC=true)
    res2D = image_to_EIS(generate_matrix((3,3), 0.0, 0.0), PRMS, return_R_RC=true)
    # the coorect answer (R_ohm, R, C) 
    the_answer = (10.33333333333333, 44.00000000000011, 0.007499999999999982)
    
    pass_3D = prod(
        isapprox.(res3D, the_answer, rtol=1e-4)
    )

    pass_2D = prod(
        isapprox.(res2D, the_answer, rtol=1e-4)
    )
    return prod([pass_2D, pass_3D])
end

@testset "testin' matrix 3x3x3 and 3x3" begin
    @test test_matrix_3x3x3_0_0()

    # TODO 
    # test for
    #           - analyticaly correct (1,2)
    #           - bigger 2D (3,3) defined "random" matrix
end
