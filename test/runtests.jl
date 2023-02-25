using ImageToEIS
using Test

function test_matrix_3x3x3_0_0()
	PRMS = ["R_YSZ" => 10.0, "R_LSM" => 1.0, "R_pol_LSM" => 66.0, "C_pol_LSM" => 0.005]
	(R_ohm, R, C) = image_to_EIS(generate_matrix((3,3,3), 0.0, 0.0), PRMS, return_R_RC=true)
    
    return prod(
        isapprox.(
            (R_ohm, R, C),
            (10.33333333333333, 44.00000000000011, 0.007499999999999982),
            rtol=1e-4
        )
    )
end

@testset "test_matrix_3x3x3_0_0" begin
    @test test_matrix_3x3x3_0_0()
end
