using Test
using ArrayInterface
using DynamicQuantities
using LinearAlgebra

@testset "ArrayInterfaceDynamicQuantitiesExt loads" begin
    ext = Base.get_extension(ArrayInterface, :ArrayInterfaceDynamicQuantitiesExt)
    @test ext !== nothing
end

@testset "Factorization instances for DynamicQuantities quantities" begin
    q = 1.0u"m"
    A = [q q; q q]

    luinst = ArrayInterface.lu_instance(A)
    @test luinst isa LinearAlgebra.LU
    @test size(luinst.factors) == (0, 0)

    qrinst = ArrayInterface.qr_instance(A)
    @test qrinst isa LinearAlgebra.QRCompactWY
    @test size(qrinst.factors) == (0, 0)

    svdinst = ArrayInterface.svd_instance(A)
    @test svdinst isa LinearAlgebra.SVD
    @test size(svdinst.U) == (0, 0)
    @test length(svdinst.S) == 0
    @test size(svdinst.V) == (0, 0)
end
