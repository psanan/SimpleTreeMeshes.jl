using Test
using MPI

push!(LOAD_PATH, "../src")
using SimpleTreeMeshes

# This test obviously does not do much!
function test1()
  refinement_function(x, y, d, level) = level < 3
  MPI.Init()
  tree = CreateSimpleTreeMesh(refinement_function)
  return true
end

@testset "basic" begin
  @test test1()

end
