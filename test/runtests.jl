using Test

push!(LOAD_PATH, "../src")
using SimpleTreeMeshes

function test1()
  refinement_function(x, y, d, level) = level < 3
  tree = CreateSimpleTreeMesh(refinement_function)
  # ..
  return true
end

@testset "basic" begin
  @test test1()

end
