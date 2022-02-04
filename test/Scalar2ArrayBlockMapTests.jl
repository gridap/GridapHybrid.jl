module Scalar2ArrayBlockMapTests
   using Test
   using Gridap
   using GridapHybrid

   bs = [8,16]
   A = rand(sum(bs),sum(bs))
   b = rand(sum(bs))

   k=Scalar2ArrayBlockMap()
   cache=Gridap.Arrays.return_cache(k,(A,b),bs)
   Ab_bb =Gridap.Arrays.evaluate!(cache,k,(A,b),bs)

   Ab = Ab_bb[1]
   bb = Ab_bb[2]
   @test Ab.array[1,1]==A[1:8,1:8]
   @test Ab.array[2,1]==A[9:24,1:8]
   @test Ab.array[1,2]==A[1:8,9:24]
   @test Ab.array[2,2]==A[9:24,9:24]
end
