module DarcyHDGTestsSeq
using PartitionedArrays
include("../DarcyHDGTests.jl")
prun(DarcyHDGTests.main,sequential,(1,1,1))
prun(DarcyHDGTests.main,sequential,(2,2,2))
end # module
