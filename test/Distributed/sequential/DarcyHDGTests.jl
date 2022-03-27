module DarcyHDGTestsSeq
using PartitionedArrays
include("../DarcyHDGTests.jl")
prun(DarcyHDGTests.main,sequential,(1,1))
end # module
