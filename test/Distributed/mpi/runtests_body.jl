module NP4
# All test running on 4 procs here

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI
using GridapPETSc

include("../DarcyHDGTests.jl")

if ! MPI.Initialized()
  MPI.Init()
end

function all_tests(parts)
  display(parts)
  t = PArrays.PTimer(parts,verbose=true)
  PArrays.tic!(t)
  options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  GridapPETSc.with(args=split(options)) do
     DarcyHDGTests.main(parts)
  end
  PArrays.toc!(t,"DarcyHDGTests")
  display(t)
end

if MPI.Comm_size(MPI.COMM_WORLD) == 8
  prun(all_tests,mpi,(2,2,2))
elseif MPI.Comm_size(MPI.COMM_WORLD) == 1
  prun(all_tests,mpi,(1,1,1))
else
  MPI.Abort(MPI.COMM_WORLD,0)
end

end #module
