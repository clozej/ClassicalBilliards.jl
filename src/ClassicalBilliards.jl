module ClassicalBilliards

using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations
using ForwardDiff
using Elliptic
using Roots
using BilliardGeometry
#using Makie
#using SavingPlotting
include("particles/particles.jl")
include("dynamics/dynamics.jl")

end
