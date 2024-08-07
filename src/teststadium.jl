using Revise
using ClassicalBilliards
#include("plottingmakie.jl")
using GLMakie
using StaticArrays
using LinearAlgebra
using Roots
using BenchmarkTools

using LaTeXStrings
using MathTeXEngine # required for texfont

textheme = Theme(fontsize = 8,
            fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)),
            Axis = (;xgridvisible = false,
            ygridvisible = false, yticksize=3,xticksize=3,yticklabelpad=2.0,xticklabelpad=2.0, xlabelpadding = 0.0),
            #lines = (;grid=false)
                        
            )
set_theme!(textheme)

stad = Stadium(0.5)

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_billiard!(ax,stad; plot_virtual = true)
display(f)
#trajectory(p,stad,100)


#scatter!(ax, p.r[1], p.r[2])
#arrows!(ax, [p.r[1]], [p.r[2]], [p.v[1]], [p.v[2]]; lengthscale = 0.05)
p = Particle(1.0,0.6,-1.0,1.2;subdomain=1)
T = 20000

f = Figure(size=(1000,1000))
ax = Axis(f[1,1])
plot_trajectory!(ax, p, stad, T; dt = 0.01)
display(f)


p = Particle(1.0,0.6,-1.0,1.2;subdomain=1)
@btime trajectory(p,stad,100000; dt = 0.01)


#type = typeof(1.0)
#typeof(type(pi))
#p = Particle(1.0,0.6,-1.0,1.2;subdomain=1)
#pts1,vel1,ts = trajectory(p,stad,100000; dt = 0.01)

#ts

#append!(pts,pts1)
