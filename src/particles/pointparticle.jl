mutable struct PointParticle{T}  <: AbsParticle where T<:Real
    r::SVector{2,T}
    v::SVector{2,T}
    m::T #mass of particle
    time::T #current time of the particle
    subdomain::Int64 #id of current subdomain
    subsegment::Int64 #id of last curve
    sym_sector::Int64 #id of current symmetry sector
    v_b::SVector{2,T} # velocity vector in local coordinate system on the boundary (v_t, v_n) v_t is tangential and v_n normal component.
end

PointParticle(x,y,vx,vy; m = 1.0, time = 0.0, subdomain=1, subsegment=0, sym_sector=1, v_b = SVector(vx,vy)) = PointParticle{typeof(x)}(SVector(x,y),SVector(vx,vy),m,time,subdomain,subsegment,sym_sector, v_b)
PointParticle(point,velocity; m = 1.0, time = 0.0, subdomain=1, subsegment=0, sym_sector=1, v_b = velocity) = PointParticle{eltype(point)}(point,velocity,m,time,subdomain,subsegment,sym_sector, v_b)

function propagate(part::P, t) where P<:PointParticle
    return @. part.r + part.v*t
end

function propagate(part::P, ts::AbstractArray) where P<:PointParticle
    return [@. part.r .+ part.v.*t for t in ts] 
end