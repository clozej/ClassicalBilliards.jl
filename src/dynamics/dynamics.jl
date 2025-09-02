include("intersections.jl")
export find_domain_exit_curves, find_intersection_times, find_intersection_angles, find_intersection, line_polar, determine_brackets
include("colissionrules.jl")
export collision_rule!

function collision!(particle::P, curve::C, time) where {P<:PointParticle, C<:AbsCurve}
    r = propagate(particle, time) #colision point
    g = domain_gradient_vector(curve, r)
    n =  g./norm(g) 
    collision_rule!(particle, curve.bc, n) #modify velocity using rule corresponding to boundary
    particle.r = r
    particle.time += time
end

function iterate_bounce!(particle::P, billiard::B; dt = 1.0, full_domain=true) where {P<:AbsParticle, B<:AbsBilliard}
    if typeof(billiard.fundamental_domain) <: AbsSimpleDomain
        domain = billiard.fundamental_domain
    else
        domain = billiard.fundamental_domain.subdomains[particle.subdomain]   
    end 
    collision_time, idx = find_intersection(particle, domain; dt)
    crv = domain.boundary[idx]
    particle.subsegment = idx
    collision!(particle, crv, collision_time)
    bc_type = typeof(crv.bc)

    if full_domain
        if bc_type <: Transparent || bc_type  <: ReflectionSymmetry
            iterate_bounce!(particle, billiard; dt, full_domain)
        end
    else
        if bc_type <: Transparent
            iterate_bounce!(particle, billiard; dt, full_domain)
        end
    end
end

function trajectory(particle::P, billiard::B, T::Int; dt = 1.0, full_domain=true) where {P<:AbsParticle, B<:AbsBilliard}
    let p = particle
        pts = Vector{typeof(p.r)}(undef,T+1)
        vel = Vector{typeof(p.v)}(undef,T+1)
        ts = Vector{typeof(p.time)}(undef,T+1)
        #sym_ids = Vector{Int64}(undef,T+1)
        pts[1] = p.r
        vel[1] = p.v
        ts[1] = p.time
        #sym_ids[1] = p.sym_sector
        #println("0, r=$(p.r)")
        for i in 1:T
            iterate_bounce!(p, billiard; dt, full_domain)
            ##println("$i, r=$(p.r), dom_idx=$(p.subdomain)")
            if full_domain
                id = p.sym_sector
                if id  == 1
                    pts[i+1] = p.r
                    vel[i+1] = p.v
                else
                    sym =  billiard.symmetries[id-1]
                    pts[i+1] = apply_symmetry(sym, p.r)
                    vel[i+1] = apply_symmetry(sym, p.v)
                end
            else
                pts[i+1] = p.r
                vel[i+1] = p.v
            end
            ts[i+1] = p.time
            #sym_ids[i+1] = p.sym_sector
        end
        return pts, vel, ts
    end
end

function symbolic_trajectory(particle::P, billiard::B, T::Int; dt = 1.0, full_domain=true) where {P<:AbsParticle, B<:AbsBilliard}
    let p = particle
        symbol = SVector{3,Int64}([particle.subsegment,particle.subdomain,particle.sym_sector])
        sym_traj = [symbol]
        for i in 1:T
            iterate_bounce!(p, billiard; dt, full_domain)
            push!(sym_traj,  SVector{3,Int64}([particle.subsegment,particle.subdomain,particle.sym_sector]))
        end
        return sym_traj
    end
end

function pb_trajectory(particle::P, billiard::B, T::Int; dt = 1.0) where {P<:AbsParticle, B<:AbsBilliard}
    let p = particle
        pb_pts = Vector{PoincareBirkhoff}(undef,T+1)
        if p.subsegment == 0
            iterate_bounce!(p, billiard; dt, full_domain=true)
        end
        pb_pt = pb_coords(billiard, p.subsegment, p.subdomain, p.sym_sector, p.r, p.v)
        pb_pts[1] = pb_pt
        for i in 1:T
            iterate_bounce!(p, billiard; dt, full_domain=true)
            pb_pt = pb_coords(billiard, p.subsegment, p.subdomain, p.sym_sector, p.r, p.v)
            pb_pts[i+1] =  pb_pt
        end
        return pb_pts
    end
end
export trajectory, symbolic_trajectory, pb_trajectory, iterate_bounce!, colission!