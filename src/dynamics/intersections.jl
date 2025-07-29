function find_domain_exit_curves(particle::P, domain::D, dt) where {P<:AbsParticle, D<:AbsSimpleDomain}
    local_time = 0.0 #time since last bounce
    pt = propagate(particle, local_time + dt) #particle position after dt
    local_time += dt #increase time by dt
    d = is_inside(domain, pt) #check if particle is inside all curves
    while all(d) #while inside all curves
        pt = propagate(particle, local_time + dt)
        d = is_inside(domain, pt)
        local_time += dt
    end
    return d, local_time
end


function find_intersection_times(particle, curve, exit_time)
    #time since last bounce
    let
        r(t) = domain_fun(curve, propagate(particle, t))
        t = find_zeros(r, (0.0+1e-14, exit_time)) #10.0*eps(exit_time)
        return t  #returnes all intersections
    end
end 


function find_intersection_times(particle::P, line::C, exit_time) where {P<:PointParticle, C<:LineSegment}
    #time since last bounce
    let
        pt1 = line.pt1
        pt0 = line.pt0
        r = particle.r
        v = particle.v
        m = pt1 .- pt0
        return (pt1[1]*pt0[2] - pt0[1]*pt1[2] + m[2]*r[1] - m[1]*r[2])/(m[1]*v[2]-m[2]*v[1])
    end
end

function find_intersection_times(particle::P, circle::C, exit_time) where {P<:PointParticle, C<:CircleSegment}
    #time since last bounce
    let R = circle.radius, cent = circle.center, r = particle.r, v = particle.v
        dx = r[1] - cent[1]
        dy = r[2] - cent[2]
        a = v[1]*v[1] + v[2]*v[2]
        b = 2*(dx*v[1] + dy*v[2])
        c = dx*dx + dy*dy - R*R
        d = sqrt(b*b - 4*a*c)
        if isapprox(c, zero(c))
            t =  Vector{Float64}([-b/a])
        else
            t = Vector{Float64}([(d-b)/(2*a)])
        end
        return t
    end
end

function find_intersection(particle::P, domain::D; dt = 0.1) where {P<:AbsParticle, D<:AbsSimpleDomain}
    #time since last bounce
    boundary = domain.boundary
    N =  length(boundary)
    inside, approx_exit_time = find_domain_exit_curves(particle, domain, dt) #d gives the curves that were exited
    
    exit_curves = boundary[.~inside]
    crv_idx = collect(1:N)[.~inside]
    idx = crv_idx[1]
    exit_time = approx_exit_time
    for (i,crv) in enumerate(exit_curves)
        times = find_intersection_times(particle, crv, approx_exit_time)
        if  isempty(times)
            #println("Trajectory is broken at time $(particle.time)")
            continue
        else
            t = times[1]
            if t < exit_time
                exit_time = t
                idx = crv_idx[i]
            end
        end    
    end 
    return exit_time, idx
end

function line_polar(r,v,theta,center)
    pt = Translation(-center)(r)
    type = eltype(v)
    if v[1] == zero(type)
        return @. pt[1] / cos(theta)
    end
    if v[2] == zero(type)
        return @. pt[2] / sin(theta)
    end
    return @. (v[2]*pt[1] - v[1]*pt[2])/(v[2]*cos(theta) - v[1]*sin(theta))
end

function determine_brackets(r,v,center; eps=1e-12) 
    pt = Translation(-center)(r)
    dir = cross(pt,v)
    theta0 = rem2pi(atan(pt[2],pt[1]), RoundNearest)
    pole = atan(v[2],v[1])
    poles = rem2pi.([-pi, pole-pi, pole, pole+pi,  pi], RoundNearest)
    poles = sort(poles)
    unique = [true for p in poles]
    for i in 1:(length(poles)-1)
        if isapprox(poles[i],poles[i+1])
            unique[i] = false
        end
    end
    poles = poles[unique]
    #brackets = [(poles[1],poles[2]),(poles[2],poles[3]),(poles[3],poles[4])]
    if poles[1] <= theta0 <= poles[2]
        if dir > zero(dir)
            brackets =  [(theta0 + eps, poles[2])]
        else
            brackets =  [(-1.0*pi, theta0 - eps), (poles[3],1.0*pi)]
        end
    elseif poles[2] <= theta0 <= poles[3]
        if dir > zero(dir)
            brackets =  [(theta0 + eps, poles[3])]
        else
            brackets =  [(poles[2], theta0 - eps)]
        end
    elseif poles[3] <= theta0 <= poles[4]
        if dir > zero(dir)
            brackets = [(theta0 + eps,1.0*pi),(-1.0*pi, poles[2])]
        else
            brackets = [(poles[3], theta0 - eps)]
        end
    end
    return brackets
end


function find_intersection_angles(particle::PointParticle{T}, curve::C ; eps=1e-12) where {T, C<:AbsPolarCurve}
    #time since last bounce
    let  r = particle.r, v = particle.v, center = curve.center
        brackets = determine_brackets(r,v,center; eps)
        fun(theta) = polar_radius(curve,theta) - line_polar(r,v,theta,center)
        angles = Vector{eltype(r)}()
        for b in brackets
            angles1 = find_zeros(fun, b) #10.0*eps(exit_time)
            append!(angles, angles1)
        end
        return angles  #returns all angles
    end
end 

function find_intersection_times(particle::PointParticle{T}, curve::C, exit_time) where {T, C<:AbsPolarCurve}
    #time since last bounce
    let  r = particle.r, v = particle.v, center = curve.center
        angles = find_intersection_angles(particle, curve)
        pt0 = particle.r
        #pts_polar = [Polar(line_polar(r,v,phi,center),phi) for phi in angles]
        #pts_polar = [Polar(polar_radius(curve,phi),phi) for phi in angles]
        pts =[Translation(center)(CartesianFromPolar()(Polar(polar_radius(curve,phi),phi))) for phi in angles]
        #pts = line_polar.(r,v,angles)
        speed = norm(particle.v)
        return sort([hypot((pt-pt0)...)/speed for pt in pts])
    end
end

function find_intersection(particle::P, domain::D; dt = 0.1) where {P<:PointParticle, D<:PolarDomain}
    #time since last bounce
    boundary = domain.boundary
    times = find_intersection_times(particle, boundary[1], 100.0)
    if isempty(times)
        bounce_time = Inf64
    else
        bounce_time = times[1]
    end
    idx  = 1
    for i in 2:length(boundary)
        crv = boundary[i]
        times = find_intersection_times(particle, crv, 100.0)
        t = times[1]
        if 1e-14 < t < bounce_time
            bounce_time = t
            idx = i
        end   
    end 
    return bounce_time, idx
end