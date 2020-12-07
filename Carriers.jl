module Carriers
include("Constants.jl")
using .Constants

export Carrier, initializecarriers, flight, drift!, scatterk!, averagevelocity, averageenergy, population

mutable struct Carrier
    force::Array{Float64}
    valley_maxima::Array{Float64}
    energy::Float64
    k::Array{Float64}
    valley_indice::Int
    #Carrier(energy, k, valley_indice) = Carrier([], [],energy, k, valley_indice)
    #Carrier(force, valley_maxima, energy, k, valley_indice) = Carrier(force::Array{Float64}, valley_maxima::Array{Float64}, energy::Float64, k::Array{Float64}, valley_indice::Int)
end

function initializecarriers(force::Array{Float64}, valleymaxima::Array{Float64})
    C = Carrier[]
    r = zero(Float64)
    for i in 1:1000
        r = rand()
        ϕ = 2*π*r
        r = rand()
        cosθ = 1 - 2*r
        sinθ = sqrt(1 - cosθ^2)

        r = rand()
        energy = -1.5*kT*log(r)
        k = sqrt(2*mass_array[1]*charge)*sqrt(energy*(1 + alpha_array[1]*energy))/ħ
        kx = k*sinθ*cos(ϕ)
        ky = k*sinθ*sin(ϕ)
        kz = k*cosθ
        k = [kx, ky, kz]
        valley_indice = rand(1:3)

        push!(C,Carrier(force, valleymaxima, energy, k, valley_indice))
    end

    return C
end

#################################################
#update wavevector and energy as time is evolved#
#################################################
function drift!(c::Carrier, τ::Float64)
    c.k[1] = c.k[1] - charge*c.force[1]*τ/ħ
    c.k[2] = c.k[2] - charge*c.force[2]*τ/ħ
    c.k[3] = c.k[3] - charge*c.force[3]*τ/ħ
    k = sqrt(c.k[1]^2 + c.k[2]^2 + c.k[3]^2)
    γₖ = (k*ħ)^2/(2*charge*mass_array[c.valley_indice])
    c.energy = (2*γₖ)/(1+sqrt(1 + 4*alpha_array[c.valley_indice]*γₖ))
end

#determine flight time based on total scattering rate
function flight(c::Carrier)
    Γ_total = c.valley_maxima[c.valley_indice]
    #print("Γ_total: ",Γ_total)
    r = rand()
    time = -log(r)/Γ_total
    return time
end

#################################################################
#scatter the wavevector depending on scattering mechanism "mech"#
#################################################################
function scatterk!(c::Carrier, mech::Int, newenergy::Float64)
    valley = c.valley_indice
    if mech == 2 || mech == 3
        ξ = 0.1
        try
            ξ = 2*sqrt(c.energy*newenergy)/(sqrt(c.energy) - sqrt(newenergy))^2
        catch e
            println("Mech is: ",mech)
            println("New energy is: ",newenergy)
        end

        r = rand()
        cosθ = (1+ξ - (1 + 2*ξ)^r)/ξ
        r = rand()
        ϕ = 2*π*r

        kp = sqrt(2*mass_array[valley]*charge)*sqrt(newenergy*(1+alpha_array[valley]*newenergy))/ħ

        sinθ = sqrt(1 - cosθ^2)
        kxp = kp*sinθ*cos(ϕ)
        kyp = kp*sinθ*sin(ϕ)
        kzp = kp*cosθ

        kxy = sqrt(c.k[1]^2 + c.k[2]^2)
        k = sqrt(kxy^2 + c.k[3]^2)

        cosθ_0 = c.k[3]/k
        sinθ_0 = kxy/k
        sinϕ_0 = c.k[2]/kxy
        cosϕ_0 = c.k[1]/kxy

        c.k[1] = kxp*cosϕ_0*cosθ_0 - kyp*sinϕ_0 + kzp*cosϕ_0*sinθ_0
        c.k[2] = kxp*sinϕ_0*cosθ_0 - kyp*cosϕ_0 + kzp*sinϕ_0*sinθ_0
        c.k[3] = -kxp*sinθ_0 + kzp*cosθ_0
    elseif mech == 10
        c.energy = newenergy #do nothing to momentum
    else
        if mech == 4 || mech == 5
            valley = 1
        elseif mech == 6 || mech == 7
            valley = 2
        elseif mech == 8 || mech == 9
            valley = 3
        end

        r = rand()
        cosθ = 1 - 2*r
        r = rand()
        ϕ = 2*π*r

        k = sqrt(2*mass_array[valley]*charge)*sqrt(newenergy*(1+alpha_array[valley]*newenergy))/ħ
        sinθ = sqrt(1 - cosθ^2)
        c.k[1] = k*sinθ*cos(ϕ)
        c.k[2] = k*sinθ*sin(ϕ)
        c.k[3] = k*cosθ
    end

    c.valley_indice = valley
    c.energy = newenergy
end

##############################################################
#return average velocity of carrier array in specified valley#
##############################################################
function averagevelocity(c::Array{Carrier}, valley::Int)
    velocity = zeros(3)
    n = 0
    for i in axes(c,1)
        if c[i].valley_indice == valley
            n += 1
            velocity += c[i].k .* ħ/(mass_array[valley]*(1 + 2*alpha_array[valley]*c[i].energy))
        end
    end

    velocity = velocity./n
    for i in velocity
        if isnan(i)
            i = 0.0
        end
    end
    return velocity
end

############################################################
#return average energy of carrier array in specified valley#
############################################################
function averageenergy(c::Array{Carrier}, valley::Int)
    n = 0
    energy = 0.0
    for i in axes(c,1)
        if c[i].valley_indice == valley
            n += 1
            energy += c[i].energy
        end
    end

    if isnan(energy/n)
        return 0
    else
        return energy/n
    end
end

########################################################
#return population of carrier array in specified valley#
########################################################
function population(c::Array{Carrier}, valley::Int)
    n = 0
    for i in axes(c,1)
       if c[i].valley_indice == valley
            n += 1
        end
    end

    return n
end

end
