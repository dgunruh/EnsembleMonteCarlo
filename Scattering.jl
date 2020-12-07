module Scattering
include("Constants.jl")
using .Constants

export createscatteringtable, choosemech, returnenergy, displaytable

#######################
#Perform initial setup#
#######################
#electron acoustic deformation potential
Ξ = [7.01*charge, 9.2*charge, 9.0*charge]

#longitudinal optical phonon energy in eV
ħω₀ = 0.03536

#longitudinal acoustic velocity
vs = 5.22e3

#Mass density
ρ = 5370

#Young's modulus
cₗ = ρ*vs^2

#matrix element squared for fermi's golden rule, use kT in J not eV
M2 = Ξ.^2.0 .*(kT*charge/cₗ)

#longitudinal optical phonon frequency
ω₀ = charge*ħω₀/ħ

#occupancy under bose-einstein statistics
n₀ = 1/(exp(ħω₀/kT) - 1)

#high-frequency dielectric constant
ϵ_infinity = 10.92

#low-frequency dielectric constant
ϵ_static = 12.9

#dielectric constant which combines the two limits. See Computational electronics pg. 119
ϵₚ = ϵ₀/(ϵ_infinity^-1.0 - ϵ_static^-1.0)

########################
#Intervalley parameters#
########################

#energy separation between valleys (eV)
Δ = [0 0.29 0.48; -0.29 0 0.19; -0.48 -0.19 0]

#intervalley phonon energy (eV)
E_matrix = [0 0.0278 0.0299; 0.0278 0.029 0.0293; 0.0299 0.0293 0.0299]

#intervalley deformation potential
Ξ_matrix = [0 1.8e10 10e10; 1.8e10 5e10 1e10; 10e10 1e10 10e10]

#total number of available final valleys for the carrier to scatter into
Z_matrix = [0 4 3; 1 3 3; 1 4 2]

#intervalley phonon frequency
ω_matrix = E_matrix .* charge/ħ
n₀_in = 1.0 ./ (exp.(E_matrix ./ kT) .- 1)

#set acoustic scattering rates for all particles based on energies
function acoustic(valley::Int, E::Float64)
    #calculate the half DOS function for nonparabolic bands
    N = ((2*mass_array[valley])^1.5 /(4* π^2 * ħ^3)) * sqrt(E*charge*(1 + alpha_array[valley]*E)) * (1 + 2*alpha_array[valley]*E)
    #now calculate transition rate using Fermi's golden rule
    Γ_ac = 2*π/ħ * M2[valley] * N
    return Γ_ac
end

#set polar optical phonon scattering rates for all particles based on energies
function polar(valley::Int, E::Float64)
    m = mass_array[valley]
    α = alpha_array[valley]
    prefactor = sqrt(m/2.0) * ω₀ * charge^2 /(4*π*ħ*ϵₚ)
    Γ_pol_em = 0.0

    #gamma factor for phonon emission
    Eₑ = E - ħω₀
    if Eₑ >= 0
        γₖₑ = Eₑ*( 1 + α*Eₑ)
        γₖ = E*(1 + α*E)
        A = (2*(1 + α*E)*(1 + α*Eₑ) + α*(γₖ + γₖₑ))^2
        B = -2*α*sqrt(γₖ*γₖₑ)*(4*(1 + α*E)*(1 + α*Eₑ) + α*(γₖ + γₖₑ))
        C = 4*(1 + α*E)*(1 + α*Eₑ) + (1 + 2*α*E)*(1 + 2*α*Eₑ)
        F = (A*log((sqrt(γₖ) + sqrt(γₖₑ))/(sqrt(γₖ) - sqrt(γₖₑ))) + B)/C

        Γ_pol_em = prefactor * (n₀ + 1) * (1 + 2*α*γₖₑ)/sqrt(γₖ) * F/sqrt(charge)
    end

    #gamma factor for phonon absorption
    Eₐ = E + ħω₀
    γₖₐ = Eₐ*( 1 + α*Eₐ)
    γₖ = E*(1 + α*E)
    A = (2*(1 + α*E)*(1 + α*Eₐ) + α*(γₖ + γₖₐ))^2
    B = -2*α*sqrt(γₖ*γₖₐ)*(4*(1 + α*E)*(1 + α*Eₐ) + α*(γₖ + γₖₐ))
    C = 4*(1 + α*E)*(1 + α*Eₐ) + (1 + 2*α*E)*(1 + 2*α*Eₐ)
    F = (A*log(abs((sqrt(γₖ) + sqrt(γₖₐ))/(sqrt(γₖ) - sqrt(γₖₐ)))) + B)/C

    Γ_pol_abs = prefactor * n₀ * (1 + 2*α*γₖₐ)/sqrt(γₖ) * F/sqrt(charge)

    if Γ_pol_em > 0.0 && Eₑ < 0.0
        println("Γ_em is: ",Γ_pol_em)
        println("New E would be: ",Eₑ)
    end

    return [Γ_pol_em, Γ_pol_abs]
end

#zeroth-order intervalley scattering with a nonparabolic band structure
function intervalley(valley::Int, E:: Float64)
    Γ_matrix = zeros((2,3))
    prefactor, weight, N, α, Eₐ, Eₑ = zeros(6)
    #define target valley. 0:Γ 1:𝐿 2:𝑋
    for target in 1:3
        prefactor = charge^2 *π * Ξ_matrix[valley, target] * Z_matrix[valley, target]/(ρ*ω_matrix[valley, target])
        weight = (2 * mass_array[target])^1.5 /(4* π^2 * ħ^3)
        N = n₀_in[valley, target]
        α = alpha_array[target]
        if valley != target
            #absorption process
            Eₐ = E + E_matrix[valley, target] - Δ[valley, target]
            if Eₐ < 0.0
                Γ_matrix[1, target] = 0.0
            else
                Γ_matrix[1, target] = prefactor * N * weight * sqrt(Eₐ*(1 + α*Eₐ))*(1 + 2*α*Eₐ)
            end

            #emission process
            Eₑ = E - E_matrix[valley, target] - Δ[valley, target]
            if Eₑ < 0.0
                Γ_matrix[2, target] = 0.0
            else
                Γ_matrix[2, target] = prefactor * (N + 1) * weight * sqrt(Eₑ*(1 + α*Eₑ))*(1 + 2*α*Eₑ)
            end
        else
            Γ_matrix[:, 2] = [0.0, 0.0]
        end

    end

    return Γ_matrix
end

function collectrates(valley::Int, E::Float64)
    Γ = zeros(9)
    Γ[1] = acoustic(valley, E)
    Γ[2:3] = polar(valley, E)
    Γ[4:9] = intervalley(valley, E) #ordered Γ, 𝐿, 𝑋
    return Γ
end

function createscatteringtable(valley::Int, n::Int)
    scatter = zeros((9,n))
    E = 0.0
    δ = .001
    for i in 1:n
        E = δ*(i-1)
        scatter[:,i] = collectrates(valley, E)
    end

    maxG = sum(scatter[:,n])
    println(scatter[:,n])
    #normalize table
    normed_table = cumsum(scatter, dims = 1) ./ maxG

    return normed_table, maxG, scatter
end

function choosemech(E::Float64, table::Array{Float64})
    if isnan(E)
        E = 0.0
    end
    index = trunc(Int, E/.001) + 1
    a = table[:,index]
    r = rand()
    mech = findfirst(x -> x > r, a)
    if mech == nothing
        mech = 10
    end
    #mech += 1

    return mech
end

function displaytable(E::Float64, table::Array{Float64})
    index = trunc(Int, E/.001) + 1
    a = table[:,index]
    println(a)
    println(table[:,index-1])
    println(table[:,index+1])
    println(index)
    println(E/.001)
    println(trunc(Int, E/.001))
end

function returnenergy(E::Float64, mech::Int, valley::Int)

    if mech == 1 || mech == 10
        return E
    elseif mech == 2
        return E - ħω₀
    elseif mech == 3
        return E + ħω₀
    elseif mech == 4
        return E + E_matrix[valley,1] - Δ[valley, 1]
    elseif mech == 5
        return E - E_matrix[valley,1] - Δ[valley, 1]
    elseif mech == 6
        return E + E_matrix[valley,2] - Δ[valley, 2]
    elseif mech == 7
        return E - E_matrix[valley,2] - Δ[valley, 2]
    elseif mech == 8
        return E + E_matrix[valley,3] - Δ[valley, 3]
    elseif mech == 9
        return E - E_matrix[valley,3] - Δ[valley, 3]
    end

    #println("mech is: ",mech)
end

end
