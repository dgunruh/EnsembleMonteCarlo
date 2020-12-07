#include("Carriers.jl")
include( "Scattering.jl")
include("Carriers.jl")
include("Constants.jl")
using PyPlot
pygui(true)
using .Scattering
using .Carriers
using .Constants

#assign constants
δt = 2e-16
t_max = 10e-12
t_steps = trunc(Int, t_max/δt)
Fx = [1e5, 2e5, 5e5, 1e6, 2e6, 3e6, 5e6, 1e7, 2e7]

#intialize scattering tables
Γ_table, Γ_max, gS = createscatteringtable(1, 23000)
L_table, L_max, lS = createscatteringtable(2, 23000)
X_table, X_max, xS = createscatteringtable(3, 23000)

carriers = initializecarriers([2e5, 0.0, 0.0], [Γ_max, L_max, X_max])
dτ = [flight(carrier) for carrier in carriers]
t_evolve = 0.0
mech = 1
newenergy = 0.0

# for c in carriers
#     println(c.k)
# end
println("done intializing")
#tables for average velocities, energies, and numbers
averagev_table = [[],[],[]]
average_e_table = [[],[],[]]
population_table = [[],[],[]]
for valley in 1:3
    push!(averagev_table[valley], averagevelocity(carriers, valley))
    push!(average_e_table[valley], averageenergy(carriers, valley))
    push!(population_table[valley], population(carriers, valley))
end

# for c in carriers
#     println(c.k)
# end
println("done averaging")

#perform evolution in time
for t in 1:t_steps
    for i in 1:length(carriers)
        if dτ[i] >= δt
            t_evolve = δt
            drift!(carriers[i], t_evolve)
            dτ[i] -= δt
        else
            t_evolve = dτ[i]
            drift!(carriers[i], t_evolve)

            while dτ[i] < δt
                if carriers[i].valley_indice == 1
                    mech = choosemech(carriers[i].energy, Γ_table)
                    newenergy = returnenergy(carriers[i].energy, mech, carriers[i].valley_indice)
                elseif carriers[i].valley_indice == 2
                    mech = choosemech(carriers[i].energy, L_table)
                    newenergy = returnenergy(carriers[i].energy, mech, carriers[i].valley_indice)
                else
                    mech = choosemech(carriers[i].energy, X_table)
                    newenergy = returnenergy(carriers[i].energy, mech, carriers[i].valley_indice)
                end

                try
                    scatterk!(carriers[i], mech, newenergy)
                catch
                    if carriers[i].valley_indice == 1
                        displaytable(carriers[i].energy, Γ_table)
                        displaytable(carriers[i].energy, gS)
                    elseif carriers[i].valley_indice == 2
                        displaytable(carriers[i].energy, L_table)
                        displaytable(carriers[i].energy, lS)
                    else
                        displaytable(carriers[i].energy, X_table)
                        displaytable(carriers[i].energy, xS)
                    end
                end
                #scatterk!(carriers[i], mech, newenergy)
                t_evolve = flight(carriers[i])
                #println("evolve time: ",t_evolve)
                δt_prime = δt - dτ[i]

                if t_evolve < δt_prime
                    drift!(carriers[i], t_evolve)
                else
                    drift!(carriers[i], δt_prime)
                end

                dτ[i] += t_evolve
            end
            dτ[i] -= δt
        end
        #print("Next carrier")
    end

    for valley in 1:3
        push!(averagev_table[valley], averagevelocity(carriers, valley))
        push!(average_e_table[valley], averageenergy(carriers, valley))
        push!(population_table[valley], population(carriers, valley))
    end

    if t%1000 == 0
        println("Time step is: ",t)
    end
    #print("time step is: ",t)
end

time = range(0, stop = t_max, step = δt)
plot(time, average_e_table[1], color="red", linewidth=2.0, linestyle="--", label = "Γ valley")
plot(time,average_e_table[2], color="blue", linewidth=2.0, linestyle="-", label = "L valley")
plot(time, average_e_table[3], color="green", linewidth=2.0, linestyle=":", label = "X valley")
title("Average energy vs time")
legend()
