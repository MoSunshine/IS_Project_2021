using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
#help functions
function calc_λs(λ_up,λ_down,λ_dt,Δ_t)
    λ_p = 0
    λ_n = 0
    if Δ_t < 0
        λ_p = λ_dt
        λ_n = max(λ_dt,λ_up)
        return λ_p,λ_n
    elseif Δ_t > 0
        λ_p = min(λ_down,λ_dt)
        println(λ_p)
        λ_n = λ_dt
        println(λ_n)
        return λ_p,λ_n
    else
        λ_p = λ_dt
        λ_n = λ_dt
        return λ_p,λ_n
    end
end
calc_λs(1,2,3,2)