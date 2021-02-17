using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
##Scenario generation
β = 0.1
α = 0.95
##Read data
scenarios = DataFrame(XLSX.readtable("C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx", "Results")...)
scenarios
##set global variables @TODO Wert muss noch gesetzt werden
Ω = size(nrow(scenarios), 1)
T = 1 
P_Max = 80
P_S_Max = 40 
P_W_Max = 40
d_t = 1 
λ_D = [[]] ##dummy
λ_A = [[]] ##dummy

N_T_one = 20 ##dummy
##build model
opt = with_optimizer(CPLEX.Optimizer)
pool = Model(opt)
#set_silent(pool)
##Add variables to model
@variable(pool,P_D[1:T,1:Ω]>=0)
@variable(pool,P_A[1:T,1:Ω]>= 0)
@variable(pool,Δ_plus[1:T,1:Ω]>=0)
@variable(pool,Δ_minus[1:T,1:Ω]>=0)
@variable(pool,η[1:Ω]>=0)
@variable(pool,ζ)
##Add expressions
#(1)
@expression(pool, R[ω in 1:Ω],    
        sum(
        λ_D[t,ω]*P_D[t,ω]*d_t+
        λ_A[t,ω]*P_A[t,ω]*d_t+
        λ_D[t,ω]*Δ_Plus[t,ω]+r_Plus[t,ω]-
        λ_D[t,ω]*Δ_Minus[t,ω]+r_Minus[t,ω]
        for t in 1:T)
)
@expression(pool, CVaR, (ζ - 1/(1-α) * sum(π[ω]*η[ω] for ω in 1:Ω)));
#simplified
@objective(pool, Max, sum(R[ω] for ω in 1:Ω) + β * CVaR);
##Add constrains to model
#(2)
@constraint(pool, Max_P_D[t in 1:T, ω in 1:Ω], 
    P_D[t,ω] <= P_Max);
#(3)
@constraint(pool,Value_P_O[t in 1:T, ω in 1:Ω],P_O[t,ω] == P_D[t,ω]+P_A[t,ω])
#(4)
@constraint(pool, Rest_Max_P_O[t in 1:T, ω in 1:Ω], 
    0<=P_O[t,ω] <= P_Max);
#(5)
@constraint(pool,Value_P_t_w[t in 1:T, ω in 1:Ω],P[t,ω]==P_S[t,ω]+P_W[t,ω])
#(6)
@constraint(pool,Value_P_max,P_Max == P_S_Max + P_W_Max)
#(7)
@constraint(pool,Value_Δ_one[t in 1:T, ω in 1:Ω],Δ[t,ω]==d_t*(P[t,ω]-P_O[t,ω]))
#(8)
@constraint(pool,Value_Δ_two[t in 1:T, ω in 1:Ω],Δ[t,ω]== Δ_Plus[t,ω]-Δ_Minus[t,ω])
#(9)
@constraint(pool,Rest_Δ_Plus[t in 1:T, ω in 1:Ω], Δ_Plus[t,ω]<= P[t,ω]*d_t) 
#(10)
@constraint(pool,Rest_Δ_Minus[t in 1:T, ω in 1:Ω],Δ_Minus[t,ω]<= P_Max[t,ω]*d_t)
#(11) @TODO hinzufügen
#NDDA = Dict()
for (t,ω,ωω) in product(1:T, 1:Ω, 1:Ω)
    if (O_D[ω,t]+1 == O_D[ωω,t])
        @constraint(pool, P_D[t,ω] - P_D[t,ωω] <= 0)
    end
end
#(12) @TODO hier nochmal checken
for (t,ω,ωω) in product(1:T,1:Ω,1:Ω)
    if(λ_D[t,ω] == λ_D[t,ωω]){
        @constraint(pool,P_D[t,ω]==P_D[t,ωω])
    }
    end
end
#(13) @TODO hinzufügen
for (t,ω,ωω) in product(1:T,1:Ω,1:Ω)
    status = false
    for (t_one) in (1:T) 
        if λ_D[t_one,ωω]!=λ_D[t_one,ω]
            status = true
            break
        end
    end
   
    for t_two in (1:N_T_one)
            if P[t_two,ωω] != P_[t_two,ω]
                status = true
                break
            end
        end
    if status == false
        @constraint(pool,P_A[t,ω] == P_A[t,ωω])
    end
end
#(14)
@constraint(pool,Not_neg[ω in 1:Ω],-sum(
    λ_D[t,ω]*P_D[t,ω]*d_t+
    λ_A[t,ω]*P_A[t,ω]*d_t+
    λ_D[t,ω]*(Δ_Plus[t,ω]+r_Plus[t,ω]-Δ_Minus[t,ω]+r_Minus[t,ω])
    for t in 1:T)+ζ-η[ω]<=0)
#(15)
@constraint(pool,Rest_η[ω in 1:Ω],η[ω]>0)

##optimization
optimize!(pool)
termination_status(pool)
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