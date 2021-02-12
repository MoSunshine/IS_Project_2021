using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
##Scenario generation
#@@TODo Scenraio generation hier einfügen
β = 0.1
α = 0.95
##set global variables @TODO Wert muss noch gesetzt werden
Ω = size(5, 1); 
T = 3
P_Max = 0;
P_S_Max = 0
d_t = 0
λ_D = [[]]
λ_A = [[]]
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
@expression(pool, CVaR, (ζ - 1/(1-α) * sum(pi[ω]*η[ω] for ω in 1:Ω)));
#simplified
@objective(pool, Max, sum(R[ω] for ω in 1:Ω) + β * CVaR);
##Add constrains to model
#(2)
@constraint(pool, Max_P_D[t in 1:T, ω in 1:Ω], 
    P_D[t,ω] <= P_Max);
#(3)
@constraint(pool,Value_P_O[t in 1:T, ω in 1:Ω],P_O[t,ω] = P_D[t,ω]+P_A[t,ω])
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
#@constraint(pool,)
#(12) @TODO hier nochmal checken
for(t,ω,ωω) in product(1:T,1:Ω,1:Ω)
    if(λ_D[t,ω]+1 == λ_D[t,ωω]){
        @constraint(pool,P_D[t,ω]==P_D[t,ωω])
    }
    end
end
@constraint(pool,)
#(13) @TODO hinzufügen
#@constraint(pool,)
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
calc_λs(1,2,3,2)
