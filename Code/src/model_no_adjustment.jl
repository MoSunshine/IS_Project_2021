using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
##Function to optimize th hybrid model 
##@param α - Alpha value for the model
##@param β - Beta value for the model
##@param scenario_path - Path to the Excle File with the scenarios
##@param P_S_Max - Production limit of solar power
##@param P_W_Max - Production limit of wind power
##@return Revenue of optimized model
function calc_optimazation(α,β,scenario_path,P_S_Max,P_W_Max)
    ##Read data
    scenarios = DataFrame(XLSX.readtable(scenario_path, "Results")...)
    ##Definition and declartion of variables
    Ω = size(scenarios, 1);
    P_Max = P_S_Max + P_W_Max
    d_t_v = 1 
    T = 1
    λ_D = scenarios[["lambda_D"]]
    if P_S_Max != 0
        P_S = scenarios[["P_tau_S"]]
    else 
        P_S = scenarios[["Zeros"]]
    end
    if P_W_Max != 0
        P_W = scenarios[["P_tau_W"]]
    else 
        P_W = scenarios[["Zeros"]]
    end
    r_Plus = scenarios[["r_plus"]]
    r_Minus = scenarios[["r_minus"]]
    O_D = scenarios[["O_D"]]
    pi = 1/Ω * ones(Ω);
    ##Create optimazation model
    opt = with_optimizer(CPLEX.Optimizer)
    pool = Model(opt)
    set_silent(pool)
    ##Add variables to model
    @variable(pool,P_D[1:Ω,1:T]>=0)
    @variable(pool,Δ_Plus[1:Ω,1:T]>=0)
    @variable(pool,Δ_Minus[1:Ω,1:T]>=0)
    @variable(pool,η[1:Ω]>=0)
    @variable(pool,ζ)
    @variable(pool,P[1:Ω,1:T]>=0)
    @variable(pool,Δ[1:Ω,1:T]>=0)
    ##Add expressions
    #(1)
    @expression(pool, R[ω in 1:Ω],    
            sum(
            λ_D[ω,t]*P_D[ω,t]*d_t_v+
            λ_D[ω,t]*Δ_Plus[ω,t]+r_Plus[ω,t]-
            λ_D[ω,t]*Δ_Minus[ω,t]+r_Minus[ω,t]
            for t in 1:T)
    )
    @expression(pool, CVaR, (ζ - 1/(1-α) * sum(pi[ω]*η[ω] for ω in 1:Ω)));
    #simplified
    @objective(pool, Max, sum(R[ω] for ω in 1:Ω) + β * CVaR);
    ##Add constrains to model
    #(2)
    @constraint(pool, Max_P_D[t in 1:T, ω in 1:Ω], 
        P_D[ω,t] <= P_Max);
    ##(3)+ (4) Not not necessary
    #(5)
    @constraint(pool,Value_P_t_w[t in 1:T, ω in 1:Ω],P[ω,t]==P_S[ω,t]+P_W[ω,t])
    #(6)
    @constraint(pool,Value_P_max,P_Max == P_S_Max + P_W_Max)
    #(7)
    @constraint(pool,Value_Δ_one[t in 1:T, ω in 1:Ω],Δ[ω,t]==d_t_v*(P[ω,t]-P_D[ω,t]))
    #(8)
    @constraint(pool,Value_Δ_two[t in 1:T, ω in 1:Ω],Δ[ω,t]== Δ_Plus[ω,t]-Δ_Minus[ω,t])
    #(9)
    @constraint(pool,Rest_Δ_Plus[t in 1:T, ω in 1:Ω], Δ_Plus[ω,t]<= P[ω,t]*d_t_v) 
    #(10)
    @constraint(pool,Rest_Δ_Minus[t in 1:T, ω in 1:Ω],Δ_Minus[ω,t]<= P_Max*d_t_v)
    #(11) 
    for (t,ω,ωω) in product(1:T, 1:Ω, 1:Ω)
        if (O_D[ω,t]+1 == O_D[ωω,t])
            @constraint(pool, P_D[ω,t] - P_D[ωω,t] <= 0)
        end
    end
    #(12)
    for (t,ω,ωω) in product(1:T,1:Ω,1:Ω)
        if λ_D[ω,t] == λ_D[ωω,t]
            @constraint(pool,P_D[ω,t] == P_D[ωω,t])
        end
    end
    #(13) Not not necessary
    #(14)
    @constraint(pool,Not_neg[ω in 1:Ω],-sum(
        λ_D[ω,t]*P_D[ω,t]*d_t_v+
        λ_D[ω,t]*(Δ_Plus[ω,t]+r_Plus[ω,t]-Δ_Minus[ω,t]+r_Minus[ω,t])
        for t in 1:T)+ζ-η[ω]<=0)
    #(15)
    @constraint(pool,Rest_η[ω in 1:Ω],η[ω]>=0)
    ##Optimazie model
    optimize!(pool)
    termination_status(pool)
    return  objective_value(pool)
end
println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",40,40))
println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",0,40))
println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",40,0))