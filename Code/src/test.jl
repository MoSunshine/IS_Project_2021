using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
##Function to optimize th hybrid model 
##@param α - Alpha value for the model
##@param β - Beta value for the model
##@param scenario_path - Path to the Excle File with the scenarios
##@param P_S_Max - Production limit of solar power
##@param P_W_Max - Production limit of wind power
##@return Revenue of optimized model
function calc_optimazation_adj(α,β,scenario_path,P_S_Max,P_W_Max)
    ##Read data
    scenarios = DataFrame(XLSX.readtable(scenario_path, "Results")...)
    ##Definition and declartion of variables
    Ω = size(scenarios, 1);
    P_Max = P_S_Max + P_W_Max
    d_t_v = 1 
    T = 1
    N_T_one = 1
    λ_D = scenarios[!,"lambda_D"]
    λ_A = scenarios[!,"lambda_A"]
    if P_S_Max != 0
        P_S = scenarios[!,"P_t_S"]
        P_S_τ = scenarios[!,"P_tau_S"]
    else 
        P_S = scenarios[!,"Zeros"]
        P_S_τ = scenarios[!,"Zeros"]
    end
    if P_W_Max != 0
        P_W = scenarios[!,"P_t_W"]
        P_w_τ = scenarios[!,"P_tau_W"]
    else 
        P_W = scenarios[!,"Zeros"]
        P_w_τ = scenarios[!,"Zeros"]
    end
    r_Plus = scenarios[!,"r_plus"]
    r_Minus = scenarios[!,"r_minus"]
    O_D = scenarios[!,"O_D"]
    pi = 1/Ω * ones(Ω);
    ##Create optimazation model
    opt = with_optimizer(CPLEX.Optimizer)
    pool = Model(opt)
    set_silent(pool)
    ##Add variables to model
    @variable(pool,P_D[1:Ω,1:T]>=0)
    @variable(pool,P_A[1:Ω,1:T]>=0)
    @variable(pool,Δ_Plus[1:Ω,1:T]>=0)
    @variable(pool,Δ_Minus[1:Ω,1:T]>=0)
    @variable(pool,η[1:Ω]>=0)
    @variable(pool,ζ)
    @variable(pool,P[1:Ω,1:T]>=0)
    @variable(pool,Pτ[1:Ω,1:T]>=0)
    @variable(pool,Δ[1:Ω,1:T]>=0)
    @variable(pool,P_O[1:Ω,1:T]>=0)
    ##Add expressions
    #(1)
    @expression(pool, R[ω in 1:Ω],    
            sum(
            λ_D[ω,t]*P_D[ω,t]*d_t_v+
            λ_A[ω,t]*P_A[ω,t]*d_t_v +
            λ_D[ω,t]*Δ_Plus[ω,t]*r_Plus[ω,t]-
            λ_D[ω,t]*Δ_Minus[ω,t]*r_Minus[ω,t]
            for t in 1:T)
    )
    @expression(pool, CVaR, (ζ - 1/(1-α) * sum(pi[ω]*η[ω] for ω in 1:Ω)));
    #simplified
    @objective(pool, Max, sum(R[ω] for ω in 1:Ω) + β * CVaR);
    ##Add constrains to model
    #(2)
    @constraint(pool, Max_P_D[t in 1:T, ω in 1:Ω], 
        P_D[ω,t] <= P_Max);
    #(3)
    @constraint(pool,Value_P_O[t in 1:T, ω in 1:Ω],P_O[ω,t] == P_D[ω,t]+P_A[ω,t])
    #(4)
    @constraint(pool, Rest_Max_P_O[t in 1:T, ω in 1:Ω], 
    0<=P_O[ω,t] <= P_Max);
    #(5)
    @constraint(pool,Value_P_t_w[t in 1:T, ω in 1:Ω],P[ω,t]==P_S[ω,t]+P_W[ω,t])
    #(5.1)
    @constraint(pool,Value_P_τ_w[t in 1:T, ω in 1:Ω],Pτ[ω,t]==P_S_τ[ω,t]+P_w_τ[ω,t])
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
    #(13)
    for (t,ω,ωω) in product(1:T,1:Ω,1:Ω)
        status = false
        for (t_one) in (1:T) 
            if λ_D[ωω,t_one]!=λ_D[ω,t_one]
                status = true
                break
            end
        end
       
        for t_two in (1:N_T_one)
                if Pτ[ωω,t_two] != Pτ[ω,t_two]
                    status = true
                    break
                end
            end
        if status == false
            @constraint(pool,P_A[ω,t] == P_A[ωω,t])
        end
    end
    #(14)
    @constraint(pool,Not_neg[ω in 1:Ω],-sum(
        λ_D[ω,t]*P_D[ω,t]*d_t_v+
        λ_A[ω,t]*P_A[ω,t]*d_t_v +
        λ_D[ω,t]*(Δ_Plus[ω,t]+r_Plus[ω,t]-Δ_Minus[ω,t]+r_Minus[ω,t])
        for t in 1:T)+ζ-η[ω]<=0)
    #(15)
    @constraint(pool,Rest_η[ω in 1:Ω],η[ω]>=0)
    ##Optimazie model
    optimize!(pool)
    termination_status(pool)
    return  objective_value(pool)
end
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
    λ_D = scenarios[!,"lambda_D"]
    if P_S_Max != 0
        P_S = scenarios[!,"P_t_S"]
    else 
        P_S = scenarios[!,"Zeros"]
    end
    if P_W_Max != 0
        P_W = scenarios[!,"P_t_W"]
    else 
        P_W = scenarios[!,"Zeros"]
    end
    r_Plus = scenarios[!,"r_plus"]
    r_Minus = scenarios[!,"r_minus"]
    O_D = scenarios[!,"O_D"]
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
            λ_D[ω,t]*Δ_Plus[ω,t]*r_Plus[ω,t]-
            λ_D[ω,t]*Δ_Minus[ω,t]*r_Minus[ω,t]
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

#println("No adjustment")
#println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",40,40))
#println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",0,40))
#println(calc_optimazation(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",40,0))
#println("Adjustment")
#println(calc_optimazation_adj(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",100,100))
#println(calc_optimazation_adj(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",0,100))
#println(calc_optimazation_adj(0.95,0.1,"C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx",100,0))
p="C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx"
α = 0.95
list=[]
beta = []
for β in collect(0.0:0.05:1) 
    rev = calc_optimazation(α,β,p,40,40)-calc_optimazation(α,β,p,0,40)-calc_optimazation(α,β,p,40,0)
    push!(list,rev)
    push!(beta,β)
end
list2=[]
beta2 = []
for β in collect(0.0:0.05:1) 
    rev = calc_optimazation_adj(α,β,p,100,100)-calc_optimazation_adj(α,β,p,0,100)-calc_optimazation_adj(α,β,p,100,0)
    push!(list2,rev)
    push!(beta2,β)
end
println("No Adjustment")
plot(beta, list, title = "Risk aversion Revenue", ylabel= "Hybrid Revenue vs. Solo Sun/Wind", xlabel= "beta")
println("Adjustment")
plot(beta2, list2, title = "Risk aversion Revenue", ylabel= "Hybrid Revenue vs. Solo Sun/Wind", xlabel= "beta")