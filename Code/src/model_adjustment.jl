using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX, CSV, Statistics, PGFPlotsX
##Function to optimize the hybrid model 
##@param α - Alpha value for the model
##@param β - Beta value for the model
##@param scenario_path - Path to the Excle file with the scenarios
##@param P_S_Max - Production limit of solar power
##@param P_W_Max - Production limit of wind power
##@return Revenue of optimized model
function calc_optimization(α,β,scenario_path,P_S_Max,P_W_Max)
    ##Read data
    scenarios = CSV.read(scenario_path,DataFrame,delim=",")
    ##Definition and declaration of variables
    Ω = size(scenarios, 1);
    P_Max = P_S_Max + P_W_Max
    D_T_V = 1 
    T = 1
    N_T_one = T
    λ_D = scenarios[!,"lambda_D"]
    λ_A = scenarios[!,"lambda_A"]
    Pi = scenarios[!,"Pi"]
    if P_S_Max != 0
        P_S = scenarios[!,"P_t_S"]
        P_S_τ = scenarios[!,"P_tau_S"]
    else 
        P_S = scenarios[!,"Zeros"]
        P_S_τ = scenarios[!,"Zeros"]
    end
    if P_W_Max != 0
        P_W = scenarios[!,"P_t_W"]
        P_W_τ = scenarios[!,"P_tau_W"]
    else 
        P_W = scenarios[!,"Zeros"]
        P_W_τ = scenarios[!,"Zeros"]
    end
    r_Plus = scenarios[!,"r_plus"]
    r_Minus = scenarios[!,"r_minus"]
    O_D = scenarios[!,"O_Rank"]
    pi = 1/Ω * ones(Ω);
    ##Create optimazation model
    opt = with_optimizer(CPLEX.Optimizer)
    pool = Model(opt)
    set_silent(pool)
    ##Add variables to model
    @variable(pool,P_D[1:Ω,1:T]>=0)
    @variable(pool,P_A[1:Ω,1:T])
    @variable(pool,Δ_Plus[1:Ω,1:T]>=0)
    @variable(pool,Δ_Minus[1:Ω,1:T]>=0)
    @variable(pool,η[1:Ω]>=0)
    @variable(pool,ζ)
    @variable(pool,P[1:Ω,1:T]>=0)
    @variable(pool,Pτ[1:Ω,1:N_T_one]>=0)
    @variable(pool,Δ[1:Ω,1:T])
    @variable(pool,P_O[1:Ω,1:T]>=0)
    ##Add expressions
    #(1)
    @expression(pool, R[ω in 1:Ω],    
            sum(
            λ_D[ω,t]*P_D[ω,t]*D_T_V+
            λ_A[ω,t]*P_A[ω,t]*D_T_V +
            λ_D[ω,t]*Δ_Plus[ω,t]*r_Plus[ω,t]-
            λ_D[ω,t]*Δ_Minus[ω,t]*r_Minus[ω,t]
            for t in 1:T)
    )
    @expression(pool, CVaR, (ζ - 1/(1-α) * sum(Pi[ω]*η[ω] for ω in 1:Ω)));
    #simplified
    @objective(pool, Max, (1-β)*sum(Pi[ω]*R[ω] for ω in 1:Ω) + β * CVaR);
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
    @constraint(pool,Value_P_τ_w[t in 1:T, ω in 1:Ω],Pτ[ω,t]==P_S_τ[ω,t]+P_W_τ[ω,t])
    #(6)
    @constraint(pool,Value_P_max,P_Max == P_S_Max + P_W_Max)
    #(7)
    @constraint(pool,Value_Δ_one[t in 1:T, ω in 1:Ω],Δ[ω,t]==D_T_V*(P[ω,t]-P_O[ω,t]))
    #(8)
    @constraint(pool,Value_Δ_two[t in 1:T, ω in 1:Ω],Δ[ω,t]== Δ_Plus[ω,t]-Δ_Minus[ω,t])
    #(9)
    @constraint(pool,Rest_Δ_Plus[t in 1:T, ω in 1:Ω], Δ_Plus[ω,t]<= P[ω,t]*D_T_V) 
    #(10)
    @constraint(pool,Rest_Δ_Minus[t in 1:T, ω in 1:Ω],Δ_Minus[ω,t]<= P_Max*D_T_V)
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
    NDDA = Dict()
    for (t,ω) in product(1:T,1:Ω)
        for ωω in 1:Ω
            status = false
            for (t_one) in (1:T) 
                if λ_D[ωω,t_one]!=λ_D[ω,t_one]
                    status = true
                    break
                end
            end
        
            for t_two in (1:N_T_one)
                    if P_S_τ[ωω,t_two]+ P_W_τ[ωω,t_two] != P_S_τ[ω,t_two]+P_W_τ[ω,t_two]
                        status = true
                        break
                    end
            end
            if status == false
                NDDA[(t,ω,ωω)] = @constraint(pool,P_A[ω,t] == P_A[ωω,t])
            end
        end
    end
    #(14)
    @constraint(pool,Not_neg[ω in 1:Ω],-sum(
        λ_D[ω,t]*P_D[ω,t]*D_T_V+
        λ_A[ω,t]*P_A[ω,t]*D_T_V +
        λ_D[ω,t]*(Δ_Plus[ω,t]*r_Plus[ω,t]-Δ_Minus[ω,t]*r_Minus[ω,t])
        for t in 1:T)+ζ-η[ω]<=0)
    #(15)
    @constraint(pool,Rest_η[ω in 1:Ω],η[ω]>=0)
    ##Optimazie model
    optimize!(pool)
    termination_status(pool)
    #Return revenue
    return objective_value(pool)
end
##Main part
##System path to scenario csv-file and declartion of variables
path = "..\\data\\outputfile E_A_negativ.csv"
α = 0.95
list=[]
beta = []
##iteration over all β between 0.0 and 1.0 in 0.05 steps
for β in collect(0.0:0.05:1) 
    rev = calc_optimization(α,β,path,40000,40000)-calc_optimization(α,β,path,0,40000)-calc_optimization(α,β,path,40000,0)
    push!(list,rev)
    push!(beta,β)
end
println("****Finished calculation start printing results****")
##printing results
pgfplotsx()
figure = @pgf Axis(    
        {
            xlabel = "beta",
            ylabel = "Hybrid Revenue vs. Solo Sun/Wind"
        }
    ,
    Plot(
        {
            no_marks
        },
        Table(beta,list)))
pgfsave("..\\output\\figure_adjustment_model.tex",figure)
println("****Finsihed printing****")