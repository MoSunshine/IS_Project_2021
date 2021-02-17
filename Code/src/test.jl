using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX
##Scenario generation
β = 0.1
α = 0.95
##Read data
scenarios = DataFrame(XLSX.readtable("C:\\Users\\wmd852\\Documents\\Doktorantenkurse\\Ketter\\Code_Gruppenabgabe\\Code\\data\\Dummy Scenarios.xlsx", "Results")...)
Ω = size(scenarios, 1);
P_Max = 80
P_S_Max = 40 
P_W_Max = 40
d_t = 1 
λ_D = scenarios[["lambda_D"]] 