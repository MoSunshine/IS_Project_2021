using Pkg
using JuMP, Plots, IterTools, CPLEX, DataFrames, XLSX, CSV, Statistics, PGFPlotsX

scenario_path = "..\\data\\outputfile E_A_negativ.csv"
scenarios = CSV.read(scenario_path,DataFrame,delim=",")
##Definition and declaration of variables
P_t_S = scenarios[!,"P_t_S"]
P_t_W = scenarios[!,"P_t_W"]
P_t = P_t_S + P_t_W

println("Mean S,W,Ges: ", mean(P_t_S), "    ", mean(P_t_W), "    ",mean(P_t))
println("Variance S,W,Ges: ",  var(P_t_S), "    ", var(P_t_W), "     ", var(P_t))
println("Std S,W,Ges: ",  std(P_t_S), "    ", std(P_t_W), "     ", std(P_t))
println(sum(mean(P_t_S)+mean(P_t_W)))
println(mean(P_t))
println(sum(var(P_t_S)+var(P_t_S)))
println(var(P_t))
println(cor(P_t_S,P_t_W))



