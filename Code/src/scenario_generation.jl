using Pkg
using Statistics, Distributions, Clustering, ARCHModels, Distances, Plots, JuMP, IterTools, CPLEX, DataFrames, XLSX, Random, CSV
path_input = "..\\data\\Historical Data.xlsx"
path_output = "..\\data\\outputfile.csv"
sheet_input ="Results"
Header=["P_tau" "P_tau_S" "P_tau_W" "P_t" "P_t_S" "P_t_W" "lambda_D" "lambda_A" "r_plus" "r_minus" "Pi" "Zeros" "O_Rank"]

price_scenarios =
[
50 40 1 1.5
50 53 1 1.5
20 10 1 1.5
20 23 1 1.5
50 40 0.9 1
50 53 0.9 1
20 10 0.9 1
20 23 0.9 1
]

price_probability=
[
0.4 0.1 0.7
0.4 0.9 0.7
0.6 0.1 0.7
0.6 0.9 0.7
0.4 0.1 0.3
0.4 0.9 0.3
0.6 0.1 0.3
0.6 0.9 0.3    
]

##Function to returns a vector with the n-moving averages for vs[n:length(vs)]
##@param vs - value vector
##@param n - number of elements to consider for the moving average
moving_average(vs,n) = [sum(@view vs[i:(i+n)])/n for i in 1:(length(vs)-(n))]

##Function to return the residual vector of a ARMA{1,1} model
##@param y - the data the ARMA model was fitted to
##@param μ - the mean that was used for the ARMA model
##@param ϕ - coefficient from the ARMA model
##@param θ - coefficient from the ARMA model
function calc_residuals(y, μ, ϕ, θ)
    ϵ = [μ-y[1]]
    y_hat = [μ]
    T = size(y, 1)
    for t in 2:T
        yt_hat = μ+ϕ*(y[t-1]-μ)+θ*ϵ[t-1]
        ϵt = y[t]-yt_hat
        push!(ϵ, ϵt)
        push!(y_hat, yt_hat)
    end
    return ϵ
end

##Function to return a clustered sample from a multivariate normal Distributions
##@param Dist - MvNormal Random Variable, from which the sample is supposed to be generated
##@param μ_S - mean offset for the first random variable
##@param μ_W - mean offset for the second random variable
##@param T - initial Sample size
##@param kmT - number of clusters to be created from the initial sample. for each cluster the kmedoid is returned based on Euclidean distance
function Residual_Sample(Dist, μ_S, μ_W, T, kmT)
    ReDo = true
    y = []
    while ReDo
        y1 = []
        y2 = []
        ykm_1= []
        ykm_2= []
        for t in 1:T
            yt = rand(Dist)
            push!(y1, yt[1])
            push!(y2, yt[2])
        end
        y = [y1 y2]
        dist = pairwise(Euclidean(), y, dims=1)
        kmed = kmedoids(dist,kmT)
        for i in 1:kmT
            push!(ykm_1, μ_S + y1[kmed.medoids[i]])
            push!(ykm_2, μ_W + y2[kmed.medoids[i]])
        end
        y = [ykm_1 ykm_2]
        if minimum(y)>0
            ReDo=false
        end
    end
    return y
end
##function to faktor the values of each row in a matrix to calculate the probability for each scenario - returns a single column
##@param prob_vector - the probability vector to be factored
function calc_prob(prob_vector)
    result=[]
    for i in 1:size(prob_vector,1)
        P=1.0
        for j in 1:size(prob_vector,2)
            P=P*prob_vector[i,j]
        end
        push!(result,P)
    end
    return result
end

##function to compute the entire scenario data based upon a MV residual distribution and the coefficients from the ARMA models
##@param Residual_Sample - Ptau_S/Ptau_W Values to be used as a starting value to generate Pt_S/Pt_W
##@param prices - a matrix of the price scenarios
##@param price_probability - the probability vector corresponding to the price scenarios 
##@param ϵ_Dist - MvNormal Random Variable, from which the Pt_S/Pt_W sample is supposed to be generated
##@param μ_S - the mean that was used for the ARMA model (first variable)
##@param ϕ_S - coefficient from the ARMA model
##@param θ_S - coefficient from the ARMA model
##@param μ_W - the mean that was used for the ARMA model (second variable)
##@param ϕ_W - coefficient from the ARMA model
##@param θ_W - coefficient from the ARMA model
##@param T - initial Sample size
##@param kmT - number of clusters to be created from the initial sample. for each cluster the kmedoid is returned based on Euclidean distance
function sample_PtauPt(Residual_Sample, prices, price_probability, ϵ_Dist, μ_S, ϕ_S, θ_S, μ_W, ϕ_W, θ_W, T, kmT)
    ReDo = true
    results = []
    while ReDo
        PT_S= []
        PT_W = []
        results = []

        P_tau = []
        P_tau_S = []
        P_tau_W = []
        P_t = []
        P_t_S = []
        P_t_W = []
        lambda_D = []
        lambda_A = []
        r_plus = []
        r_minus = []
        Pi =[]

        for i in 1:size(Residual_Sample, 1)
            yt_S = []
            yt_W = []
            ykm_S = []
            ykm_W = []
            for t in 1:T
                ϵt = rand(ϵ_Dist)
                #y0 = μ
                #y1 = μ + ϵt
                #y2 = μ+ϕ*(y[t-1]-μ)+ϵt+θ*ϵ[t-1]
                yt_S = μ_S+ϕ_S*(Residual_Sample[i,1] - μ_S)+ϵt[1]+θ_S*(Residual_Sample[i,1] - μ_S)
                yt_W = μ_W+ϕ_W*(Residual_Sample[i,2] - μ_W)+ϵt[2]+θ_W*(Residual_Sample[i,2] - μ_W)
                push!(PT_S, yt_S)
                push!(PT_W, yt_W)
            end
            y = [PT_S PT_W]
            dist = pairwise(Euclidean(), y, dims=1)
            kmed = kmedoids(dist,kmT)
            for j in 1:kmT
                push!(ykm_S, PT_S[kmed.medoids[j]])
                push!(ykm_W, PT_W[kmed.medoids[j]])
            end

            for a in 1:size(prices,1)
                for b in 1:kmT

                    push!(P_tau,round(Residual_Sample[i,1] + Residual_Sample[i,2],digits=2))
                    push!(P_tau_S,round(Residual_Sample[i,1],digits=2))
                    push!(P_tau_W,round(Residual_Sample[i,2],digits=2))
                    push!(P_t,round(ykm_S[b] + ykm_W[b],digits=2))
                    push!(P_t_S,round(ykm_S[b],digits=2))
                    push!(P_t_W,round(ykm_W[b],digits=2))
                    push!(lambda_D,round(price_scenarios[a,1],digits=2))
                    push!(lambda_A,round(price_scenarios[a,2],digits=2))
                    push!(r_plus,round(price_scenarios[a,3],digits=2))
                    push!(r_minus,round(price_scenarios[a,4],digits=2))
                    push!(Pi,price_probability[a]/(size(Residual_Sample, 1)*kmT))
                end
            end
        end
        results = [P_tau P_tau_S P_tau_W P_t P_t_S P_t_W lambda_D lambda_A r_plus r_minus Pi]
        if minimum(results)>0
            ReDo=false
        end
    end
    return results
end

function ordinal_rank(input)
    output=[]
    sort_unique_input = sort(unique(input))
    rank_input = ordinalrank(sort_unique_input)
    for i in 1:size(input,1)
       for j in 1:size(sort_unique_input,1)
            if  input[i] == sort_unique_input[j]
                push!(output,rank_input[j])
            end
        end
    end
    return output
end

# -------------------------------------------------- main part -------------------------------------------------

scenarios = DataFrame(XLSX.readtable(path_input, sheet_input)...)
P_t_W = scenarios[!,"Wind - 12"]
P_t_S = scenarios[!,"Solar - 12"]
μ_W = mean(P_t_W)
μ_S = mean(P_t_S)
MA_W = moving_average(P_t_W,30)
MA_S = moving_average(P_t_S,30)
P_t_W=P_t_W[31:2100]
P_t_S=P_t_S[31:2100]
P_t_W=P_t_W.-MA_W
P_t_S=P_t_S.-MA_S
P_t_W=P_t_W.+μ_W
P_t_S=P_t_S.+μ_S

μ_W = mean(P_t_W)
μ_S = mean(P_t_S)

arma = selectmodel(ARCH{0}, P_t_S .- μ_S; meanspec=ARMA{1,1}, criterion=bic)
ϕ_S, θ_S = arma.meanspec.coefs[2:3]
normal_residuals_S = calc_residuals(P_t_S, μ_S, ϕ_S, θ_S)

arma = selectmodel(ARCH{0}, P_t_W .- μ_W; meanspec=ARMA{1,1}, criterion=bic)
ϕ_W, θ_W = arma.meanspec.coefs[2:3]
normal_residuals_W = calc_residuals(P_t_W, μ_W, ϕ_W, θ_W)

COV_WS = cov(normal_residuals_W,normal_residuals_S)
Var_S = var(normal_residuals_S)
Var_W = var(normal_residuals_W)
Cov_mat = [[Var_S,COV_WS] [COV_WS,Var_W]]

D = MvNormal([mean(normal_residuals_S),mean(normal_residuals_W)],Cov_mat)

Residuals = Residual_Sample(D, μ_S, μ_W, 100, 4)

Results = sample_PtauPt(Residuals, price_scenarios, calc_prob(price_probability), D, μ_S, ϕ_S, θ_S, μ_W, ϕ_W, θ_W, 100, 4)
O_Rank = ordinal_rank(Results[:,7])
Z = zeros(size(Results,1),1)
Results=[Results Z O_Rank]

df = convert(DataFrame, [Header;Results])
CSV.write(path_output,df)
print("DONE")