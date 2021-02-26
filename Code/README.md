AUTHORS:
=======
David Könen, Daniel Piersig, Laura Poreschack, Moritz Wegener

REQUIREMENTS:
=======
Code was tested and is runable using the following software versions
Microsoft Windows 10 Education - Version 10.0.18363 Build 10.0.18363
* julia version 1.5.3
* ARCHModels v1.3.0
* Clustering v0.14.2
* CPLEX v0.7.6
* CSV v0.8.3
* DataFrames v0.22.5
* Distances v0.10.2
* Distributions v0.24.14
* IterTools v1.3.0
* JuMP v0.21.6
* PGFPlotsX v1.2.1
* Plots v1.10.3
* XLSX v0.7.6


RUN
=======

* run the file model_adjustment.jl to calculate the model with the adjustment market and print out the corresponding figure
* run the file model_no_adjustment.jl to calculate the model without the adjustment market and print out the corresponding figure
* run the file variancecheck.jl to calculate and print out the basic statistical values
* run the file scenario generation.jl to generate scenarios for the model
* if you want to reproduce the results from the paper use the scenarios ./data/Example_Data_Sets/outputfile_used.csv
