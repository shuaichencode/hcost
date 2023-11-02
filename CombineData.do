/*
An example to show how to combine the separate survival and cost data into the reqired data format.

Input format:
1. surv.dat: id trt delta surv
Each row represents a person.
2. costs.dat: id start stop cost
Each row represnets a cost entry of a person. 

Output format: id start stop cost trt delta surv
Each row represnets a cost entry of a person.   
*/

cd "Your Path"

*read survival data of group 0 
infile id trt delta surv using "surv0.dat",clear
save surv_d0,replace
*read cost data of group 0 
infile id start stop cost using "costs0.dat",clear
save cost0,replace
*merge survival and cost data
merge m:1 id using surv_d0
keep if _merge == 3
save costtotal0,replace

*read survival data of group 1, make sure id is unique for each person after combining two groups
infile id trt delta surv using "surv1.dat",clear
save surv_d1,replace
*read cost data of group 1 
infile id start stop cost using "costs1.dat",clear
save cost1,replace
*merge survival and cost data
merge m:1 id using surv_d1
keep if _merge == 3
save costtotal1,replace

*combine data of group 0 and group 1
append using costtotal0
describe
sort trt id start stop
drop _merge

save example, replace


