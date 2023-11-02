/*
An example to show how to use hcost function.

Each row represnets a cost entry of a person in the input data example.dta. To see how to generate the data 
from separate survival and cost data files, please see CombineData.do as an axample. 

*/

* copy hcost.ado to "C:\ado\personal"

adopath ++ .

cd "Your Path"
use example, clear

* Declare the data as survival data by stset before using hcost
stset surv, failure(delta)

* Estimate the mean survival time and the mean costs for the whole population using the default ZT estimator without discounting
hcost id cost, start(start) stop(stop) l(1461) 

* Example of improper L
hcost id cost, start(start) stop(stop) l(3000)

* Estimate the mean discounted survival time and the mean discounted costs for the whole population using the default ZT estimator
hcost id cost, start(start) stop(stop) l(1461) drsurv(0.03) drcost(0.03)

* Estimate the mean discounted survival time and the mean discounted costs for the two treatments using the BT estimator
hcost id cost, start(start) stop(stop) l(1461) group(trt) method(0) drsurv(0.03) drcost(0.03) 

* Estimate the mean discounted survival time and the mean discounted costs for the two treatments using the default ZT estimator and the option icer=1 for calculating the iCER and its CI.
hcost id cost, start(start) stop(stop) l(1461) group(trt) drsurv(0.03) drcost(0.03) icer(1)

* Saved results
ereturn list

* Mean costs and mean survival time for Group 1
matrix list e(b1)

* Covariance matrix between the mean costs and mean survival time for Group 1
matrix list e(V1)
