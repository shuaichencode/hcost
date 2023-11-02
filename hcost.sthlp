{smcl}
{* *! version 1.3.1  Jul 2019}{...}
{cmd:help hcost}{right: ({browse "http://www.stata-journal.com/article.html?article=st0399":SJ15-3: st0399})}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{hi: hcost} {hline 2 }}Estimation of mean health care costs within a time horizon with possibly censored data{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 15 2}
{cmdab:hcost}
{it:idvar}
{it:costvar}
{ifin}{cmd:,}
{opth start(varname)}
{opth stop(varname)}
{opt l(#)}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opth start(varname)}}specify variable that contains the starting time of cost{p_end}
{p2coldent:* {opth stop(varname)}}specify variable that contains the end time of cost{p_end}
{p2coldent:* {opt l(#)}}specify time limit for calculating mean cost and mean survival time{p_end}
{synopt :{opth group(varname)}}specify variable that contains the treatment group information if more than one{p_end}
{synopt :{opt drsurv(#)}}specify annual discount rate for survival time; default is {cmd:drsurv(0)}{p_end}
{synopt :{opt drcost(#)}}specify annual discount rate for cost; default is {cmd:drcost(0)}{p_end}
{synopt :{opt icer(#)}}specify option for calculation of the incremental cost-effectiveness ratio (ICER); {cmd:1} for calculating ICER and confidence interval, {cmd:0} for not; default is {cmd:icer(0)}{p_end}
{synopt :{opt k(#)}}specify number of days in a time interval for discounting costs; default is {cmd:k(90)}{p_end}
{synopt :{opt level(#)}}specify confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt method(#)}}specify method of mean cost calculation; {cmd:1} for using cost history, {cmd:0} for only using total cost; default is {cmd:method(1)}{p_end}
{synoptline}{...}
{p2colreset}
{p 4 6 2}
* {opt start()} is required unless {cmd:method()} is {cmd:0} and {cmd:drcost()} is {cmd:0}.{p_end}
{p 4 6 2}
* {opt stop()} is required unless {cmd:method()} is {cmd:0} and {cmd:drcost()} is {cmd:0}.{p_end}
{p 4 6 2}
* {opt l()} is required.{p_end}
{p 4 6 2}
You must {cmd:stset} your data before using {cmd:hcost}; see {manhelp stset ST}. {p_end}


{title:Description}

{pstd}
The {opt hcost} command estimates mean cost with possibly censored data.  The
methods used are given in the {it:{help hcost##references:Reference}} section.
Your data are assumed to be {cmd:stset}; see {manhelp stset ST}.


{title:Options}

{phang}
{opth start(varname)} specifies the variable that contains the start time of
accumulated costs (in days).  {cmd:start()} is required unless {cmd:method()}
is {cmd:0} and {cmd:drcost()} is {cmd:0}.

{phang}
{opth stop(varname)} specifies the variable that contains the end time of
accumulated costs (in days), which needs to be within the survival time. {cmd:stop()} is required unless {cmd:method()}
is {cmd:0} and {cmd:drcost()} is {cmd:0} 

{phang}
{opt l(#)} specifies the time limit L (in days) for calculating mean cost and
mean survival time.  {cmd:l()} is required.

{phang}
{opth group(varname)} specifies the variable that contains the treatment group
information if there is more than one group.  If there is only one group, this
option need not be specified.  There cannot be more than two groups.

{phang}
{opt drsurv(#)} specifies the annual discount rate for survival time.  The
value must be between 0 and 1.  The default is {cmd:drsurv(0)}; that is,
survival time is already discounted in the data.

{phang}
{opt drcost(#)} specifies the annual discount rate for cost.  The value must
be between 0 and 1.  The default is {cmd:drcost(0)}; that is, cost is already
discounted in the data.  Each cost entry is discounted according to its start
date.  If the time interval for the cost entry is too big, it will be divided
into {cmd:k()}-day subintervals first.

{phang}
{opt icer(#)} specifies the option for performing cost-effectiveness analysis:
{cmd:1} for calculating the ICER and its confidence interval and {cmd:0}
otherwise.  The default is {cmd:icer(0)}.

{phang}
{opt k(#)} specifies the number of days in a time interval for discounting
costs.  The default is {cmd:k(90)}.  When the time interval between starting
and stopping times is smaller than {cmd:k()}, the costs accumulated are
discounted using the start date; if the time interval is bigger than
{cmd:k()}, then the interval is divided into {cmd:k()}-day subintervals first,
and the costs at each subinterval are discounted at the starting time of the
subinterval.

{phang}
{opt level(#)} specifies the level for confidence intervals.  The value must
be between 10 and 99.  The default is {cmd:level(95)}.

{phang}
{opt method(#)} specifies the method for mean costs estimation: {cmd:1} for
using cost history (Zhao and Tian [2001] estimator), and {cmd:0} for using
only total costs (Bang and Tsiatis [2000] estimator).  The default is
{cmd:method(1)}.  If {cmd:method()} is {cmd:0} and {cmd:drcost()} is {cmd:0},
{cmd:start()} and {cmd:stop()} are not required.  When {cmd:start()} and
{cmd:stop()} are not provided, the costs are assumed to be collected within
the time limit L.


{title:Remarks}

{pstd}
Selection of {cmd:l()}: Although it is allowed to use large {cmd:l()} (for example, setting l() to be the largest follow-up time), 
such large {cmd:l()} is not recommended. Large {cmd:l()} will lead to huge variance, and hence the results will be unstable and less meaningful. 
It is recommended to select {cmd:l()} when there are still enough observations being followed. 

{pstd}
L-truncated cost: If {cmd:stop()} exceeds L-truncated follow-up time, the costs will be prorated. If {cmd:stop()} is not provided, the cost will be assumed to be evenly spreading from 0 to follow-up time, and will be prorated by L-truncated follow-up time.

{title:Examples}

{pstd}
Setup{p_end}
{phang2}{cmd:. use example}{p_end}
{phang2}{cmd:. stset surv, failure(delta)}{p_end}

{pstd}
One-group estimation{p_end}
{phang2}{cmd:. hcost id cost, start(start) stop(stop) l(1461)}{p_end}

{pstd}
One-group estimation with discounting{p_end}
{phang2}{cmd:. hcost id cost, start(start) stop(stop) l(1461) drsurv(0.03) drcost(0.03)}{p_end}

{pstd}
Two-group estimation using default Zhao and Tian (2001) estimator{p_end}
{phang2}{cmd:. hcost id cost, start(start) stop(stop) l(1461) group(trt)}{p_end}

{pstd}
Two-group estimation using Bang and Tsiatis (2000) estimator{p_end}
{phang2}{cmd:. hcost id cost, start(start) stop(stop) l(1461) group(trt) method(0)}{p_end}


{title:Stored results}

{p 4 4 2}
{cmd:hcost} stores the following in {cmd:e()} (when {cmd:group()} is not
specified):

{synoptset 15 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(n)}}number of patients{p_end}
{synopt:{cmd:e(nobs)}}number of cost observations{p_end}
{synopt:{cmd:e(pcensor)}}censoring rate{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimator{p_end}

{pstd}If {cmd:group()} is specified, similar results are stored in
{cmd:e(n1)}, {cmd:e(n2)}, {cmd:e(nobs1)}, {cmd:e(nobs)2}, {cmd:e(pcensor1)},
{cmd:e(pcensor2)}, {cmd:e(b1)}, {cmd:e(b2)}, {cmd:e(V1)}, and {cmd:e(V2)} for
groups 1 and 2.{p_end}


{marker references}{...}
{title:References}

{phang}
Bang, H., and A. A. Tsiatis. 2000. Estimating medical costs with censored
data.  {it:Biometrika} 87: 329-343.

{phang}
Zhao, H., and L. Tian. 2001.  On estimating medical cost and
incremental cost-effectiveness ratios with censored data.  {it:Biometrics} 57: 1002-1008.


{title:Authors}

{pstd}Shuai Chen{p_end}
{pstd}Division of Biostatistics{p_end}
{pstd}Department of Public Health Sciences{p_end}
{pstd}University of California{p_end}
{pstd}Davis, CA{p_end}
{pstd}shschen@ucdavis.edu{p_end}

{pstd}Jennifer Rolfes{p_end}
{pstd}College of Health Solutions{p_end}
{pstd}Arizona State University{p_end}
{pstd}Phoenix, Arizona{p_end}

{pstd}Hongwei Zhao{p_end}
{pstd}Department of Epidemiology and Biostatistics{p_end}
{pstd}School of Public Health{p_end}
{pstd}Texas A&M University{p_end}
{pstd}College Station, TX{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 3: {browse "http://www.stata-journal.com/article.html?article=st0399":st0399}
{p_end}
