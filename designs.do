* Self-contained Stata exercises; all observations below are teaching data.
version 16
clear all
set more off

* OLS, FWL, and conditional comparisons.
input byte(x d) double y
0 0 10
0 0 10
0 1 14
1 0 20
1 1 24
1 1 24
end
regress y d
assert abs(_b[d]-(62/3-40/3)) < 1e-8
regress y d i.x
scalar adjusted = _b[d]
assert abs(adjusted-4) < 1e-8
quietly regress d i.x
predict double rd, residuals
quietly regress y i.x
predict double ry, residuals
regress ry rd, noconstant
assert abs(_b[rd]-adjusted) < 1e-8

* Cluster-robust versus conventional variance needs many credible clusters
* for inference. The two-state table here checks the point estimate only.
clear
input byte(group post) double y
0 0 65
0 1 73
1 0 70
1 1 85
end
regress y i.group##i.post
assert abs(_b[1.group#1.post]-7) < 1e-8

* Local linear RD on a constructed piecewise linear conditional mean.
clear
set obs 40
gen double score = (_n-20.5)/10
gen byte eligible = score>=0
gen double y = 10+score+4*eligible+0.5*score*eligible
regress y i.eligible##c.score if abs(score)<1
assert abs(_b[1.eligible]-4) < 1e-8

* Assignment, take-up, reduced form and the Wald/2SLS identity.
clear
input byte(z d) double y
0 0 10
0 0 10
0 0 10
0 1 14
1 0 10
1 1 14
1 1 14
1 1 14
end
quietly regress y z
scalar itt = _b[z]
quietly regress d z
scalar first_stage = _b[z]
ivregress 2sls y (d=z)
assert abs(_b[d]-itt/first_stage)<1e-8
display "ITT=" itt ", first stage=" first_stage ", Wald=" itt/first_stage
