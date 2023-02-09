** Set directory
cd "C:/Users/giaco/Dropbox/NYU/TA Work/Quant II Spring 2023/Lab material/lab3"

***** Bootstrap regression *****

** Import the data
import delimited "mtcars.csv", clear

** Regression
reg mpg qsec

** Regression with bootstrap SEs
set seed 000
reg mpg qsec, vce(bootstrap, reps(1000))


***** WCB *****
import delimited "min_wage_teen_hours2.csv", clear

** Regression first
reg hours2 mw black female i.educ i.age i.year i.state,  robust

** Run as post-estimation command
boottest mw, cluster(styear) noci reps(999) seed(123)
