### R code from vignette source 'BayesRGMM.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt = "R> ", digits = 4, show.signif.stars = FALSE)


###################################################
### code chunk number 2: Simulation 1： MCD Correlation Structure
###################################################
# Simulation study for MCD correlation structure
library(BayesRGMM)
rm(list=ls(all=TRUE))

Fixed.Effs = c(-0.2, -0.3, 0.8, -0.4) 
P = length(Fixed.Effs)
q = 1
T = 5
N = 100
num.of.iter = 100

HSD.para = c(-0.5,  -0.3)
a = length(HSD.para)
w = array(runif(T*T*a), c(T, T, a))

for(time.diff in 1:a)
	w[, , time.diff]=1*(as.matrix(dist(1:T, 1:T, method="manhattan"))
	                    ==time.diff)

HSD.sim.data = SimulatedDataGenerator(Num.of.Obs = N, Num.of.TimePoints = T, 
      Fixed.Effs = Fixed.Effs, Random.Effs = list(Sigma = 0.5*diag(1), df=3), 
      Cor.in.DesignMat = 0., Missing = list(Missing.Mechanism = 2, 
      RegCoefs = c(-1.5, 1.2)), Cor.Str = "HSD", 
      HSD.DesignMat.para = list(HSD.para = HSD.para, DesignMat = w))

hyper.params = list(
        sigma2.beta = 1,
        sigma2.delta = 1,
        v.gamma = 5,
        InvWishart.df = 5,
        InvWishart.Lambda = diag(q) )

HSD.output = BayesRobustProbit(fixed = as.formula(paste("y~-1+", 
             paste0("x", 1:P, collapse="+"))), data=HSD.sim.data$sim.data, 
             random = ~ 1, HS.model = ~IndTime1+IndTime2, subset = NULL, 
             na.action='na.exclude', hyper.params = hyper.params, 
             num.of.iter = num.of.iter)

original = options(digits = 4)
Model.Estimation = BayesRobustProbitSummary(HSD.output)

cat("\nCoefficients:\n")
print(Model.Estimation$beta.est.CI)

cat("\nParameters in HSD model:\n")
print(Model.Estimation$delta.est.CI)

cat("\nRandom effect: \n")
print(Model.Estimation$random.cov)

cat("\nModel Information:\n")
print(Model.Estimation$model.info)

cat("\nEstimate of Ri: \n")
print(Model.Estimation$Ri, quote = FALSE)

options(original)



###################################################
### code chunk number 3: Simulation 2： ARMA Correlation Structure
###################################################
library(BayesRGMM)
rm(list=ls(all=TRUE))


Fixed.Effs = c(-0.2,-0.8, 1.0, -1.2)
P = length(Fixed.Effs)
q = 1
T = 10
N = 100
num.of.iter = 100

ARMA.sim.data = SimulatedDataGenerator(Num.of.Obs = N, Num.of.TimePoints = T, 
  Fixed.Effs = Fixed.Effs, Random.Effs = list(Sigma = 0.5*diag(1), df=3), 
  Cor.in.DesignMat = 0., list(Missing.Mechanism = 2, RegCoefs = c(-1.5, 1.2)), 
  Cor.Str = "ARMA", ARMA.para=list(AR.para = 0.4, MA.para=0.2))

ARMA.output = BayesRobustProbit(fixed = as.formula(paste("y~-1+", 
  paste0("x", 1:P, collapse="+"))), data=ARMA.sim.data$sim.data, random = ~ 1, 
  subset = NULL, na.action='na.exclude', arma.order = c(1, 1), 
  num.of.iter = num.of.iter)

original = options(digits = 4)

Model.Estimation = BayesRobustProbitSummary(ARMA.output)

cat("\nCoefficients:\n")
print(Model.Estimation$beta.est.CI)

cat("\nAMRA parameters:\n\n")
print(Model.Estimation$arma.est)

cat("\nRandom effect: \n")
print(Model.Estimation$random.cov)

cat("\nModel Information:\n")
print(Model.Estimation$model.info)

options(original)




