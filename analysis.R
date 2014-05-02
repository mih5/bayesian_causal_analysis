

#------------------------LOAD DATA-----------------------------------------
 
fev.org <- read.csv("C:/Users/MaoHu/Dropbox/Junior Year/STA 320/Project/fev.csv")

#remove missing cases
fev = fev.org[complete.cases(fev.org),]

#decision: I will restrict analysis to persons whose
# "proportion of covered days" using medication is above .8
# the reason for this decision is a desire to control for
# prescription medicine.  Unfortunately, proportion of covered days
# cannot be used as a covariate - it is not clear if it is collected before
# the "treatment assignment"- use of alternative medicine
# one way around this is to restrict analysis to only one group
# those who take controller medicine for 0.8 of the last 12 months or more


summary(fev)

fevdata = subset(fev,fev$controller_medicine_use=="PDC>=0.8",)

#remove outcomes
W = fevdata$cam_use=="Yes"
Y = fevdata$fev1
X = data.frame(cbind(fevdata$age,fevdata$ethnicity=="White",
                     fevdata$household_income==">=$6,000",
                     fevdata$sex=="Male",fevdata$education!="Less than college"))
names(X)=c("age","ethnicity","income","sex","education")


source("http://stat.duke.edu/courses/Spring14/sta320.01/CausalInference.R")
require(mvtnorm)

#------------------------FUNCTIONS-------------------------------------

plotPropensityScores = function(ps,W,main=NULL){
  require(ggplot2)
  ##Plot Propensity Scores using ggplot2 
  ##dataplot not found?
  dataplot=data.frame(cbind(ps,W))
  names(dataplot)=c("ps","treat")
  print(summary(dataplot))
  p=ggplot(dataplot)
  p=p+geom_density(aes(x=ps,color=as.factor(treat),group=as.factor(treat)))
  if(!is.null(main)){
    p=p+labs(title=main) 
  }
  else{
    p=p+labs(title="Density Plot of Propensity Scores")  
  }  
  if (sum(ps<0)>0 | sum (ps>1) >0){
    p=p+labs(x="Linear Propensity Score (log odds)")
  }
  else{
    p=p+labs(x="Propensity Score")
  }
  p=p+labs(color="Treatment Group")
  p
}

set.seed(835892)


Gibbs_Sampler = function(steps,y_obs,W,X,g,nu_t0,sigma_t0,nu_c0,sigma_c0,beta_c,sigma_c,beta_t,sigma_t){
  # Gibbs_Sampler: a function which samples from the 
  # posterior conditional distributions of the parameters iteratively
  # and from the posterior predictive of the missing potential outcomes
  # Input arguments:
  #
  #   steps: number of steps sampler should run
  #
  #   y_obs: observed outcomes
  #   W: assignment vector
  #   X: matrix of covariates
  #
  #   g: hyperparameter
  #   nu_t0: hyperparameter
  #   sigma_t0: hyperparameter
  #   nu_c0: hyperparameter
  #   sigma_c0: hyperparameter
  #   
  #   beta_c: initial value for regression coefficients for control units
  #   sigma_c: intial value for standard deviation for control units
  #   beta_t: initial value for regression coefficients for treated units
  #   sigma_t: intial value for standard deviation for treated units
  #
  # Output arguments:
  #   SAMPLES: a matrix where each row is one "round" of draws (one iteration)
  #
  
  require(mvtnorm)
  SAMPLES = NULL
  
  X.mean = apply(X,MARGIN=2,mean)
  
  n = length(W)
  n_c = sum(W==0)
  n_t = n-n_c
  
  #subset the outcomes by treatment assignment
  y_c = y_obs[W==0,drop=FALSE]
  y_t = y_obs[W==1,drop=FALSE]
  
  #subset the covariates by treatment assignment
  X_c = X[W==0,,drop=FALSE]
  X_t = X[W==1,,drop=FALSE]
  
  #updated hyperparameters for coefficients b_c, b_t
  
  lambda_c0_star = g/(g+1)*solve(t(X_c)%*%X_c)
  beta_c0_star = lambda_c0_star%*%t(X_c)%*%y_c
  
  lambda_t0_star = g/(g+1)*solve(t(X_t)%*%X_t)
  beta_t0_star = lambda_t0_star%*%t(X_t)%*%y_t
  
  for (s in 1:steps){
    print(s)
    #sample control parameters
    
    #sample beta_c
    beta_c = t(rmvnorm(1,beta_c0_star,lambda_c0_star*sigma_c^2))
    
    #sample sigma_c
    SSR_c = t(y_c)%*%y_c - 2*t(beta_c)%*%t(X_c)%*%y_c+t(beta_c)%*%t(X_c)%*%X_c%*%beta_c
    sigma_c = sqrt(as.numeric(1/rgamma(1,(nu_c0+n_c)/2,(nu_c0*sigma_c0^2+SSR_c)/2)))
    
    #sample treatment parameters
    
    #sample beta_t
    beta_t = t(rmvnorm(1,beta_t0_star,lambda_t0_star*sigma_t^2))
    
    #sample sigma_t_star
    SSR_t = t(y_t)%*%y_t - 2*t(beta_t)%*%t(X_t)%*%y_t+t(beta_t)%*%t(X_t)%*%X_t%*%beta_t
    sigma_t = sqrt(as.numeric(1/rgamma(1,(nu_t0+n_t)/2,(nu_t0*sigma_t0^2+SSR_t)/2)))
    
    #sample from conditional distribution of y_miss - imputation
    y_miss_mean = W*(X)%*%beta_c+(1-W)*(X)%*%beta_t
    y_miss_var = (W*sigma_c^2+(1-W)*sigma_t^2)*diag(1,n)
    y_miss = rmvnorm(1,y_miss_mean,y_miss_var)
    
    tau=X.mean%*%beta_t-X.mean%*%beta_c
    
    store = cbind(tau,t(beta_c),t(beta_t),sigma_c,sigma_t)
    SAMPLES = rbind(SAMPLES,store) 
    
  }
  
  
  return(SAMPLES)
  
}

##trimIterate
##function that trims the dataset so that there are no controls above treated and no treated below controls
##using linear propensity scores
trimIterate = function(model,data,W){
  ##trimIterate takes the following parameters:
  ##model: a model object from variable section
  ##data: the data without the outcomes
  ##W: the assignment vector
  
  condition = TRUE
  iter = 0
  
  while(condition){
    
    iter = iter +1
    print(iter)
    
    ## get linear propensity scores
    ps=predict(model)
    
    overlap(ps,W)
    plot.ps(ps,W)
    
    ##get which observations are controls below treated or treated above controls
    controls.below.treated = ps < min(ps[W==1])
    treated.above.controls = ps > max(ps[W==0])
    
    
    ##check condition
    if((sum(controls.below.treated)+sum(treated.above.controls)<1)||iter>10){
      condition=FALSE
    }
    else{
      ##trim data and outcomes
      exclude=(1-controls.below.treated)&(1-treated.above.controls)
      data = data[exclude,]
      W = W[exclude]
      
      ##refit model
      model = glm(formula(model),data=data,family="binomial")
    }
    
  }
  
  return(data)
}



#--------------------"Design" Phase---------------------------
summary(X)

#visualize covariate balance
cov.balance(X,W)
boxplot(X$age~W)
#covariate balance is not terrible
#t-statistics are less than 2 in absolute magnitude
#however, we can still improve

# ESTIMATE PROPENSITY SCORES
#logistic regression with higher order terms for quantitative variables
#and second-order interactions
age_sq=X$age^2
designmatrix = data.frame(cbind(W,X,age_sq))

#not substantial evidence of strong predictors of assignment
model.full = glm(W~.+.^2,family="binomial",data=designmatrix)
summary(model.full)

#fit a stepwise regression to obtain propensity score model
model.step = step(model.full,scope=list(lower=~age+ethnicity+income+sex+education,upper=~.))
summary(model.step)

#I will use linearized propensity scores
ps = predict(model.step)
ps.r = predict(model.step,type="response")

#the balance doesn't look that bad
#not suprising since covariate balance was alright
plotPropensityScores(ps,W,main="Before Trimming")
overlap(ps,W)
plotPropensityScores(ps.r,W)
overlap(ps,W)
# will need to trim the data

#trimming data
fevdata2 = data.frame(cbind(Y,designmatrix,W))

fevdata_trimmed = trimIterate(model.step,fevdata2,W)

#data trimmed
Y=fevdata_trimmed$Y
X=fevdata_trimmed[,c(3:7)]
W=fevdata_trimmed$W

age_sq=X$age^2
designmatrix2 = data.frame(cbind(W,X,age_sq))

model.step2=glm(formula(model.step),data=designmatrix2,family="binomial")
ps = predict(model.step2)
plotPropensityScores(ps,W,main="After Trimming")
#still, for the sake of improvement, I will subclassify
#based on propensity scores


#break units by propensity score quantiles
breaks.2=quantile(ps,c(0.5))
subclass.2 = subclasses(ps,breaks.2)
cov.balance(X,W,subclass.2)

breaks.3=quantile(ps,c(1/3,2/3))
subclass.3 = subclasses(ps,breaks.3)
cov.balance(X,W,subclass.3)

breaks.4=quantile(ps,c(0.25,0.5,0.75))
subclass.4 = subclasses(ps,breaks.4)
cov.balance(X,W,subclass.4)

#best covariance balanced obtained by 2 breaks

#----------------------------ANALYSIS-----------------------------------

pairs(cbind(Y,X))

#Neyman Super Population Analysis

#estimate for population average treatment effect
tau_neyman_1=(mean(Y[subclass.2==1 & W==1])-mean(Y[subclass.2==1 & W==0]))
tau_neyman_2=(mean(Y[subclass.2==2 & W==1])-mean(Y[subclass.2==2 & W==0]))


tau_neyman=tau_neyman_1*96/191+tau_neyman_2*95/191
tau_neyman

#estimate for variance of population average treatment effect
var_y_1_control = 1/(sum(1-W[subclass.2==1])*(sum(1-W[subclass.2==1])-1))*sum((Y[subclass.2==1& W==0]-mean(Y[subclass.2==1 & W==0]))^2)
var_y_1_treat = 1/(sum(W[subclass.2==1])*(sum(W[subclass.2==1])-1))*sum((Y[subclass.2==1& W==1]-mean(Y[subclass.2==1 & W==1]))^2)
var_y_1_tau = var_y_1_control+var_y_1_treat

var_y_2_control = 1/(sum(1-W[subclass.2==1])*(sum(1-W[subclass.2==1])-1))*sum((Y[subclass.2==2& W==0]-mean(Y[subclass.2==2 & W==0]))^2)
var_y_2_treat = 1/(sum(W[subclass.2==1])*(sum(W[subclass.2==1])-1))*sum((Y[subclass.2==2& W==1]-mean(Y[subclass.2==2 & W==1]))^2)
var_y_2_tau = var_y_2_control+var_y_2_treat

var_tau = (96/191)^2*var_y_1_tau+(95/191)^2*var_y_2_tau
var_tau

tau_neyman+2*sqrt(var_tau)
tau_neyman-2*sqrt(var_tau)

#bayesian model imputation analysis
model=lm(fev1~cam_use+age+ethnicity+household_income+sex+education,data=fevdata)
sigmaols=as.numeric(summary(model)[6])

#model imputation, no subclassification
samples = Gibbs_Sampler(1000,Y,W,as.matrix(X),
                        g=205,
                        nu_t0=1,
                        sigma_t0=sigmaols,
                        nu_c0=1,
                        sigma_c0=sigmaols,
                        beta_c=t(t(c(0,0,0,0,0))),
                        sigma_c=1,
                        beta_t=t(t(c(0,0,0,0,0))),
                        sigma_t=1)

#check for good mixing
impute.mcmc = mcmc(samples,start=5)
plot(impute.mcmc[,2])
autocorr.plot(impute.mcmc[,2:6])
summary(impute.mcmc)

#estimate of population average treatment effect
mean(impute.mcmc[,1])
sd(impute.mcmc[,1])

#model imputation, with subclassification

samples1 = Gibbs_Sampler(1000,Y[subclass.2==1],W[subclass.2==1],as.matrix(X[subclass.2==1,]),
                        g=95,
                        nu_t0=1,
                        sigma_t0=sigmaols,
                        nu_c0=1,
                        sigma_c0=sigmaols,
                        beta_c=t(t(c(0,0,0,0,0))),
                        sigma_c=1,
                        beta_t=t(t(c(0,0,0,0,0))),
                        sigma_t=1)

samples2 = Gibbs_Sampler(1000,Y[subclass.2==2],W[subclass.2==2],as.matrix(X[subclass.2==2,]),
                         g=96,
                         nu_t0=1,
                         sigma_t0=sigmaols,
                         nu_c0=1,
                         sigma_c0=sigmaols,
                         beta_c=t(t(c(0,0,0,0,0))),
                         sigma_c=1,
                         beta_t=t(t(c(0,0,0,0,0))),
                         sigma_t=1)


impute.mcmc.1 = mcmc(samples1)
autocorr.plot(impute.mcmc.1[,2:6])
summary(impute.mcmc.1)
tau_impute_1=mean(impute.mcmc.1[,1])
tau_var_impute_1=sd(impute.mcmc.1[,1])^2

impute.mcmc.2 = mcmc(samples2)
autocorr.plot(impute.mcmc.2[,2:6])
summary(impute.mcmc.2)
tau_impute_2=mean(impute.mcmc.2[,1])
tau_var_impute_2=sd(impute.mcmc.2[,1])^2

#distribution of the population average treatment effect
tau_impute =impute.mcmc.1[,1]*96/191+impute.mcmc.2[,1]*95/191

#population average treatment effect
mean(tau_impute)

#variance of PATE
sqrt(var(tau_impute))


qplot(as.numeric(tau_impute),geom="auto",xlab="Estimate of tau")
hist(tau_impute)
quantile(tau_impute,c(.025,.975))

#distribution of treatment effect in subclass 1
#taking alternative medicine increases capacity?
qplot(as.numeric(impute.mcmc.1[,1]))

#distribution of treatment effect in subclass 2
#taking alternative decreases capacity?
qplot(as.numeric(impute.mcmc.2[,1]))



summary(X[subclass.2==1,])
summary(X[subclass.2==2,])


summary(Y)
sd(Y)
summary(X)
summary(as.numeric(W))

sd(X$age)
