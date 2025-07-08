require(frailtypack) #requires frailtypack version >=3.6.4

#1st illustration

data(gastadj)
gastadj$timeS<-gastadj$timeS/365
gastadj$timeT<-gastadj$timeT/365
##"statusS" corresponds now
##to a time-to-relapse event
gastadj[gastadj$timeS==gastadj$timeT &
          gastadj$statusS==1,c("statusS")]<-0
##select 20% of the original dataset
set.seed(1)
n<-nrow(gastadj)
subset<-gastadj[sort(sample(1:nrow(gastadj),
                            round(n*0.2),replace = F)),]
## Model fitting
mod.gast<-jointSurroPenal(subset,n.knots = 4,
                          indicator.zeta = 0,
                          indicator.alpha = 0,
                          mediation=TRUE,g.nknots=1,
                          pte.times=seq(1.5,2,length.out=30),pte.nmc=10000,
                          pte.boot=TRUE,pte.nboot=1000,
                          pte.boot.nmc=1000)

##Results
summary(mod.gast)
plot(mod.gast,type.plot="Hazard",plot.mediation="All")

#2nd illustration
data(colorectal)
data(colorectalLongi)

colorectalSurv <- subset(colorectal, new.lesions == 0)
colorectalSurv$treatment<-sapply(colorectalSurv$treatment,
                                 function(t) ifelse(t=="S",1,0))
colorectalLongi$treatment<-sapply(colorectalLongi$treatment,
                                  function(t) ifelse(t=="S",1,0))

##Model fitting
mod.col=longiPenal(Surv(time1, state) ~ age+treatment,
                   tumor.size ~ age+year*treatment,
                   data = colorectalSurv, data.Longi = colorectalLongi,
                   random = c("1", "year"), id = "id",
                   link = "Current-level", timevar = "year",
                   method.GH = "Pseudo-adaptive",
                   mediation = TRUE,
                   med.trt = colorectalSurv$treatment,
                   med.center = NULL,
                   n.knots = 7,
                   kappa = 2,
                   pte.times = seq(1,2,length.out = 30),
                   pte.boot = TRUE,
                   pte.nboot = 2000,
                   pte.nmc = 1000)


##Results
print(mod.col)
plot(mod.col,plot.mediation='All',conf.bands = T)
