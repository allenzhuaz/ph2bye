library(haven)
mydata <- read_sas("adrs.sas7bdat")
unique(mydata$PARAM)
## there are 53 kinds of parameter
response <- subset(mydata, PARAM=="Recist Best Overall Response")
unique(response$AVALC)
p.PD <- sum(response$AVALC=="PD")/length(response$AVALC)
p.SD <- sum(response$AVALC=="SD")/length(response$AVALC)
p.uPR <- sum(response$AVALC=="uPR")/length(response$AVALC)
p.PR <- sum(response$AVALC=="PR")/length(response$AVALC)
p.CR <- sum(response$AVALC=="CR")/length(response$AVALC)
# p.NA <- sum(response$AVALC=="NA")/length(response$AVALC)
r01 <- ifelse(response$AVALC%in% c("uPR","PR","CR"), 1, 0)
RR <- mean(resp)
RR <- sum(response$AVALC %in% c("uPR","PR","CR") )/length(response$AVALC)

KT.mean <- 0.236; KT.n <- 173; KT.var <- KT.mean*(1-KT.mean)/KT.n
OD.mean <- 0.317; OD.n <- 120; OD.var <- OD.mean*(1-OD.mean)/OD.n

betaprior <- function(m,v){
  a <- m^2*(1-m)/v - m; b <- m*(1-m)^2/v - (1-m)
  return(c(a=a,b=b))
}

binom.test(x=round(KT.n*KT.mean),n = KT.n, p=0.2)
binom.test(x=round(OD.n*OD.mean),n = OD.n, p=0.2)

KT.prior <- betaprior(KT.mean, KT.var)
OD.prior <- betaprior(OD.mean, OD.var)


ani.options(convert='C:\\Users\\yalin.zhu\\Documents\\ImageMagick-7.0.1-10-portable-Q16-x64\\convert.exe')

library(ph2bye)
BB.aniplot(a = KT.prior[1], b = KT.prior[2], r = r01, N = 1, time.interval = 0.01, output = F)
BB.aniplot(a = OD.prior[1], b = OD.prior[2], r = r01, N = 1, time.interval = 0.01, output = F)
bayes.design(a = KT.prior[1], b = KT.prior[2], r = r01, stop.rule = "futility", p0 = 0.2, time.interval = 0.01)
bayes.design(a = KT.prior[1], b = KT.prior[2], r = r01, stop.rule = "efficacy", p0 = 0.2, time.interval = 0.01)

bayes.design(a = OD.prior[1], b = OD.prior[2], r = r01, stop.rule = "futility", p0 = 0.2, time.interval = 0.01)
bayes.design(a = OD.prior[1], b = OD.prior[2], r = r01, stop.rule = "efficacy", p0 = 0.2, time.interval = 0.01)

