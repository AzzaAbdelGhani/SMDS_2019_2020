## ----global_options, include=FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10), fig.path = 'figs/')
library(MASS)

## ----setup, include=FALSE--------------------------------------------------------------------------------------------
library(knitr)
local({
  hook_plot = knit_hooks$get('plot')
  knit_hooks$set(plot = function(x, options) {
    paste0('\n\n----\n\n', hook_plot(x, options))
  })
})


## ----data, echo =FALSE-----------------------------------------------------------------------------------------------
library(DAAG)
library(kableExtra)

options(knitr.table.format = "html") 
n <- length(nihills[1:8,1])
data_frame <- round(nihills[1:n,],3)

nameshills <- row.names(nihills)
dimnames(data_frame)<- list(nameshills[1:n], c( "dist","climb", 
                            "time","timef"))
kable(data_frame, "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)



## ----plot------------------------------------------------------------------------------------------------------------
dist <- nihills$dist
time <- nihills$time
n <- length(time)
plot(dist, time)


## ----fit one---------------------------------------------------------------------------------------------------------
nihills.lm.one <- lm(time ~ dist, data = nihills)
summary(nihills.lm.one)
par(mfrow=c(2,2))
plot(nihills.lm.one)
par(mfrow=c(1,1))
plot(dist, time)
lines(dist, as.vector(nihills.lm.one$fitted.values), col ="red")
text(13,3, "y = b0+b1*dist", col="red")

regr_ss <- sum( (as.vector(nihills.lm.one$fitted.values)-mean(time))^2)
tot_ss <- sum( (time-mean(time))^2)
R_sq <- regr_ss/tot_ss
R_sq


## ----diagn, echo =TRUE-----------------------------------------------------------------------------------------------

library(DAAG)
library(lattice)

par(mfrow=c(1,2))
#scatterplot matrix with original data
splom( ~nihills[, c("dist", "climb", "time")], cex.labels=1.2,
      varnames= c("dist\n(miles)", "climb\n(feet", "time\n(hours"))



## ----first-----------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2), oma =c(0,0,0,0))
nihills.lm <- lm(time ~  dist + climb, data = nihills)
plot(nihills.lm)


## ----summary---------------------------------------------------------------------------------------------------------
summary(nihills.lm, corr=TRUE)
coef <- summary(nihills.lm)$coef



## ----second----------------------------------------------------------------------------------------------------------

nihills.lm.noint <- lm(time ~ -1 + dist + climb, data = nihills)
summary(nihills.lm.noint)


## ----hist------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
hist(time, breaks=30, xlab="time", probability =TRUE, xlim =c(0,4))
hist(log(time), breaks=8, xlab="time", probability =TRUE, xlim =c(-1.6,1.6))



## ----log lm----------------------------------------------------------------------------------------------------------
splom( ~log(nihills[, c("dist", "climb", "time")]), cex.labels=1.2,
      varnames= c("dist\n(miles)", "climb\n(feet", "time\n(hours"))
nihills.log.lm <- lm(log(time)~log(dist)+log(climb), data =nihills)
par(mfrow=c(2,2))
plot(nihills.log.lm)
summary(nihills.log.lm)




## ----criterion-------------------------------------------------------------------------------------------------------
#AIC
c(extractAIC(nihills.lm.one)[2], extractAIC(nihills.lm)[2],
extractAIC(nihills.log.lm)[2])
#BIC
c(extractAIC(nihills.lm.one, k =log(n))[2], extractAIC(nihills.lm, k =log(n))[2],
extractAIC(nihills.log.lm, k=log(n))[2])




## ----pred------------------------------------------------------------------------------------------------------------

coverage <- predict(nihills.lm, interval ="confidence")
predict(nihills.lm, interval ="confidence")[1:5,]
freq_coverage <- sum( time >= coverage[,2] & time <= coverage[,3])
freq_coverage <- freq_coverage/n
freq_coverage

#log model
coverage_log <- exp(predict(nihills.log.lm, interval ="confidence"))
exp(predict(nihills.log.lm, interval ="confidence"))[1:5,]
freq_coverage_log <- sum( time>=coverage_log[,2] & time <=coverage_log[,3])
freq_coverage_log <- freq_coverage_log/n
freq_coverage_log



## ----anova-----------------------------------------------------------------------------------------------------------
anova(nihills.log.lm)



## ----vif-------------------------------------------------------------------------------------------------------------
vif(nihills.log.lm)



## ----vif sum---------------------------------------------------------------------------------------------------------
climb <- nihills$climb
dist <- nihills$dist
logsum <- 4 + 3*log(dist)-2*climb 
nihills.log.lm.sum <- lm(log(time)~log(dist)+log(climb)+ logsum,
  data=nihills)
summary(nihills.log.lm.sum)
vif(nihills.log.lm.sum)


## ----vif ind---------------------------------------------------------------------------------------------------------
logsum <- log(dist)+log(climb)
logdiff <- log(dist)-log(climb)
nihills.log.lm.ind <- lm(log(time)  ~ logsum + logdiff, data=nihills)
vif(nihills.log.lm.ind)


## ----ridge-----------------------------------------------------------------------------------------------------------
library(MASS)

#lambda =0 

nihills.ridge <- lm.ridge( log(time)~log(dist)+log(climb), 
  data=nihills, lambda =0 )
nihills.ridge

#select lambda in terms of GCV error

select(lm.ridge(log(time)~log(dist)+log(climb)+logsum, data=nihills,
              lambda = seq(0,1,0.001)))

nihills.ridge.sel <- lm.ridge( log(time)~log(dist)+log(climb), 
  data=nihills, lambda =0.371  )
nihills.ridge.sel



## ----lasso-----------------------------------------------------------------------------------------------------------
library(lasso2)
nihills.lasso <- l1ce(log(time)~log(dist)+log(climb), data =nihills,
  absolute.t =  TRUE)
summary(nihills.lasso)$coefficients



## ----poiss-----------------------------------------------------------------------------------------------------------
y <- c(4, 27, 2, 21, 8, 3, 3, 0, 4, 4, 5, 2, 2, 5, 7, 1, 5, 0, 15, 36)
age <- c(54 ,26,41,23, 41, 75, 28, 78, 44, 52, 50, 28, 42, 50, 31, 59, 33, 35, 26, 22)
family <- c(1,1,0,1,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1)
earning <- c( 1.2, 5.3, 2.4, 6.3, 4.2, 1.5, 1.1, 0.9, 1.2, 1.1, 1.3, 1.5, 4.2, 1.3, 1.2, 1.3, 1.9, 2.7, 7.3, 5)
n <- length(y)
calls.glm <- glm(y ~ age+family+earning, family = poisson)
summary(calls.glm)
par(mfrow=c(1,2))
plot(calls.glm, which=c(1,3))


## ----quasi poiss-----------------------------------------------------------------------------------------------------

calls.quasi.glm <- glm(y ~ age+family+earning, family = quasipoisson)
summary(calls.quasi.glm)
par(mfrow=c(1,2))
plot(calls.quasi.glm, which=c(1,3))



