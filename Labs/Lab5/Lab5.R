## ----global_options, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10), fig.path = 'figs/',options(width=400))
library(MASS)

## ----setup, include=FALSE--------------------------------------------------------------------------------
library(knitr)
local({
  hook_plot = knit_hooks$get('plot')
  knit_hooks$set(plot = function(x, options) {
    paste0('\n\n----\n\n', hook_plot(x, options))
  })
})


## ----load data-------------------------------------------------------------------------------------------
data <- read.csv2("Wells.csv", header=TRUE)
str(data)
summary(data)


## ----explorative plots, echo=FALSE-----------------------------------------------------------------------
par(mfrow=c(1,2))
boxplot(data$dist ~ data$switch, horizontal = TRUE)
boxplot(data$arsenic ~ data$switch, horizontal = TRUE)
par(mfrow=c(1,1))
barplot(prop.table(table(data$switch, data$educ), margin = 1),
        beside=TRUE, legend.text = c("no switch", "switch"))


## ----logistic one predictor------------------------------------------------------------------------------
## Logistic regression with one predictor

fit.1 <- glm(switch ~ dist, family=binomial(link="logit"), data=data)
summary(fit.1)


## ----logistic one coeff intercept------------------------------------------------------------------------
invlogit <- function (x) {1/(1+exp(-x))}
invlogit(fit.1$coefficients%*%c(1,0))


## ----logistic coeff interpretation-----------------------------------------------------------------------
# difference on the probability to switch well due to an increase of one unit on the predictor from the mean
invlogit(fit.1$coefficients%*%c(1,mean(data$dist)+1))-
  invlogit(fit.1$coefficients%*%c(1,mean(data$dist)))


## ----logistic coeff interpretation2----------------------------------------------------------------------
# difference on the probability to switch well due to an increase of one unit on the predictor from the 99 percentile
invlogit(fit.1$coefficients%*%c(1,quantile(data$dist,.99)+1))-
  invlogit(fit.1$coefficients%*%c(1,quantile(data$dist,.99)))


## ----logistic one predictor transf-----------------------------------------------------------------------
## Repeat the regression above with distance in 100-meter units
data$dist100 <- data$dist/100
fit.2 <- glm (switch ~ dist100, family=binomial(link="logit"), data=data)
summary(fit.2)

# difference on the probability scale due to an increase of 1 unit on the predictor (dist100) from the mean
invlogit(fit.2$coefficients%*%c(1,mean(data$dist100)+1))-
  invlogit(fit.2$coefficients%*%c(1,mean(data$dist100)))

#notice that
# difference on the probability scale due to an increase of 100 unit on the predictor from the mean (fit.1)
invlogit(fit.1$coefficients%*%c(1,mean(data$dist)+100))-
  invlogit(fit.1$coefficients%*%c(1,mean(data$dist)))



## ----logistic one predictor graph, echo=FALSE------------------------------------------------------------
## Graphing the fitted model with one predictor

jitter.binary <- function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

switch.jitter <- jitter.binary(data$switch)

plot(data$dist, switch.jitter, xlab="Distance (in meters) to nearest safe well", 
     ylab="Pr (switching)", type="n", xaxs="i", yaxs="i", mgp=c(2,.5,0))
curve (invlogit(coef(fit.1)[1]+coef(fit.1)[2]*x), lwd=1, add=TRUE)
points (data$dist, jitter.binary(data$switch), pch=20, cex=.1)


## ----logistic arsenic, echo=FALSE------------------------------------------------------------------------
## Histogram on arsenic levels (Figure 5.10)

hist (data$arsenic, breaks=seq(0,.25+max(data$arsenic[!is.na(data$arsenic)]),.25), 
      freq=TRUE, xlab="Arsenic concentration in well water", 
      ylab="", main="", mgp=c(2,.5,0))

## Logistic regression with second input variable

fit.3 <- glm(switch ~dist100 + arsenic, family=binomial(link="logit"), data=data)
summary(fit.3)


## ----logistic arsenic graph, echo=FALSE------------------------------------------------------------------
## Graphing the fitted model with two predictors (Figure 5.11)

plot(data$dist, switch.jitter, xlim=c(0,max(data$dist)), 
     xlab="Distance (in meters) to nearest safe well", ylab="Pr (switching)", 
     type="n", xaxs="i", yaxs="i", mgp=c(2,.5,0))
curve(invlogit(cbind(1, x/100, .5) %*% coef(fit.3)), lwd=.5, add=TRUE)
curve(invlogit(cbind(1, x/100, 1.0) %*% coef(fit.3)), lwd=.5, add=TRUE)
points(data$dist, jitter.binary(data$switch), pch=20, cex=.1)
text (50, .27, "if As = 0.5", adj=0, cex=.8)
text (75, .50, "if As = 1.0", adj=0, cex=.8)

plot(data$arsenic, switch.jitter, xlim=c(0,max(data$arsenic)), 
     xlab="Arsenic concentration in well water", ylab="Pr (switching)", 
     type="n", xaxs="i", yaxs="i", mgp=c(2,.5,0))
curve(invlogit(cbind(1, 0, x) %*% coef(fit.3)), lwd=.5, add=TRUE)
curve(invlogit(cbind(1, 0.5, x) %*% coef(fit.3)), lwd=.5, add=TRUE)
points(data$arsenic, jitter.binary(data$switch), pch=20, cex=.1)
text(1.5, .78, "if dist = 0", adj=0, cex=.8)
text(2.2, .6, "if dist = 50", adj=0, cex=.8)



## ----logistic with interaction---------------------------------------------------------------------------
# Logistic regression with interactions

fit.4 <- glm (switch ~ dist100 + arsenic + dist100:arsenic, #switch ~ dist100*arsenic equivalently
  family=binomial(link="logit"), data=data)
summary(fit.4)


## --------------------------------------------------------------------------------------------------------
# constant term
invlogit(c(1,mean(data$dist100),mean(data$arsenic),
           mean(data$dist100)*mean(data$arsenic))%*%coefficients(fit.4))


## --------------------------------------------------------------------------------------------------------
# an increse of 100 meters in the distance, helding arsenic at its mean
invlogit(c(1,mean(data$dist100)+1,mean(data$arsenic),
           (mean(data$dist100)+1)*mean(data$arsenic))%*%coefficients(fit.4))-
invlogit(c(1,mean(data$dist100),mean(data$arsenic),
           mean(data$dist100)*mean(data$arsenic))%*%coefficients(fit.4))

# an increse of 1 in the arsenic, helding distance at its mean
invlogit(c(1,mean(data$dist100),mean(data$arsenic)+1,
           (mean(data$dist100))*(mean(data$arsenic)+1))%*%coefficients(fit.4))-
invlogit(c(1,mean(data$dist100),mean(data$arsenic),
           mean(data$dist100)*mean(data$arsenic))%*%coefficients(fit.4))


## ----logistic with interaction coeff---------------------------------------------------------------------
coef(fit.4)[2]+coef(fit.4)[4]*mean(data$arsenic)
coef(fit.4)[3]+coef(fit.4)[4]*mean(data$dist100)


## ----logistic with interaction center--------------------------------------------------------------------
## Centering the input variables

data$c.dist100 <- scale(data$dist100, scale = FALSE)
data$c.arsenic <- scale(data$arsenic, scale = FALSE)

## Refitting the model with centered inputs

fit.5 <- glm (switch ~ c.dist100 + c.arsenic + c.dist100:c.arsenic,
  family=binomial(link="logit"), data=data)
summary(fit.5)


## ----logistic with interaction graph, echo=FALSE---------------------------------------------------------
## Graphing the model with interactions

plot(data$dist, switch.jitter, xlim=c(0,max(data$dist)), xlab="Distance (in meters) to nearest safe well", 
   ylab="Pr (switching)", type="n", xaxs="i", yaxs="i", mgp=c(2,.5,0))
curve (invlogit(cbind (1, x/100, .5, .5*x/100) %*% coef(fit.4)), lwd=.5, add=TRUE)
curve (invlogit(cbind (1, x/100, 1.0, 1.0*x/100) %*% coef(fit.4)), lwd=.5, add=TRUE)
points (data$dist, jitter.binary(data$switch), pch=20, cex=.1)
text (50, .37, "if As = 0.5", adj=0, cex=.8)
text (75, .50, "if As = 1.0", adj=0, cex=.8)

plot(data$arsenic, switch.jitter, xlim=c(0,max(data$arsenic)), xlab="Arsenic concentration in well water",
   ylab="Pr (switching)", type="n", xaxs="i", yaxs="i", mgp=c(2,.5,0))
curve (invlogit(cbind (1, 0, x, 0*x) %*% coef(fit.4)), lwd=.5, add=TRUE)
curve (invlogit(cbind (1, 0.5, x, 0.5*x) %*% coef(fit.4)), lwd=.5, add=TRUE)
points (data$arsenic, jitter.binary(data$switch), pch=20, cex=.1)
text (1.5, .8, "if dist = 0", adj=0, cex=.8)
text (2.2, .6, "if dist = 50", adj=0, cex=.8)


## ----logistic with other predictors----------------------------------------------------------------------
## Adding social predictors

data$c.educ4 <- scale(data$educ/4, scale=FALSE)
fit.7 <- glm (switch ~ c.dist100 + c.arsenic + c.educ4 + c.dist100:c.arsenic, family=binomial(link="logit"), data=data)
summary(fit.7)


## ----logistic with other predictors and interactions-----------------------------------------------------
## Adding further interactions (centering education variable)

fit.8 <- glm (switch ~ c.dist100 + c.arsenic + c.educ4 + c.dist100:c.arsenic +
                c.dist100:c.educ4 + c.arsenic:c.educ4, family=binomial(link="logit"), data=data)
summary(fit.8)


## ----residuals-------------------------------------------------------------------------------------------
## Residual Plot

pred.8 <- fit.8$fitted.values

plot(c(0,1), c(-1,1), xlab="Estimated  Pr (switching)", ylab="Observed - estimated", type="n", main="Residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
points (pred.8, data$switch-pred.8, pch=20, cex=.2)


## ----binned residuals------------------------------------------------------------------------------------
### Binned residual Plot 

 ## Defining binned residuals

binned.resids <- function (x, y, nclass=sqrt(length(x))){
  breaks.index <- floor(length(x)*(1:(nclass-1))/nclass)
  breaks <- c (-Inf, sort(x)[breaks.index], Inf)
  output <- NULL
  xbreaks <- NULL
  x.binned <- as.numeric (cut (x, breaks))
  for (i in 1:nclass){
    items <- (1:length(x))[x.binned==i]
    x.range <- range(x[items])
    xbar <- mean(x[items])
    ybar <- mean(y[items])
    n <- length(items)
    sdev <- sd(y[items])
    output <- rbind (output, c(xbar, ybar, n, x.range, 2*sdev/sqrt(n)))
  }
  colnames (output) <- c ("xbar", "ybar", "n", "x.lo", "x.hi", "2se")
  return (list (binned=output, xbreaks=xbreaks))
}

 ## Binned residuals vs. estimated probability of switching

br.8 <- binned.resids (pred.8, data$switch-pred.8, nclass=40)$binned
plot(range(br.8[,1]), range(br.8[,2],br.8[,6],-br.8[,6]), xlab="Estimated  Pr (switching)", ylab="Average residual", type="n", main="Binned residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
lines (br.8[,1], br.8[,6], col="gray", lwd=.5)
lines (br.8[,1], -br.8[,6], col="gray", lwd=.5)
points (br.8[,1], br.8[,2], pch=19, cex=.5)


## ----binned residuals vs predictors----------------------------------------------------------------------
## Plot of binned residuals vs. inputs of interest

  # distance

br.dist <- binned.resids(data$dist, data$switch-pred.8, nclass=40)$binned
plot(range(br.dist[,1]), range(br.dist[,2],br.dist[,6],-br.dist[,6]), xlab="Distance to nearest safe well", ylab="Average residual", type="n", main="Binned residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
lines (br.dist[,1], br.dist[,6], col="gray", lwd=.5)
lines (br.dist[,1], -br.dist[,6], col="gray", lwd=.5)
points (br.dist[,1], br.dist[,2], pch=19, cex=.5)

  # arsenic

br.arsenic <- binned.resids(data$arsenic, data$switch-pred.8, nclass=40)$binned
plot(range(0,br.arsenic[,1]), range(br.arsenic[,2],br.arsenic[,6],-br.arsenic[,6]), xlab="Arsenic level", ylab="Average residual", type="n", main="Binned residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
lines (br.arsenic[,1], br.arsenic[,6], col="gray", lwd=.5)
lines (br.arsenic[,1], -br.arsenic[,6], col="gray", lwd=.5)
points (br.arsenic[,1], br.arsenic[,2], pch=19, cex=.5)


## ----error rate------------------------------------------------------------------------------------------
mean(data$switch)
error.rate <- mean((rep(1,length(data$switch))>0.5 & data$switch==0) | 
                     ((rep(1,length(data$switch))<0.5 & data$switch==1)))
error.rate


## ----error rate model------------------------------------------------------------------------------------
error.rate <- mean((pred.8>0.5 & data$switch==0) | (pred.8<0.5 & data$switch==1))
error.rate


## ----poiss-----------------------------------------------------------------------------------------------
y <- c(4, 27, 2, 21, 8, 3, 3, 0, 4, 4, 5, 2, 2, 5, 7, 1, 5, 0, 15, 36)
age <- c(54 ,26,41,23, 41, 75, 28, 78, 44, 52, 50, 28, 42, 50, 31, 59, 33, 35, 26, 22)
family <- c(1,1,0,1,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1)
earning <- c( 1.2, 5.3, 2.4, 6.3, 4.2, 1.5, 1.1, 0.9, 1.2, 1.1, 1.3, 1.5, 4.2, 1.3, 1.2, 1.3, 1.9, 2.7, 7.3, 5)
n <- length(y)
calls.glm <- glm(y ~ age+family+earning, family = poisson)
summary(calls.glm)



## ----poiss overd-----------------------------------------------------------------------------------------
pred.pois <- calls.glm$fitted.values
res.st <- (y-pred.pois)/sqrt(pred.pois)

plot(pred.pois, res.st)
abline(h=0,lty=3,col="gray75")

est.overdispersion <- sum(res.st^2)/16
est.overdispersion


## ----quasi poiss-----------------------------------------------------------------------------------------
calls.quasi.glm <- glm(y ~ age+family+earning, family = quasipoisson)
summary(calls.quasi.glm)
par(mfrow=c(1,2))



