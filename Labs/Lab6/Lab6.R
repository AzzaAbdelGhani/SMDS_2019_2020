## ----global_options, include=FALSE-----------------------------------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10))#, fig.path = 'figs/')
library(MASS)

## ----setup, include=FALSE--------------------------------------------------------------------------
library(knitr)
local({
  hook_plot = knit_hooks$get('plot')
  knit_hooks$set(plot = function(x, options) {
    paste0('\n\n----\n\n', hook_plot(x, options))
  })
})


## ----pairs, eval=TRUE, echo=FALSE------------------------------------------------------------------
library(mgcv)
panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
pairs(trees, panel = panel.smooth, main = "trees data",diag.panel = panel.hist)


## ----glm.1 vs gam.1--------------------------------------------------------------------------------
glm.1 <- glm(Volume ~ Height, family = Gamma(link=log), data=trees)
gam.1 <- gam(Volume ~ s(Height), family=Gamma(link=log), data=trees)

summary(gam.1)

plot(trees$Height, trees$Volume, xlab="Height", ylab="Fitted values")
points(trees$Height, glm.1$fitted.values, col="blue", bg=3, pch=19)
points(trees$Height, gam.1$fitted.values, col="red", bg=2, pch=19, cex=.5)
legend("topright", c("GLM","GAM"), pch=19, col=c("blue", "red"))


## ----gam.1 residuals-------------------------------------------------------------------------------
plot(gam.1,residuals=TRUE,pch=19)


## ----second----------------------------------------------------------------------------------------
glm.2<-glm(Volume ~ Girth + Height, family = Gamma(link=log), data=trees)
gam.2 <- gam(Volume ~ s(Height) + s(Girth), family=Gamma(link=log), data=trees)
summary(gam.2)


## ----plot, echo=FALSE------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(trees$Height, trees$Volume, xlab="Height", ylab="Fitted values")
points(trees$Height, glm.2$fitted.values, col="blue", pch=22, bg="blue")
points(trees$Height, gam.2$fitted.values, col="red", pch=23, bg="red")
legend("topright", c("GLM","GAM"), pch=c(22,23), col=c("blue", "red"))

plot(trees$Girth, trees$Volume, xlab="Girth", ylab="Fitted values")
points(trees$Girth, glm.2$fitted.values, col="blue", pch=22, bg="blue")
points(trees$Girth, gam.2$fitted.values, col="red", pch=23, bg="red")
legend("topright", c("GLM","GAM"), pch=c(22,23), col=c("blue", "red"))


par(mfrow=c(1,2))
plot(gam.2, residuals =TRUE, pch =19)


## ----vis.gam, warning=FALSE------------------------------------------------------------------------
par(mfrow=c(1,2))
vis.gam(gam.2,theta=-45,type = "response", color="terrain")
vis.gam(gam.2,theta=-45,type = "link", color="terrain")
par(mfrow=c(1,2))
vis.gam(gam.2,type = "response", color="terrain", plot.type = "contour")
vis.gam(gam.2,type = "link", color="terrain", plot.type = "contour")



## ----aic-------------------------------------------------------------------------------------------
AIC(glm.1, gam.1, glm.2, gam.2)


## ----gam.2 playing with basis dimension and penalization-------------------------------------------
#unpenalized regression spline
gam.2.1 <- gam(Volume ~ s(Height, k=10, fx=TRUE)+
                 s(Girth, k=10, fx=TRUE),
               family=Gamma(link=log), data=trees)
gam.2.1$sp
summary(gam.2.1)


## ----gam.2 playing with basis dimension and penalization-graph, warning=FALSE----------------------
par(mfrow=c(1,2))
plot(gam.2.1,residuals=TRUE,pch=19)
par(mfrow=c(1,1))
vis.gam(gam.2.1,theta=-45,type = "link", color="terrain")



## ----lambda_aic------------------------------------------------------------------------------------
#extract the values
sp <- gam.2$sp

tuning.scale<-c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5) 
scale.exponent <- log10(tuning.scale) 
n.tuning <- length(tuning.scale) 

minus2loglik <- rep(NA,n.tuning) 
edf <- rep(NA,n.tuning)  
aic <- rep(NA,n.tuning)  


for (i in 1:n.tuning) {   
  gamobj <- gam(Volume ~ s(Height) + s(Girth), family=Gamma(link=log),
                data=trees, sp=tuning.scale[i]*sp) 
  minus2loglik[i] <- -2*logLik(gamobj) 
  edf[i]          <- sum(gamobj$edf)+1  
  aic[i]          <- AIC(gamobj)
}

par(mfrow=c(2,2)) 
plot(scale.exponent,minus2loglik,type="b",main="-2 log likelihood") 
plot(scale.exponent,edf,ylim=c(0,70),type="b",main="effective number of parameters")  
plot(scale.exponent,aic,type="b",main="AIC build-in function") 
plot(scale.exponent,minus2loglik+2*edf,type="b",main="AIC") 

# find the minimum
opt.tuning.scale <- tuning.scale[which.min(aic)]
opt.tuning.scale
opt.sp<-opt.tuning.scale*sp

# fitting the final model with the optimal level of smoothing

gam.2.opt <- gam(Volume ~ s(Height)+s(Girth), family=Gamma(link=log),data=trees,
                 sp=opt.sp)

AIC(gam.2, gam.2.opt)



## ----tree load data--------------------------------------------------------------------------------
library(reshape2)
library(rpart)
library(rpart.plot)
spam <- read.table("spambase.data", sep=",")
spam <- spam[,c(58,#"yesno"
                57,#"crl.tot"
                53,#"dollar" 
                52,#"bang"
                24,#"money" 
                23,#"n000"
                1  #"make"
                )]
colnames(spam) <- c("yesno","crl.tot","dollar","bang", "money","n000", "make")
spam$yesno <- factor(spam$yesno, levels = c(0, 1))
levels(spam$yesno) <- c("n","y")

par(mfrow=c(2,3))
for(i in 2:7){
  boxplot(spam[,i]~yesno, data=spam, ylab = colnames(spam)[i])
}


## ----logistic fails--------------------------------------------------------------------------------
#summary(glm(yesno ~ crl.tot + dollar + bang + money + n000 + log(make+.5),
#            family=binomial, data=spam))


## ----tree------------------------------------------------------------------------------------------
library(rpart)

spam.rpart <- rpart(yesno ~ crl.tot + dollar + bang + money + n000 + make,
                    method="class", data=spam) 
plot(spam.rpart) # Draw tree
text(spam.rpart) # Add labeling


## --------------------------------------------------------------------------------------------------
printcp(spam.rpart)


## --------------------------------------------------------------------------------------------------
spam.rpart0001 <- rpart(yesno ~ crl.tot + dollar + bang + money + n000 + make,
                    method="class", data=spam, cp = 0) 
printcp(spam.rpart0001)


## --------------------------------------------------------------------------------------------------
plotcp(spam.rpart0001)


## --------------------------------------------------------------------------------------------------

best.cp <- spam.rpart0001$cptable[which.min(spam.rpart0001$cptable[,"xerror"]),]
best.cp

sd.rule <- best.cp["xerror"]+best.cp["xstd"]
cptable.sd.rule <- spam.rpart0001$cptable[spam.rpart0001$cptable[,"xerror"]<=sd.rule,]
best.cp.sd <- cptable.sd.rule[which.min(cptable.sd.rule[,"nsplit"]),]
best.cp.sd

tree.pruned <- prune(spam.rpart0001, cp=best.cp[1])
rpart.plot(tree.pruned, extra=104, box.palette="GnBu",branch.lty=3, shadow.col="gray", nn=TRUE)

tree.pruned.sd <- prune(spam.rpart0001, cp=best.cp.sd[1])
printcp(tree.pruned.sd)
rpart.plot(tree.pruned.sd, extra=104, box.palette="GnBu",branch.lty=3, shadow.col="gray", nn=TRUE)



## --------------------------------------------------------------------------------------------------
library(randomForest)
spam.rf <- randomForest(yesno ~ ., data=spam, importance=TRUE) 
print(spam.rf)
importance(spam.rf)

