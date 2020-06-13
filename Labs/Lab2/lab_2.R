set.seed(123)
R <- 1000
n <- 30
#generate n values R times
samples <- array(0, c(3,R,n))
for (i in 1:R){
  samples[1, i, ] <- rnorm(n, 0, 1)
  samples[2, i, ] <- rt(n, df = 3)
  samples[3, i, ] <- runif(n, 0, 1)
}
#compute the sample statistics

samples_stat <- array(0, c(3, 2, R))
for (j in 1:3){
  #sample mean over R replications 
  samples_stat[j, 1, ] <- apply(samples[j, , ], 1 , mean )
  #sample variance over R replications
  samples_stat[j, 2, ] <- apply(samples[j, , ], 1 , var )
}

#visualize the results
par (mfrow=c(1,3), oma=c(0,0,0,0))
hist(samples_stat[1, 1, ], breaks= 40, probability = TRUE, 
     xlab="y", main= "N(0,1)", cex.main=1.5)
#overlap the true distribution for the sample mean
curve(dnorm(x, 0, sqrt(1/n)), add = TRUE, col = "red", lwd = 2)

hist(samples_stat[2, 1, ], breaks= 40, probability = TRUE, 
     xlab="y", main= expression(t[3]), cex.main =1.5)
#overlap the asymptotic distribution for the sample mean
curve(dnorm(x, 0, sqrt(  (3/((3-2)*n)))), add = TRUE, col = "red", lwd = 2)

hist(samples_stat[3, 1, ], breaks= 40, probability = TRUE, 
     xlab="y", main=  "U(0,1)", cex.main = 1.5)
#overlap the asymptotic distribution for the sample mean
curve(dnorm(x, 1/2, sqrt(1/(12*n))), add = TRUE, col = "red", lwd = 2)

######################################################################################################

#visualize the results: variance
sigma <- 1
par (mfrow=c(1,1), oma=c(0,0,0,0))
hist(samples_stat[1, 2, ], breaks= 40, probability = TRUE, 
     xlab=expression(s^2), main= bquote(s^2), cex.main=1.5)
curve(((n-1)/sigma^2) * dchisq(x * ((n-1)/sigma^2), df = n - 1),
      add = TRUE, col="red", lwd=2, main="N(0,1)")

######################################################################################################


set.seed(4)
# number of replications
R <- 1000
#number of samples for each replication
n <- 10
mu <- 5
sigma <- 2

estimator_1 <- c()
estimator_2 <- c()
estimator_3 <- c()
estimator_4 <- c()
y <- matrix(NA, R,n)

for (i in 1:R){
  y[i, ] <- rnorm(n, mu, sigma)
  estimator_1[i] <- mean(y[i, ])
  estimator_2[i] <- median(y[i, ])
  estimator_3[i] <- (max(y[i, ])+min(y[i, ]))/2
  estimator_4[i] <- sum(sort(y[i, -c(1,n)]))/(n-2)
}

par(mfrow=c(1,1), xaxt="n")
boxplot(cbind(estimator_1, estimator_2, estimator_3, estimator_4), main="Comparison between four estimators")
par(xaxt="s")
axis(1, 1:4, c(expression(hat(mu)[1]), expression(hat(mu)[2]), expression(hat(mu)[3]), expression(hat(mu)[4])) )
abline(h=5, lwd=2, col="blue")


######################################################################################################

estimators <- cbind(estimator_1, estimator_2, estimator_3, estimator_4)
variances <- c()
estimators
for (g in 1:4){
  variances[g] <- var (estimators[,g])
}

variances


######################################################################################################

#normal case

estimators_cons <- matrix(NA, R, 4)
n <- 200


for (i in 1:R){
  y <- rnorm(n, mu, sigma)
  estimators_cons[i , 1] <- mean(y)
  estimators_cons[i , 2] <- median(y)
  estimators_cons[i , 3] <- (max(y)+min(y))/2
  estimators_cons[i , 4] <- sum(sort(y[ -c(1,n)]))/(n-2)
}


# n =10
par(mfrow=c(2,4))
for (g in 1:4){
  hist(estimators[,g],  probability = TRUE, 
       breaks=40, main=substitute(hat(mu)[g] ,list(g = g)),
       xlab="", xlim=c(0,10), cex.main = 1.5)
  abline(v=mu, col="blue", lwd=2)
}

# n= 200

for (g in 1:4){
  hist(estimators_cons[,g],  probability = TRUE, 
       breaks=40, main=substitute(hat(mu)[g] ,list(g = g)), 
       xlab="", xlim=c(0,10), cex.main = 1.5)
  abline(v=mu, col="blue", lwd=2)
}



######################################################################################################
#normal case
set.seed(345)

mu <- 5
R <- 100
n <- 10

plot(0,0,xlim=c(0,10),ylim=c(0,11), type="n", xlab=expression(mu), ylab="",
     main = paste("100 IC for the mean (known variance)"), cex.main=1.2)
abline(v=mu)
alpha <- 0.05
inf <- vector(mode="numeric", length=R)
sup <- vector(mode="numeric", length=R)
l <- vector(mode="numeric",   length=R)
d <- 0
y <- matrix(NA, R,n)
sigma <- 2
for (i in 1:R)
{
  y <- rnorm(n, mu, sigma)
  inf[i] <- mean(y)-qnorm(1-alpha/2)*sigma/sqrt(n)
  sup[i] <- mean(y)+qnorm(1-alpha/2)*sigma/sqrt(n)
  d <- d + 0.1 
  l[i] <- (mu > inf[i] & mu < sup[i]) 
  lines(seq(inf[i],sup[i],length=100), rep(d, 100), col=(l[i]+1))
}

######################################################################################################

sum(l)

######################################################################################################
install.packages("kableExtra")
library(kableExtra)
library(DAAG)
options(knitr.table.format = "html") 
pair_data_frame <- cbind(pair65, pair65[,1]-pair65[,2])
pair_data_frame <- cbind(c(1:9), pair_data_frame)
dimnames(pair_data_frame)<- list(1:9, c("pair","heated", "ambient", "difference"))

kable(pair_data_frame, "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


######################################################################################################

#compute quantiles for student t with n-1 =8 degrees of freedom
alpha <- 0.05
q_inf_95 <- qt(alpha/2, df=8)
q_sup_95 <- qt(1-alpha/2, df=8)

alpha <- 0.01
q_inf_99 <- qt(alpha/2, df=8)
q_sup_99 <- qt(1-alpha/2, df=8)


#compute the confidence intervals

library(RColorBrewer)
plotclr <- brewer.pal(6,"YlOrRd")
d <- mean(pair_data_frame[,4])
s <- sd(pair_data_frame[,4] )
n <- length(pair_data_frame[,4])

conf_int_95 <- c(d-q_sup_95*s/sqrt(n),d-q_sup_95*s/sqrt(n))
conf_int_99 <- c(d-q_sup_99*s/sqrt(n),d-q_sup_99*s/sqrt(n))

par(mfrow=c(1,2), oma=c(0,0,0,0))

# 95% confidence interval

curve(dt(x,8),xlim=c(-6,6), ylim=c(0,0.4),
      main="", col = "blue", lwd = 2, xlab="y", ylab=expression(t[8]),  yaxs="i")
cord.x <- c(q_inf_95,seq(q_inf_95,q_sup_95,0.01),q_sup_95)
cord.y <- c(0,dt(seq(q_inf_95,q_sup_95,0.01),8),0)
polygon(cord.x,cord.y,col=plotclr[3], border = NA )
curve(dt(x,8),xlim=c(-6,6),main=expression(t[8]), col = "blue", lwd = 2, add = TRUE, yaxs="i")

# 99% confidence interval

curve(dt(x,8),xlim=c(-6,6), ylim=c(0,0.4),
      main="", col = "blue", lwd = 2, xlab="y",  ylab=expression(t[8]),  yaxs="i")
cord.x2 <- c(q_inf_99,seq(q_inf_99,q_sup_99,0.01),q_sup_99)
cord.y2 <- c(0,dt(seq(q_inf_99,q_sup_99,0.01),8),0)
polygon(cord.x2,cord.y2,col=plotclr[5], border = NA )
curve(dt(x,8),xlim=c(-6,6),main=expression(t[8]), col = "blue", lwd = 2, add = TRUE, yaxs="i")


######################################################################################################

Musicians <- c( "Selena Gomez",  "Ariana Grande", "Beyonce" , "Taylor Swift",  "Justin Bieber", "Nicki Minaj")
Others <- c("Cristiano Ronaldo", "Kim Kardashian",  "Kylie Jenner", "Dwayne Johnson", "Neymar", "Lionel Messi", 
            "Kendall Jenner", "Kourtney Kardashian", "Kevin Hart")

n1 <- length(Musicians)
n2 <- length(Others)

Followers_M<- c(135, 118, 113, 107,  98,  86 )
Followers_O <- c(123, 110, 106, 103,  91,  89,  89,  62,  58)


######################################################################################################

t <- t.test(Followers_M, Followers_O, mu = 0, alternative ="greater", var.equal =TRUE) 
t

######################################################################################################

curve(dt(x,13),xlim=c(-5,5), ylim=c(0,0.4),
      main="p-values and rejection region", col = "blue", lwd = 2, xlab="x-y",  ylab=expression(t[13]),  yaxs="i")
cord.x <- c(qt(0.95,13),seq(qt(0.95,13), 5, 0.01), 5)
cord.y <- c(0,dt(seq(qt(0.95,13), 5, 0.01),13),0)
polygon(cord.x,cord.y,col=plotclr[3], border = NA )
curve(dt(x,13),xlim=c(-5,5),main=expression(t[13]), col = "blue", lwd = 2, add = TRUE, yaxs="i")
abline(v =t$statistic, lty=2, lwd=2, col="red")
text (0,0.2, paste("Accept", expression(H0)))
text (2.7,0.08, paste("Reject", expression(H0)))
text(as.double(t$statistic)-0.15, 0.02, "t", col="red", cex=1.2)

######################################################################################################
set.seed(101)
n <- 50
K <- 4
# generate the values
y <- sample( 1:K, n, replace=TRUE, prob =c( 7/16, 5/16, 3/16, 1/16))
y
observed <- table(y)
observed
expected <- c( n*(7/16), n*(5/16), n*(3/16), n*(1/16))
expected
x2 <- sum( (observed-expected)^(2)/expected)
x2
#manually compute the p-value
pchisq(x2, df =K-1, lower.tail =FALSE )

######################################################################################################
#same result with the chisq.test function
chisq.test( observed, p = c( 7/16, 5/16, 3/16, 1/16))

######################################################################################################

y2 <- sample( 1:K, n, replace=TRUE, prob=c(5/16, 3/16, 6/16, 2/16))
new_observed <- table(y2)
#new test after the energetic drink
chisq.test( new_observed, p = c( 7/16, 5/16, 3/16, 1/16)  )

######################################################################################################
table(y)
table(y2)

######################################################################################################
Owners <- c( "Katy Perry", "Justin Bieber", "Taylor Swift", "Cristiano Ronaldo",
             "Kim Kardashian", "Ariana Grande", "Selena Gomez", "Demi Lovato")
Instagram <- c( 69, 98,107, 123, 110, 118, 135, 67)
Twitter <- c( 109, 106, 86, 72, 59, 57, 56, 56)
plot( Instagram, Twitter, pch=21, bg=2, xlim=c(60, 150), ylim=c(40, 120) )
text( Instagram[-6], Twitter[-6]+5, Owners[-6], cex=0.8 )
text( Instagram[6], Twitter[6]-5, Owners[6], cex=0.8 )

######################################################################################################

log_lik_weibull <- function( data, param){
  -sum(dweibull(data, shape = param[1], scale = param[2], log = TRUE))
}
######################################################################################################
y <- c(155.9, 200.2, 143.8, 150.1,152.1, 142.2, 147, 146, 146,
       170.3, 148, 140, 118, 144, 97)
n <- length(y)

#define parameters grid
gamma <- seq(0.1, 15, length=100)
beta <- seq(100,200, length=100)
parvalues <- expand.grid(gamma,beta)
llikvalues <- apply(parvalues, 1, log_lik_weibull, data=y)
llikvalues <- matrix(-llikvalues, nrow=length(gamma), ncol=length(beta),
                     byrow=F)
conf.levels <- c(0,0.5,0.75,0.9,0.95,0.99)

#contour plot
contour(gamma, beta, llikvalues-max(llikvalues),
        levels=-qchisq(conf.levels, 2)/2,
        xlab=expression(gamma),
        labels=as.character(conf.levels),
        ylab=expression(beta))
title('Weibull relative log likelihood')


#image
image(gamma,beta,llikvalues-max(llikvalues),zlim=c(-6,0),
      col=terrain.colors(20),xlab=expression(gamma),
      ylab=expression(beta))
title('Weibull relative log likelihood')


######################################################################################################


gammahat<-uniroot(function(x) n/x+sum(log(y))-n*
                    sum(y^x*log(y))/sum(y^x),
                  c(1e-5,15))$root
betahat<- mean(y^gammahat)^(1/gammahat)
weib.y.mle<-c(gammahat,betahat)
#first element is the MLE for the shape gamma, second element the MLE for the scale beta
weib.y.mle

######################################################################################################

#observed information matrix
jhat<-matrix(NA,nrow=2,ncol=2)
jhat[1,1]<-n/gammahat^2+sum((y/betahat)^gammahat*
                              (log(y/betahat))^2)
jhat[1,2]<-jhat[2,1]<- n/betahat-sum(y^gammahat/betahat^(gammahat+1)*
                                       (gammahat*log(y/betahat)+1))
jhat[2,2]<- -n*gammahat/betahat^2+gammahat*(gammahat+1)/
  betahat^(gammahat+2)*sum(y^gammahat)
solve(jhat)


#se of the mle
mle.se<-sqrt(diag(solve(jhat)))
mle.se
