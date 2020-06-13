## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10), fig.path = 'figs/')
library(MASS)

## ----setup, include=FALSE------------------------------------------------
library(knitr)
local({
  hook_plot = knit_hooks$get('plot')
  knit_hooks$set(plot = function(x, options) {
    paste0('\n\n----\n\n', hook_plot(x, options))
  })
})

## ----basic, echo=TRUE----------------------------------------------------
c(2,5,6,4,5)
c("Roma", "Milano", "Trieste")
c(T,F,T,T,F)

#standard way for storing a vector, with <- assignment
x <- c( 0.2, 0.4, 0.5, 0.6)


## ----basic 2, echo=TRUE--------------------------------------------------
# the vector length
length(x)
# test the condition of"greater than"
x > 0.5
# test the condition "not"
x != 0.2

## ----basic 3, echo=TRUE--------------------------------------------------
# extract the first element of the vector
x[1]
# extract the first two items of the vector
x[c(1,2)]
# omit the first two items of the vector
x[-c(1,2)]
# extract given a condition
x[ x > 0.5 ]
# concatenate vectors
y <- c(1,3,5,7)
z <- c(x,y)
z
w<-(x+y)^2
w


## ----basic 4, echo=TRUE--------------------------------------------------
seq(from=4, to=8, by=0.5)
seq(from=4, to=8, length.out=9)
rep(c("Bello", "Brutto"), 4)

c(rep("Data scientists",1), rep("are creasy",1))



## ----basic 5, echo=TRUE--------------------------------------------------
# define a matrix of 1 values with 3 rows and 5 columns
mat<-matrix(1, nrow=3, ncol=5)
mat
# define a matrix where each row cell represents the row number
mat2<-matrix( c(rep(1,5), rep(2,5), rep(3,5)),nrow=3,ncol=5, byrow=TRUE)
mat2
# concatenate matrices
mat_conc<-rbind(mat, mat2)
# dimension of the matrix
dim(mat_conc)
mat_conc
# extract the first columns of the matrix
mat_conc[,1]
# extract the first row of the matrix
mat_conc[1,]

## ----basic 6, echo=TRUE--------------------------------------------------
# compute the total sum of each row
apply(mat_conc,1, sum)
# compute the mean of each row
apply(mat_conc,1, mean)
# compute the total sum of each column
apply(mat_conc,2, sum)


## ----basic 7, echo=TRUE--------------------------------------------------
mat%*%t(mat2)


## ----basic 8, echo=TRUE--------------------------------------------------
year<-c(1800, 1850, 1900, 1950, 2000)
carbon<-c(8, 54, 534, 1630, 6611)
plot(carbon~year, pch=16)
plot(year, carbon, pch=3, col="red")
fossil_frame<-data.frame(year=year, carbon=carbon)
fossil_frame

## ----basic 9, echo=TRUE--------------------------------------------------
# with this function I access general information about the dataset
summary(cars)
is.data.frame(cars)
# display the first 5 items
head(cars, n=5)
# extract column names
colnames(cars)

## ----basic 10, echo=TRUE-------------------------------------------------
# extract the first five rows
cars[1:5,]
# extract the first five rows and the first column
cars[1:5,]

## ----basic 11, echo=TRUE-------------------------------------------------
# extract the units with a speed greater than 50 mph
subset(cars, speed>20)
# extract the units with the dist minor than 10 or greater than 25
subset(cars, dist<10 | dist>75)

## ----data frame,echo=TRUE------------------------------------------------
# for .txt files

studenti <- read.table( file="file_text.txt", header=TRUE, sep="", dec=".")
studenti

# for .csv files
teams <- read.csv(file="file_excel.csv", header=TRUE, sep=";")
teams



## ----data frame 2,echo=TRUE----------------------------------------------
# for .txt files
studenti$voto
teams$shots.home


## ----basic 12, echo=TRUE-------------------------------------------------
s <- c(20,30,45,30, 34)

if (mean(s)<=20){
  print("mean is minor than 20")
}else if(mean(s)>20 & mean(s)>=30){
  print("mean is bounded between 20 and 30")
}else{
  print("mean is greater than 30")
}

## ----basic 13, echo=TRUE-------------------------------------------------
n<-20
fibonacci_seq<-c()
fibonacci_seq[1]<-1
fibonacci_seq[2]<-1
for (i in 3: n){
  fibonacci_seq[i] <- fibonacci_seq[i-2]+fibonacci_seq[i-1]
}
fibonacci_seq


## ----basic 14, echo=TRUE-------------------------------------------------
# generic function syntax
foo <- function(param1, param2){
  output=f(param1, param2)
return(output)
}


## ----basic 15, echo=TRUE-------------------------------------------------
parabola <- function(x,a,b,k){
 y=a*x^{2}+b*x+k
return(y)
}
parabola(2,-2,3,1)
curve(parabola(x, a=1, b=1, k=1), xlim=c(-8,8), ylab="y")
curve(parabola(x, a=2, b=2, k=1), xlim=c(-8,8),col="red", lty=2, ylab="y", add=TRUE)


## ----examples, echo=TRUE-------------------------------------------------
#Compute the density of a normal(0,1) in 1.96
dnorm(1.96)


## ----examples_2, echo=TRUE-----------------------------------------------
set.seed(1)
#Compute the density of a normal(1,2) in 1.96
dnorm(1.96, mean=1, sd=sqrt(2))
#Compute the distribution function of a gamma(2,2) in 1
pgamma(1,shape=2, rate=2, lower.tail=TRUE)
#generate 100 values from a Poisson distribution with
#rate 2 and plot the correspondent histogram
z<-rpois(100,2)
barplot(table(z), ylab="frequency")




## ----binom, echo=TRUE----------------------------------------------------
#graphical setting for margins and type of points
par(mfrow=c(1,2),mar=c(5,4,4,1), oma=c(0,0.2,0.2,0), pty="s", pch = 16)
#plot the binomial distributions with different input
plot(0:20, dbinom(0:20, 20, 0.2), 
     xlab = "x", ylab = "f(x)", cex.lab=2, main="n=20", cex.main=2)
plot(0:50, dbinom(0:50, 50, 0.5), xlab ="x", ylab = "f(x)",
      cex.lab=2, main= "n=50", cex.main=2)

## ----geom,  echo=TRUE----------------------------------------------------
X=0:10
# use dgeom(x,p), with x number of failures and probability of success p 
plot(X, dgeom(X,0.6), type="h", ylim=c(0,1), main="Geometric distribution for p=0.6", ylab="f(x)", xlab="X=Number of failures before first success")

## ----binom_neg,  echo=TRUE-----------------------------------------------
k<-3;r<-20; p<-0.08; 
dnbinom(x=r, size=k,prob=p )
plot(0:1000, dnbinom(0:1000, size=r, prob=p), xlab="k",
  ylab="f(x)",pch=21, bg=1)

## ----gamma, fig.height=2, fig.wifth=3, echo=TRUE-------------------------
n<-1000; alpha<-2; beta<-2; sample_rep<-1000
X<-matrix(NA, n, sample_rep)
for (h in 1:n){
  X[h,]<-rexp(sample_rep, beta)
}
Y<-apply(X,1,sum)
hist(Y, breaks=40, probability=TRUE)
curve(dgamma(x, n, beta), col="red", lwd=2, add=TRUE)


## ----mix, fig.height=5, fig.width=6, echo=TRUE---------------------------
t_mixture <- function(mean,  df, n){
  
  Z=rgamma(n, df/2, df/2)
  X=rnorm(n, mean, sd=sqrt(1/sqrt(Z)))
  return(X)
  
}
df<-5
hist(t_mixture(0,df,10000), probability=TRUE, breaks=40, 
  main=paste("Histogram for a t with", df, "d.f"), xlab="x")
curve(  dt(x,5), col="red", lwd=2, add=TRUE )

## ----ecdf, echo=TRUE-----------------------------------------------------
set.seed(2)
par(mfrow=c(1,2))
n<-50
y<-rbeta(n, 3,4)
edf_beta<-ecdf(y)
tt<-seq(from=0, to=1, by=0.01)
plot(edf_beta, verticals=TRUE, do.p=FALSE, main="ECDF and CDF: n=50")
lines(tt, pbeta(tt,3,4), col=2, lty=2, lwd=2)
n2<-500
y2<-rbeta(n2, 3,4)
edf_beta2<-ecdf(y2)
tt<-seq(from=0, to=1, by=0.01)
plot(edf_beta2, verticals=TRUE, do.p=FALSE, main="ECDF and CDF: n=500")
lines(tt, pbeta(tt,3,4), col=2, lty=2, lwd=2)



## ----qq, echo=TRUE-------------------------------------------------------
par(mfrow=c(1,2))
qqplot(qbeta(ppoints(n),3,4),y,
  xlab="True quantiles", ylab="Sample quantiles",
  main = "Q-Q plot for Beta(3,4): n=50")
qqline(y, distribution = function(p) qbeta(p, 3,4),
       prob = c(0, 1), col = 2)

qqplot(qbeta(ppoints(n2),3,4),y2,
  xlab="True quantiles", ylab="Sample quantiles",
  main = "Q-Q plot for Beta(3,4): n=500")
qqline(y2, distribution = function(p) qbeta(p, 3,4),
       prob = c(0, 1), col = 2)


## ----binom_appr_10, echo=TRUE--------------------------------------------
#set the seed
set.seed(123)
par(mfrow=c(1,2), mar=c(5,4,4,1), oma=c(1,1,1,1))
p<-0.5; size1<-10; size2<-50
#normal approximation
prob_bin <- dbinom(seq(0:size1), size1,p)
curve(dnorm(x, size1*p, sqrt(size1*p*(1-p))), xlim=c(0,10), xlab="x", ylab="f(x)",
  main=paste("n=", size1), cex.main=1.5)
lines(prob_bin, type="h", main="n=10", xlim=c(1,10), 
     cex.lab=1.5, col="red")

prob_bin <- dbinom(seq(0:size2), size2,p)
curve(dnorm(x, size2*p, sqrt(size2*p*(1-p))), xlim=c(0,size2), xlab="x", ylab="f(x)", main=paste("n=", size2), cex.main=1.5)
lines(prob_bin, type="h", main="n=10", xlim=c(1,size2), 
     cex.lab=1.5, col="red")


## ----binom_water, echo=TRUE,  results='hold', fig.keep='high', fig.height=8, fig.width=5----
p<-0.5; q<-0.7; n<-20; m<-20
Prob_posillipo=pnorm(0, mean=p*n-q*m, 
               sd=sqrt( n*p*(1-p)+m*q*(1-q)  ),  
               lower.tail=FALSE)
Prob_posillipo

## ----binom_water2, echo=TRUE,  results='hold', fig.keep='high', fig.height=8, fig.width=5----
conf_int<-function(alpha,n,p,m,q){
  mean=p*n-q*m; sd=sqrt( n*p*(1-p)+m*q*(1-q));
  return(c(qnorm(alpha/2, mean=mean, sd=sd  ),
           qnorm( 1-alpha/2, mean=mean, sd=sd  ) ))
}

conf_int(0.05, n,p,m,q)
conf_int(0.5, n,p, m,q)

## ----binom_water3, echo=TRUE---------------------------------------------
par(mar=c(5,4,4,1))
curve(dnorm(x,p*n-q*m, sqrt( n*p*(1-p)+m*q*(1-q)) ),
      xlim=c(-12,20), ylim=c(0,0.3), ylab="f(x)", cex.lab=2)
segments(conf_int(0.05,n,p,m,q)[1],-0.01,
        conf_int(0.05,n,p,m,q)[1],
        dnorm(conf_int(0.05,n,p,m,q)[1],p*n-q*m,
              sqrt( n*p*(1-p)+m*q*(1-q)) ),
        col="black",lty=4, lwd=2)
segments(conf_int(0.05,n,p,m,q)[2],-0.01,
        conf_int(0.05,n,p,m,q)[2],
        dnorm(conf_int(0.05,n,p,m,q)[2],p*n-q*m,
              sqrt( n*p*(1-p)+m*q*(1-q)) ), 
        col="black",lty=4, lwd=2)
segments(0, -0.01, 0, dnorm(0, p*n-q*m,
              sqrt( n*p*(1-p)+m*q*(1-q)) ), lwd=2 )
points(dbinom(0:n, n, p), pch=21, bg=1)
points(dbinom(0:m, m, q), pch=21, bg=2)
text(15, 0.1, "Y", cex=2, col=2)
text(10,0.1, "X", cex=2, col=1)
text(-5, 0.15, "X-Y", cex=2, col=1)



## ----binom_law2, echo=TRUE,  results='hold', fig.keep='high', fig.height=8, fig.width=5----
par(mar = c(4, 4, 0.5, 0.5),mfrow=c(1,1))
#I write the function depending on two arguments: n and p
law_large_numb<-function(n,p){
  x<-rbinom(n, 1,p)
  return(mean(x))
}
law_large_numb=Vectorize(law_large_numb)

## ----binom_law3, echo=TRUE,  results='hold', fig.keep='high', fig.height=4, fig.width=4----
set.seed(12345);par(mar = c(6.2, 4, 0, 0.5),mfrow=c(1,1))
curve(law_large_numb(x, p=0.5), 1,1000, xlab="n", 
      ylab="frequency")
abline(a=0.5, b=0, lwd=2, col="red")
legend(80, 0.9, lty=1, col=c("black", "red"), 
       lwd=c(1,2), c( "X/n", "p"))

