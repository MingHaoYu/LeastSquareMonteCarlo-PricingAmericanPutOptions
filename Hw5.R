#title : 405-Computational Methods HW5
#author: Ming-Hao Yu
#date  : 2018-02-08
#install.packages("corpcor")
library(corpcor)

LSMC <- function(s0, k, sigma, r, t, n, paths, ktype, polynomials) { #ktype: A matrix with k row k col
    dt <- t/n
    set.seed(0)
    w1 <- sqrt(dt)*rnorm(paths/2*n)
    w <- c(w1, -w1)
    dwTable <- matrix(nrow=paths, ncol=n)
    stockPriceTable <- matrix(nrow=paths, ncol=n+1, s0)
    
    for(i in 1:(paths)) {
        dwTable[i, ] <- w[((i-1)*n+1):(i*n)]
    }
    # use simulated random walks to generate simulated stock prices
    for(i in 1:n) {
        stockPriceTable[, i+1] <- stockPriceTable[, i] + r*stockPriceTable[, i]*dt + sigma*stockPriceTable[, i]*dwTable[, i]
    }
    
    if(toupper(polynomials) == "LAGUERRE") {
        stockPriceTable = stockPriceTable/k
        k0 <- k
        k <- 1
    }
      
    cashFlow <- 0
    price <- 0
    for(i in 1:paths) { cashFlow[i] <- max(k-stockPriceTable[i, n+1], 0)}
    
    for(i in n:2) {
        A <- matrix(nrow=ktype, ncol=ktype)
        B <- matrix(nrow=ktype, ncol=1)
        EV <- 0
        for(j in 1:paths) {EV[j] <- max(k-stockPriceTable[j, i], 0)}
        
        mask <- which(EV>0)
        x <- stockPriceTable[mask, i]
        y <- exp(-r*dt)*cashFlow[mask]
        
        if(toupper(polynomials) == "LAGUERRE") {
            L1 <- exp(-x/2)
            L2 <- exp(-x/2)*(1-x)
            L3 <- exp(-x/2)*(1-2*x+x^2/2)
            L4 <- exp(-x/2)*(1-3*x+3*x^2/2-x^3/6)
        }
        else if(toupper(polynomials) == "HERMITE") {
            L1 <- rep(1,length(x))
            L2 <- 2*x
            L3 <- 4*x^2-2
            L4 <- 8*x^3-12*x
        }
        else if(toupper(polynomials) == "SIMPLE") {#simple monomials
            L1 <- rep(1,length(x))
            L2 <- x
            L3 <- x^2
            L4 <- x^3
        }
        else {
            cat("[No:", polynomials, " polynomials type.")
            return(0)
        }
        
        if(ktype==2) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[2,2] <- sum(L2*L2)
            B <- matrix(c(sum(y*L1),sum(y*L2)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2
        }
        else if(ktype==3) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[1,3] <- A[3,1] <- sum(L1*L3)
            A[2,3] <- A[3,2] <- sum(L2*L3)
            A[2,2] <- sum(L2*L2)
            A[3,3] <- sum(L3*L3)
            B <- matrix(c(sum(y*L1),sum(y*L2),sum(y*L3)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2 + a[3]*L3
        }
        else if(ktype==4) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[1,3] <- A[3,1] <- sum(L1*L3)
            A[1,4] <- A[4,1] <- sum(L1*L4)
            A[2,3] <- A[3,2] <- sum(L2*L3)
            A[2,4] <- A[4,2] <- sum(L2*L4)
            A[3,4] <- A[4,3] <- sum(L3*L4)
            A[2,2] <- sum(L2*L2)
            A[3,3] <- sum(L3*L3)
            A[4,4] <- sum(L4*L4)
            B <- matrix(c(sum(y*L1),sum(y*L2),sum(y*L3),sum(y*L4)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2 + a[3]*L3 + a[4]*L4
        }  
        else {
            cat("polynomials term ranges from 2-4")
            return(0)
        }
        
        newCashFlow <- cashFlow*exp(-r*dt)
        newCashFlow[mask] <- ifelse(EV[mask]>ECV, EV[mask], y)
        cashFlow <- newCashFlow
        price <- mean(cashFlow)
    }
    if(toupper(polynomials) == "LAGUERRE") {price = price*k0}
    return(price)
}

LSMCForwardStart <- function(s0, sigma, r, t, t2, n, paths, ktype) { #ktype: A matrix with k row k col
    dt <- t/n
    set.seed(0)
    w1 <- sqrt(dt)*rnorm(paths*n)
    w2 <- -w1[(paths/2*n+1):(paths*n)]
    dwTable <- matrix(nrow=paths, ncol=n)
    stockPriceTable <- matrix(nrow=paths, ncol=n+1, s0)
    
    for(i in 1:(paths)) {
        dwTable[i,1:(n/2)] <- w1[((i-1)*n/2+1):(i*n/2) ]
        dwTable[i,(n/2+1):n] <- w2[((i-1)*n/2+1):(i*n/2) ]
    }
    # use simulated random walks to generate simulated stock prices
    for(i in 1:n) {
        stockPriceTable[, i+1] <- stockPriceTable[, i] + r*stockPriceTable[, i]*dt + sigma*stockPriceTable[, i]*dwTable[, i]
    }
    
    cashFlow <- 0
    price <- 0
    k <- stockPriceTable[,(t2/t*n+1)]
    stockPriceTable <- stockPriceTable/k
    for(i in 1:paths) { cashFlow[i] <- max(1-stockPriceTable[i, n+1], 0)}
    
    for(i in n:(t2/t*n+2)) {
        A <- matrix(nrow=ktype, ncol=ktype)
        B <- matrix(nrow=ktype, ncol=1)
        EV <- 0
        for(j in 1:paths) {EV[j] <- max(1-stockPriceTable[j, i], 0)}
        
        mask <- which(EV>0)
        x <- stockPriceTable[mask, i]
        y <- exp(-r*dt)*cashFlow[mask]
        L1 <- rep(1,length(x))
        L2 <- x
        L3 <- x^2
        L4 <- x^3
        
        if(ktype==2) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[2,2] <- sum(L2*L2)
            B <- matrix(c(sum(y*L1),sum(y*L2)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2
        }
        else if(ktype==3) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[1,3] <- A[3,1] <- sum(L1*L3)
            A[2,3] <- A[3,2] <- sum(L2*L3)
            A[2,2] <- sum(L2*L2)
            A[3,3] <- sum(L3*L3)
            B <- matrix(c(sum(y*L1),sum(y*L2),sum(y*L3)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2 + a[3]*L3
        }
        else if(ktype==4) {
            A[1,1] <- sum(L1*L1)
            A[1,2] <- A[2,1] <- sum(L1*L2)
            A[1,3] <- A[3,1] <- sum(L1*L3)
            A[1,4] <- A[4,1] <- sum(L1*L4)
            A[2,3] <- A[3,2] <- sum(L2*L3)
            A[2,4] <- A[4,2] <- sum(L2*L4)
            A[3,4] <- A[4,3] <- sum(L3*L4)
            A[2,2] <- sum(L2*L2)
            A[3,3] <- sum(L3*L3)
            A[4,4] <- sum(L4*L4)
            B <- matrix(c(sum(y*L1),sum(y*L2),sum(y*L3),sum(y*L4)), ncol=1)
            a <- pseudoinverse(A)%*%B
            ECV <- a[1]*L1 + a[2]*L2 + a[3]*L3 + a[4]*L4
        }  
        else {
            cat("polynomials term ranges from 2-4")
            return(0)
        }
        
        newCashFlow <- cashFlow*exp(-r*dt)
        newCashFlow[mask] <- ifelse(EV[mask]>ECV, EV[mask], y)
        cashFlow <- newCashFlow
        price <- mean(cashFlow*k)
    }
    return(price*exp(-r*t2))
}


#Problem 1 
sigma <- 0.2
r <- 0.06
k <- 40
paths <- 100000
ktype <- c(2,3,4)
s0 <- c(36,40,44)
laguerre1 <- matrix(nrow=3,ncol=3)
laguerre2 <- matrix(nrow=3,ncol=3)
laguerre3 <- matrix(nrow=3,ncol=3)
hermite1  <- matrix(nrow=3,ncol=3)
hermite2 <- matrix(nrow=3,ncol=3)
hermite3  <- matrix(nrow=3,ncol=3)
simple1 <- matrix(nrow=3,ncol=3)
simple2 <- matrix(nrow=3,ncol=3)
simple3 <- matrix(nrow=3,ncol=3)

#(a) Laguerre 
for(i in 1:3) {
    for(j in 1:3) {
        laguerre1[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=0.5, n=0.5*252, paths=paths, ktype=ktype[i], "Laguerre")
        laguerre2[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=1, n=252, paths=paths, ktype=ktype[i], "Laguerre")
        laguerre3[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=2, n=2*252, paths=paths, ktype=ktype[i], "Laguerre")
    }
}
colnames(laguerre1) <- c(36, 40, 44)
rownames(laguerre1) <- c("k=2", "k=3", "k=4")
colnames(laguerre2) <- c(36, 40, 44)
rownames(laguerre2) <- c("k=2", "k=3", "k=4")
colnames(laguerre3) <- c(36, 40, 44)
rownames(laguerre3) <- c("k=2", "k=3", "k=4")

#(b) Hermite 
for(i in 1:3) {
    for(j in 1:3) {
        hermite1[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=0.5, n=0.5*252, paths=paths, ktype=ktype[i], "Hermite")
        hermite2[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=1, n=252, paths=paths, ktype=ktype[i], "Hermite")
        hermite3[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=2, n=2*252, paths=paths, ktype=ktype[i], "Hermite")
    }
}
colnames(hermite1) <- c(36, 40, 44)
rownames(hermite1) <- c("k=2", "k=3", "k=4")
colnames(hermite2) <- c(36, 40, 44)
rownames(hermite2) <- c("k=2", "k=3", "k=4")
colnames(hermite3) <- c(36, 40, 44)
rownames(hermite3) <- c("k=2", "k=3", "k=4")

#(c) Simple Monomials 
for(i in 1:3) {
    for(j in 1:3) {
        simple1[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=0.5, n=0.5*252, paths=paths, ktype=ktype[i], "simple")
        simple2[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=1, n=252, paths=paths, ktype=ktype[i], "simple")
        simple3[i,j] <- LSMC(s0=s0[j], k=k, sigma=sigma, r=r, t=2, n=2*252, paths=paths, ktype=ktype[i], "simple")
    }
}
colnames(simple1) <- c(36, 40, 44)
rownames(simple1) <- c("k=2", "k=3", "k=4")
colnames(simple2) <- c(36, 40, 44)
rownames(simple2) <- c("k=2", "k=3", "k=4")
colnames(simple3) <- c(36, 40, 44)
rownames(simple3) <- c("k=2", "k=3", "k=4")

ans1 <- list(Laguerre1 = laguerre1, Laguerre2 = laguerre2, Laguerre3 = laguerre3, 
             Hermite1 = hermite1, Hermite2 = hermite2, Hermite3 = hermite3,
              SimpleMonomial1 = simple1, SimpleMonomial2 = simple2, SimpleMonomial3 = simple3)

#Problem 2
#a
sigma <- 0.2
s0 <- 65
k <- 65
t2 <- 0.2
t <- 1
r <- 0.06

paths <- 10000
n <- 100

set.seed(0)
dw <- sqrt(t/n)*rnorm(n*paths)
dwTable <- matrix(nrow=paths, ncol=n)
stockPriceTable <- matrix(nrow=paths, ncol=n+1, s0)

for(i in 1:(paths)) {
    dwTable[i,] <- cumsum(dw[((i-1)*n+1):(i*n)])
    #dwTable[i,] <- dw[((i-1)*n+1):(i*n)]
}

# use simulated random walks to generate simulated stock prices
for(i in 1:n) {
    stockPriceTable[, i+1] <- s0*exp((r-sigma^2/2)*(t*i/n)+sigma*dwTable[, i])
    #stockPriceTable[, i+1] <- stockPriceTable[, i] + r*stockPriceTable[, i]*t/n + sigma*stockPriceTable[, i]*dwTable[, i]
}

price <- 0
for(i in 1:paths) {
    price[i] <- exp(-r*t)*max(stockPriceTable[i,(t2/t*n+1)] - stockPriceTable[i, n+1], 0)
}
EuroPutPrice <- mean(price)

#b
AmerPutPrice <- LSMCForwardStart(s0=s0, sigma=sigma, r=r, t=t, t2=t2, paths= paths, ktype=3, n=n)

ans2 <- list(EuropeanPutPrice = EuroPutPrice, AmericanPutPrice = AmerPutPrice)


ans1
ans2
