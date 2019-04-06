library(stats)
library(stats4)

x <- 1:4
x[1]<-1
x[2:4]<-c(0,0,0)
fft(x)
fft(fft(x), inverse = TRUE)/length(x)

## Slow Discrete Fourier Transform (DFT) - e.g., for checking the formula
fft0 <- function(z, inverse=FALSE) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  ff <- (if(inverse) 1 else -1) * 2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(ff*(h-1))), complex(1))
}

relD <- function(x,y) 2* abs(x - y) / abs(x + y)
n <- 2^8
z <- complex(n, rnorm(n), rnorm(n))
## relative differences in the order of 4*10^{-14} :
summary(relD(fft(z), fft0(z)))
summary(relD(fft(z, inverse=TRUE), fft0(z, inverse=TRUE)))


# FFT initial testing -----------------------------------------------------
a <- c(1, 0, 0, 0)

a.hat <- fft(a,inverse=FALSE)
a.hat
# either dft=sqrt(n)*F  or  dft=sqrt(n)F-1

b<-c(0,1,0,0)
b.hat<-fft(fft(b,inverse=TRUE)/length(b))
b.hat
#--------------------------------------------------------------------------------


# Binary Search Using library(data.tables) --------------------------------
set.seed(2L)
N = 2e7L
DT = data.table(x = sample(letters, N, TRUE),
                y = sample(1000L, N, TRUE),
                val = runif(N), key = c("x", "y"))
print(object.size(DT), units = "Mb")
key(DT)

## (1) Usual way of subsetting - vector scan approach
t1 <- system.time(ans1 <- DT[x == "g" & y == 877L])
t1
head(ans1)
dim(ans1)

## (2) Subsetting using keys
t2<-system.time(ans2<-DT[.("g", 877L)])
t2
head(ans2)
dim(ans2)

## Test if == to each other
identical(ans1$val,ans2$val)
