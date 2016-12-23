
x <- 1:100
y1 <- c()
for(j in 1:100){ y1 <- cbind(y1, 2*x+3 + rnorm(100)) }
y2 <- c()
for(j in 1:100){ y2 <- cbind(y2, -1*x^2+30 + runif(100)) }
y1 <- cbind(y1,y2,runif(100))
pr <- prcomp(t(y1), scale = T, center=T)
plot(pr$x[,1], pr$x[,2])
which(pr$x[,2] <= -6) ## how to find outliers

pr <- prcomp(t(y1), scale = F, center=F)
max(abs(t(y1)%*%pr$rotation - pr$x))  ### relations
