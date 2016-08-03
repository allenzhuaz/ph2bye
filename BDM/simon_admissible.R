simon.power <- function(r1,n1,r,n,p){
  pmf <- 0
  for (x1 in (r1+1):n1){
      pmf <- pmf+dbinom(x1, n1, p)*pbinom(r-x1, n-n1, p, lower.tail = F)
  }
  return(pmf)
}



output.type <- c("minimax","optimal","maxpower","admissible")
binom.design <- function(method = design.method, p0,p1,signif.level=0.05,power.level=0.85, nmax=100, output = output.type, plot.out = FALSE){
  zz <- clinfun::ph2simon(pu=p0, pa=p1, ep1 = signif.level, ep2 = 1-power.level)[[5]]
  r1 <- zz[,1]
  n1 <- zz[,2]
  r  <- zz[,3]
  n  <- zz[,4]
  result <- data.frame(cbind(zz,
                             error=mapply(simon.power, r1,n1,r,n,p=p0), power=mapply(simon.power, r1,n1,r,n,p=p1)))

  y <- result[,4:5];
  con.ind <- chull(y)[chull((y))==cummin(chull((y)))]

  plot.ph2simon <- function(x, ...) {
    xout <- x$out
    n <- nrow(xout)
    nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
    nopt1 <- min(nopt+5,n)
    nadm <- setdiff(con.ind, c(1, nopt))
    npow <- ((1:n)[result[,8]==max(result[,8])])[1]
    plot(xout[1:nopt1,4],xout[1:nopt1,5],type="l",xlab="Maximum Sample Size N" ,ylab=expression(paste("E( N | ",p[0], " )")), main = "Two-stage Designs")
    points(xout[1,4],xout[1,5],pch="M")
    points(xout[nopt,4],xout[nopt,5],pch="O")
    points(xout[nadm,4],xout[nadm,5],pch="A")
    points(xout[npow,4],xout[npow,5],pch="P")
 }


  x <- switch (output,
    minimax = {subset(result , n == min(n))},
    optimal = {subset(result , EN.p0. == min(EN.p0.))},
    maxpower  = {subset(result , power == max(power))},
    admissible = {
 #     subset(result , n >= min(n) & n <= subset(result, EN.p0. == min(EN.p0.))$'n')[con.ind,]
 #     result[con.ind,]
        subset(result[con.ind,],n >= min(n) & n <= subset(result, EN.p0. == min(EN.p0.))$'n')
      }

    )
  rownames(x) <- make.names(c("Optimal",rep("Admissible",nrow(x)-2),"Minimax"),unique=T)
  
  
if (plot.out==TRUE){
  plot.ph2simon(clinfun::ph2simon(pu=p0, pa=p1, ep1 = signif.level, ep2 = 1-power.level,nmax = nmax))
}
  return(x)
}

