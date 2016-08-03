

library(animation)
ani.options(convert='C:\\Users\\yalin.zhu\\Documents\\ImageMagick-7.0.1-10-portable-Q16-x64\\convert.exe')

  
  BB.sim <- function(a =0.5, b =0.5, M, N, p=0.5, alpha=0.05, seed=123,  output = FALSE ){
    # a,b  shape parameters of beta prior
  
    set.seed(seed)
    options(digits = 3)
    r <- rbinom(n = M,size = N,prob = p)
# r=rep(0,12)
 #   saveGIF({
   
      cat("Prior: Beta(",a ,",", b, ") \n\n", sep = "")
      for ( i in 1:M){
      x <- seq(0,1,0.01)
      ## calculate the prior and posterior density
      y.prior <- dbeta(x,a,b)
      mu0 <- a/(a+b)
      a1 <- a+r[i]; b1 <- b+N-r[i]
      mu1 <- a1/(a1+b1)
      y.post <- dbeta(x,a1, b1)
      
      lh <-  dbinom(r[i], N, x) * (N+1)
      if (output == TRUE) {
      LB <- qbeta(alpha/2,a1,b1); UB <- qbeta(alpha/2,a1,b1,lower.tail = F)   
      # Updating prior parameters
      cat("======== Cohort Number: ",i, " ======== \n", sep = "")
      cat("Observations -- Sample Size: ", i*N, "(",N,")", "  ||  Number of Response: ", sum(r[1:i]), "(",r[i],")", "  ||  Number of Failure: ", i*N-sum(r[1:i]), "(",N-r[i],")","\n", sep = "")
      cat("  Observed Response Rate: ", sum(r[1:i])/(i*N), "\n", sep = "")
      cat("Posterior: Beta(",a1 ,",", b1, ") \n", sep = "")
      cat("  Posterior Mean: ", mu1, ",  Difference between posterior and true response rate: ", abs(mu1-sum(r[1:i])/(i*N)),"\n", sep = "")
      cat("  ",(1-alpha)*100,"% Credible Interval: (", LB, ",", UB,") \n\n", sep="")
      }
          
      ylim <- range(0, min(20,max(y.prior)), min(20,max(y.post)), max(lh))
      plot(x, y.prior, type="l",  col="blue",xlim=c(0,1), ylim = ylim, xlab="Response Rate", ylab="Density" , main=paste(N*i,"patients are included",sep =" "))
      ## add a vertical line for the mean prior response probability
   #   abline(v=mu0,lty=1,col="blue")
      ## add comment text
   #   mtext( bquote(E(pi[prior]) == .(mu0)), side = 1, line = 0.5,at = mu0,col="blue")
     

      ## plot the prior density
      
      ## plot the posterior density
      lines(x,y.post, col = "red")
      ## add a vertical line for the estimation of the mean posterior response probability
   #   abline(v=mu1,lty=2, col = "red" )
      ## add comment text
#      mtext( bquote(E(pi[post]) == .(mu1)), side = 3, line = 0,at =mu1,col="red")
#      abline(v=qbeta(alpha/2,a1,b1),lty=2, col='dark green')
#      abline(v=qbeta(alpha/2,a1,b1,lower.tail = F),lty=2, col='dark green')
      
      lines(x,lh, col="dark green")
#      legend("topright",c("Prior", "Posterior","Cred.Int","Likelihood"), lty=c(1, 2, 2,1), col=c("blue","red", "dark green","dark green"))
      legend("topright",c("Prior", "Posterior","Likelihood"), lty=1, col=c("blue","red", "dark green"))
      
      
      a <- a1; b <- b1
      }
 # }, interval = 3)
  }
  
  
  
