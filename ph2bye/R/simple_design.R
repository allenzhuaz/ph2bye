## create a function using parameters: prior mean (m) and varaince (v)
## n is the number of patients sequentially enrolled in trials; s is the number of response after inclusion of each new patient
bayes.desgin <- function(mu0,sigma0,n=0,r=0,stop.rule="futility", ymax, add.size=5, alpha=0.05,
                         Rmax=NULL,Rmin=NULL,tau1=0.9,tau2=0.9,tau3=0.9,tau4=0.9){

  ############################
  ## Part 0: Basic settings ##
  ############################

  # if (dist="beta")......
  a <- mu0^2*(1-mu0)/sigma0^2 - mu0; b <- mu0*(1-mu0)^2/sigma0^2 - (1-mu0)
  # calculate the posterior mean according to known information
  postm <- (a+r)/(a+b+n)
  result <- list("para.a"=a, "para.b"=b, "poterior mean"= postm)
  print(result)



  ################################################################
  ## Part 1: Plot the prior and posterior distribution denstity ##
  ################################################################
  ## layout of the plots
  library(animation)
  ani.options(interval=1) ## set speed of the animation

  #saveHTML({
    #saveGIF({
    for ( i in seq_along(n)){
      x <- seq(0,1,0.01)
      ## calculate the prior and posterior density
      y.prior <- dbeta(x,a,b)
      y.post <- dbeta(x,a+r[i], b+n[i]-r[i])
      ############## Adding animation code part###################
      options(digits = 3)
      ## plot the prior density

      plot(x, y.prior, type="l", lty=1,lwd=2, col="blue",xlim=c(0,1),ylim = c(0,ymax) ,xlab="Response Rate", ylab="Density" , main=paste(n[i],"patients are included",sep =" "))
      ## add a vertical line for the mean prior response probability
      abline(v=mu0,lty=1,lwd=2,col="blue")
      ## add comment text
 #     text(mu0+0.1,ymax-0.1, bquote(E(pi[prior]) == .(mu0)), col="blue")
      mtext( bquote(E(pi[prior]) == .(mu0)), side = 1, line = 0.5,at = mu0,col="blue")
      ## plot the posterior density
      lines(x,y.post,lty=2,lwd=2, col = "red")
      ## add a vertical line for the estimation of the mean posterior response probability
      abline(v=(a+r[i])/(a+b+n[i]),lty=2, lwd=2,col = "red" )
      ## add comment text
      mtext( bquote(E(pi[post]) == .((a+r[i])/(a+b+n[i]))), side = 3, line = 0,at = (a+r[i])/(a+b+n[i]),col="red")
      abline(v=qbeta(alpha/2,a+n[i],b+n[i]-r[i]),lty=2, col='dark green')
      abline(v=qbeta(alpha/2,a+n[i],b+n[i]-r[i],lower.tail = F),lty=2, col='dark green')
      legend("topright",c("Prior", "Posterior"), lty=c(1, 2), col=c("blue","red"),lwd=2)
    }
#  })

  ani.pause()

  ################################################################
  ## Part 2: Sample size determination rules                 #####
  ################################################################
  ## Rule 1
  # obj <- function(N){
  #  return(tau1-(pbeta(0.76,a+N,b+N-N*Rmin)-pbeta(0.42,a+N,b+N-N*Rmin)))
  #  }
  # out = optimize(3,fn =obj )

  ################################################################
  ## Part 3: Stopping rule based on pre-specified threshold  #####
  ################################################################

  library(VGAM)
  # Rmax <- mu0 #APL data

  if (stop.rule=="efficacy"){
    ### rule 1 for efficacy
    n1 <- min(which(pbeta(Rmax,a+r, b+n-r, lower.tail = F) > tau1))
    print(paste("Stop the trial for efficacy after the inclusion of",n1, "patients."))

    ## could add a gradually animation
    par(mfrow=c(1,2))
    plot(n,pbeta(Rmax,a+r, b+n-r, lower.tail=F),type="o", xlab="Number of patients", ylab="Posterior Probability")
    abline(h=tau1, lty=2)
    legend("bottomright",c("Posterior tial prob","Stopping threshold"), lty=c(1,2),cex = 0.8)
    ## rule 3 for efficacy

    cdf <- pbetabinom.ab(q = c(0:add.size),size = add.size, shape1 = a+max(r), shape2 = b+max(n)-max(r) )
    plot(c(0:add.size),cdf,type="S")
    ## rule 5
  }

  else if (stop.rule=="futility"){
    ## rule 2 for futility
    # Rmax = mu0  MM data
    n2 <- min(which(pbeta(Rmin,a+r, b+n-r) > tau2))
    print(paste("Stop the trial for futility after the inclusion of",n2, "patients."))

    par(mfrow=c(1,2))
    plot(n,pbeta(Rmin,a+r, b+n-r),type="o", xlab="Number of patients", ylab="Posterior Probability")
    abline(h=tau2, lty=2)
    legend("bottomright",c("Posterior tial probability","Stopping threshold"), lty=c(1,2),cex=0.8)
    ## rule 4 for futility
    cdf <- pbetabinom.ab(q = c(0:add.size),size = add.size, shape1 = a+max(r), shape2 = b+max(n)-max(r) ) ## rule 3
    plot(c(0:add.size),cdf,type="S")


  }


  else if (stop.rule=="satisfactory"){
    ## rule 5
  }

  else{
    print("Warning: please assign a stopping rule!")
  }

}


