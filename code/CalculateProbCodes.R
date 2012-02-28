barcode=function(y,l=NULL,gc=NULL,conv=.01) {
	X=1.*cbind(y>median(y),y<=median(y),l,gc)

	pi_h=.5
	beta=solve(t(X)%*%X)%*%t(X)%*%y
	sigma2_1=var(y[y>median(y)])
	sigma2_2=var(y[y<=median(y)])
	x_1=cbind(matrix(c(1,0),length(y),2,byrow=T),X[,-c(1:2)])
	x_2=cbind(matrix(c(0,1),length(y),2,byrow=T),X[,-c(1:2)])
	
	theta=c(pi_h,beta,sigma2_1,sigma2_2)
	thetaold=theta+1000
	i<- 0

	while(max(abs((theta-thetaold)/thetaold))>conv)
    {
      thetaold=theta	

      #E-step
      gamma <- (pi_h*dnorm(y,x_1%*%beta,sqrt(sigma2_1)))/(pi_h*dnorm(y,x_1%*%beta,sqrt(sigma2_1))+(1-pi_h)*dnorm(y,x_2%*%beta,sqrt(sigma2_2)))
      x=cbind(gamma,1-gamma,X[,-c(1:2)])

      #M-step
      Sigma <- (gamma*sigma2_1+(1-gamma)*sigma2_2)^(-1)
      beta <- solve(t(X)%*%(Sigma*X))%*%t(X)%*%(Sigma*y)
	
      sigma2_1 <- (t(y-X%*%beta)%*%(gamma*(y-X%*%beta)))/(sum(gamma))
      sigma2_2 <- (t(y-X%*%beta)%*%((1-gamma)*(y-X%*%beta)))/(sum(1-gamma))
	
      pi_h <- mean(gamma)

      theta=c(pi_h,beta,sigma2_1,sigma2_2)
	
      i <- i+1
      print(paste("Iteration", i))
      print(max(abs((theta-thetaold)/thetaold)))

      if (i == 1000)
      {
        print("Stopped at 1000 iterations")
        break
      }
    }

    return(gamma)
}

inFilePath = commandArgs()[7]
outFilePath = commandArgs()[8]
if (length(commandArgs()) > 8)
{
  conv = as.numeric(commandArgs()[9])
} else
{
  conv = 0.01
}

print("Reading data")
data = read.table(inFilePath, sep="\t", header=FALSE, stringsAsFactors=FALSE, row.names=1, quote="\"")

expr = data[,1]

lengths = NULL
gc = NULL

if (ncol(data) >= 3)
{
  print("Getting length and gc information")
  lengths = as.numeric(data[,2])
  gc = as.numeric(data[,3]) / lengths

  # Don't adjust for length if they are all the same
  if (length(unique(lengths))==1)
  {
    print("Lengths are all the same")
    lengths = NULL
  }
}

print("Calculating bar codes")
barcodes = barcode(expr, l=lengths, gc=gc, conv=conv)
outData = cbind(rownames(data), round(barcodes, 6))

write.table(outData, outFilePath, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
