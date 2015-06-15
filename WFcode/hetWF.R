hetWF<-function(mu=.001,initPop,G=10000,k=50){
  popsize<-sum(initPop) ## integer value of initial population size
  pop<-initPop ## vector of population size of each allelic class
    
    
  alleleV<-1:k ## vector of all possible alleles

  ## function for mulivariate hypergeometric distribution
  rmhyper <- function(B,n) ## # of balls of each color, sample size
      {
          sB <- sort(B, decreasing=TRUE) ## sort B
          ord <- abs(rank(B, ties.method="first") - (length(B)+1)) ## preserve rank
          draws <- rep(0, length(B)) ## initials vector to keep track of draws of each color

          i=1
          samp <- n
          while(sum(draws)<n)
              {
                  w <- sB[i] ## # of white balls in urn
                  b <- sum(sB[-(1:i)])
                  w+b
                  sum(draws)
                  draws[i] <- rhyper(1,w,b,samp)
                  sum(draws)
                  samp <- samp - draws[i]
                  i <- i+1
              }
          return(draws[ord])
      }

  
    
    ##  for each allelic class, deterimine the number of mutants in each group
  for(i in 1:G)
      {
          offSpring<-rmultinom(1,popsize,pop) ## sample from multinomial dist
    
          mutNum<-rpois(1,popsize*mu) ## number of mutations

          existPop <- offSpring[offSpring>0] ## populations > 1 indiv

          muts <- rmhyper(existPop,mutNum)

          mutsC <- rep(0,k)
          
          mutsC[offSpring>0] <- muts  ## mutations arising from each class
        
                    
          noMut<-offSpring-mutsC ## remove mutants from population
        
        
          if(mutNum>=1)
              {
              mutGroup<-alleleV[mutsC>0] ## find allelic classes with mutatations  
            

              ## for each allele class with mutations, determine the new allele class
              mutType<-rep(0,k)
              for(j in mutGroup)
                  {
                      mutDist <- rmultinom(1,mutsC[j], rep(1/(k-1),k-1))
                      if(j==1)
                          {
                              mutDist2 <- c(0,mutDist)
                          }
                      else if(j==k)
                          {
                              mutDist2 <- c(mutDist,0)
                          }
                      else
                          {
                              mutDist2 <- c(mutDist[1:(j-1)],0,mutDist[j:(k-1)])
                          }
                      mutType <- mutType+mutDist2
                  }
          
              newPop<-noMut+mutType 
          
              pop<-newPop
          }
      else
          {
              pop<-offSpring
          }
      }

  

  res<-as.vector(pop)                   
  return(res)
    
}

 
###


