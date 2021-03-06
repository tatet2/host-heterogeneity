## generate n samples from a K allle (halploid)population with k possible alleles



coal<-function(n, theta, k) 
{

  ## function to make stucture of tree
  coaltree<-function(n){
     

    ## sample size (# of genes)
    ## function replace sample function for when length(x) = 1 (from help file)
    resample <- function(x, ...) x[sample.int(length(x), ...)]

  
    ##define things
    headers<-c("time", "desc1", "desc2", "ancestor")
    nodes<-array(NA, c(2*n-1,4))
    dimnames(nodes)<-list(paste("node[",1:(2*n-1), "]" ,sep=""),c("time", "desc1", "desc2", "ancestor"))
    nodes[1:n,1:3]<-0


    ## generate the times of the internal nodes
    t=0
    for(k in n:2){
    t=t+rexp(1, rate=choose(k,2))
    nodes[2*n-k+1,1]<-t

  }

    ##generate topology

    nodeList<-1:n
    for(k in n:2){
      pick <- resample(nodeList,1) ## pick a lineage
      nodeList <- nodeList[nodeList!=pick]  ## remove pick from  nodeList
      nodes[pick,4] <- 2*n-k+1 ## assigns ancestor to randomly choosen lineage
      nodes[2*n-k+1,2] <- pick  ## assigns decendent to 1st node
      pick2<-resample(nodeList,1) ## pick a 2nd lineage to coalesce with 1st
      nodeList <- nodeList[nodeList!=pick2] ## remove pick2 from  nodeList
      nodes[pick2,4] <- 2*n-k+1 ## assigns ancestor to 2nd lineage
      nodes[2*n-k+1,3] <- pick2 ## assigns 2nd decendent to 1st node
      nodeList <- c(nodeList, 2*n-k+1) ## add ancestor to nodeList
    
    }

    return(nodes)
  }


  tree<-as.data.frame(coaltree(n=n))

  tree$state<-c(rep(NA,2*n-2),sample(1:k,1))
  K <- 1:k
  p <- 1/k
  dm <- matrix(0,nrow=k, ncol=k)
  dm[col(dm) != row(dm)] <- p
  matrix.power <- function(mat, n)
    {
    # test if mat is a square matrix
    # treat n < 0 and n = 0 -- this is left as an exercise
    # trap non-integer n and return an error
    if (n == 1) return(mat)
    result <- diag(1, ncol(mat))
    while (n > 0) {
      if (n %% 2 != 0) {
        result <- result %*% mat
        n <- n - 1
      }
      mat <- mat %*% mat
      n <- n / 2
    }
    return(result)
  }

  for(i in (2*n-1):(n+1))
    {

      ##calculate state for 1st branch
      branch <- tree$time[i] - tree$time[tree$desc1[i]] ## length of branch
      Mt <- rpois(1,(theta/2)*branch)## number of mstat mutations along branch
      W <- matrix.power(dm, Mt)
      tree$state[tree$desc1[i]] <- sample(K,1, prob=W[,tree$state[i]]) ## update tree state




    
##calculate state for 2nd branch, same as above
      branch2 <- tree$time[i] - tree$time[tree$desc2[i]]
      Mt2 <- rpois(1,(theta/2)*branch2)
      W <- matrix.power(dm, Mt2)
      tree$state[tree$desc2[i]] <- sample(K,1,prob=W[,tree$state[i]]) ## update tree state


      
    }

##  finalState <- count(tree$state[1:n])
##  pop <- rep(0,k)
##  pop[finalState$x] <- finalState$freq  

  res <- tabulate(tree$state[1:n], nbins=k)
  
  return(res)

}



## test code not run

## coaltest <- function(ntest){

##     res <- matrix(NA, nrow=ntest, ncol=k)
##     for(i in 1:ntest){
##         res[i,] <- coal(n=10,theta=1,k=50)
##     }
##     return(res)
## }
