# HMM example e.g. two state GC-rich detection

# Simulate data

state <- 1
N <- 1000
data <- sapply(1:N, FUN=function(i){
  # transition
  if( state ==  1){
    state <<- ifelse(runif(1) > 0.95, 2, 1)
  } else{
    state <<- ifelse(runif(1) > 0.95, 1, 2)
  }
  
  # emission
  if( state == 1){
    ifelse(runif(1)>0.7, 1, 2)
  } else{
    ifelse(runif(1)>0.7, 2, 1)
  }
})

data <- data.frame(matrix( c(rep(NA,N),data) , 1, 2*N))

# Train model
varDim <- rep(2, 2*N)
facPot <- c(list(matrix(c(0.4,0.6),1,2)),
            list(matrix(c(0.6,0.4,0.45,0.55), 2, 2, byrow = T)),
            list(matrix(c(0.6,0.4,0.7,0.3), 2, 2, byrow = T)))
potMap <- c( 1, rep(2, N-1), rep(3, N))
facNbs <- c(list(c(1L)), lapply(2:N, FUN=function(i){c(i-1, i)}), lapply(1:N, FUN=function(i){c(i, i+N)}) )
mydfg <- dfg(varDim, facPot, facNbs, potMap = potMap)
train(data, mydfg)
potentials(mydfg)
