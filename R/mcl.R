mcl <-
function(x, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = FALSE, max.iter = 100, ESM = FALSE ){
    if (addLoops) diag(x) <- 1.0

    # normalize the weights in adjacency matrix
    adj.norm <- x %*% diag(1/colSums(x))

    # do MCL iterations
    a <- 1
    repeat{

      # expansion and inflation of the current adjacency matrix
      expans <- adj.norm %^% expansion
      infl <- expans ^ inflation
      # normalize the new adjacency matrix
      infl.norm <- infl %*% diag(1/colSums(infl))

      if(identical(infl.norm,adj.norm)) {
        ident <- TRUE
        break
      }

      if(a==max.iter) {
        ident <- FALSE
        a <- a+1
        break
      }
      adj.norm <- infl.norm
      a<- a+1
    }

    if(!is.na(infl.norm[1,1]) & ident){

      # remove rows containing only zero elements
      neu <- infl.norm[rowSums(abs(infl.norm)) > 0.0,]

      for(i in 1:nrow(neu)){
        for(j in 1:ncol(neu)) {
          if((neu[i,j] < 1) & (neu[i,j] > 0)){
            neu[,j] <- 0
            neu[i,j] <- 1
          }
        }
      }

      for(i in 1:nrow(neu)){
        for (j in 1:ncol(neu)){
          if(neu[i,j] != 0){
            neu[i,j] <- i
          }
        }
      }

      # assign cluster indexes to each node
      ClusterNummern <- sum(neu[,1])
      for(j in 2:ncol(neu)){
        ClusterNummern <- c(ClusterNummern,sum(neu[,j]))
      }
    }

    ifelse(!(!is.na(infl.norm[1,1]) & ident), output <- paste("An Error occurred at iteration", a-1),
      {
      if(!allow1){
        # collapse all size 1 clusters into one with index 0
        dub <- duplicated(ClusterNummern) + duplicated(ClusterNummern,fromLast = T)
        ClusterNummern[!dub] <- 0
      }

      #### dimnames for infl.norm
      dimnames(infl.norm) <- list(1:nrow(infl.norm), 1:ncol(infl.norm))

      output <- list()
      output[[1]] <- length(table(ClusterNummern))
      output[[2]] <- a-1 
      output[[3]] <- ClusterNummern
      output[[4]] <- infl.norm

      names(output) <-c("K", "n.iterations","Cluster",
                        "Equilibrium.state.matrix")
    }
    )
  ifelse(ESM==TRUE,return(output),return(output[-4]))
}
