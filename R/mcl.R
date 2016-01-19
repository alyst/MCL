mcl <-
function(x, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = FALSE, max.iter = 100L, tol=1E-5, ESM = FALSE ){
    if (addLoops) diag(x) <- 1.0

    # normalize the weights in adjacency matrix
    adj.norm <- x %*% diag(1/colSums(x))

    # do MCL iterations
    rel_err <- NA
    for (niter in 1:max.iter) {
      # expansion and inflation of the current adjacency matrix
      expans <- adj.norm %^% expansion
      infl <- expans ^ inflation
      # normalize the new adjacency matrix
      infl.norm <- infl %*% diag(1/colSums(infl))

      rel_err = norm(infl.norm - adj.norm, "F")/norm(adj.norm, "F")
      if(rel_err <= tol) {
        break
      }
      adj.norm <- infl.norm
    }

    output <- list()
    if (is.na(infl.norm[1,1])) {
      output$status <- paste0("Error: matrix norm is NA")
    } else if (rel_err > tol) {
      output$status <- paste0("Error: the algorithm did not converge")
    } else {
      # converged, prepare the results
      output$status <- "OK"
    }
    output$relative.error <- rel_err

    #### dimnames for infl.norm
    dimnames(infl.norm) <- list(1:nrow(infl.norm), 1:ncol(infl.norm))

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

    if(!allow1){
      # collapse all size 1 clusters into one with index 0
      dub <- duplicated(ClusterNummern) + duplicated(ClusterNummern,fromLast = T)
      ClusterNummern[!dub] <- 0
    }

    output$K <- length(table(ClusterNummern))
    output$n.iterations <- niter
    output$Cluster <- ClusterNummern
    if (ESM) {
      output$Equilibrium.state.matrix <- infl.norm
    }

    return(output)
}
