mcl <-
function(x, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = FALSE, max.iter = 100L, tol=1E-5, ESM = FALSE, verbose=FALSE){
    if (addLoops) diag(x) <- 1.0

    # normalize the weights in adjacency matrix
    adj.norm <- x %*% diag(1/colSums(x))

    # do MCL iterations
    if (verbose) message("Starting MCL iterations...")
    rel_err <- NA
    for (niter in 1:max.iter) {
      # expansion and inflation of the current adjacency matrix
      expans <- adj.norm %^% expansion
      infl <- expans ^ inflation
      # normalize the new adjacency matrix
      infl.norm <- infl %*% diag(1/colSums(infl))

      rel_err = norm(infl.norm - adj.norm, "F")/norm(adj.norm, "F")
      if (verbose) {
        message("MCL iteration #", niter, " rel.error=", rel_err, "...")
      }
      if(rel_err <= tol) {
        if (verbose) message("MCL iterations converged")
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

    if (verbose) message("Generating MCL clusters...")
    # remove rows containing only zero elements
    neu <- infl.norm[rowSums(abs(infl.norm)) > 0.0,]

    # assign cluster indexes to each node
    ClusterNummern <- rep.int(-1, ncol(neu))
    for(i in 1:nrow(neu)){
      for(j in 1:ncol(neu)) {
        if(neu[i,j] > 0){
          neu[,j] <- 0 # exclude the node from further cluster assignment
          ClusterNummern[j] <- i
        }
      }
    }

    if(!allow1){
      # collapse all size 1 clusters into one with index 0
      dub <- duplicated(ClusterNummern) + duplicated(ClusterNummern,fromLast = T)
      ClusterNummern[!dub] <- 0
    }
    # recode clusters numbers to be in 1:N range (or 0:N if there's collapsed cluster)
    ClusterNummern <- match(ClusterNummern, sort(unique(ClusterNummern))) - ifelse(0 %in% ClusterNummern, 1, 0)

    output$K <- length(table(ClusterNummern))
    output$n.iterations <- niter
    output$Cluster <- ClusterNummern
    if (ESM) {
      output$Equilibrium.state.matrix <- infl.norm
    }

    return(output)
}
