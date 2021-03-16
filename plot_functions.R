
#' @param assignments df of variant assignments and cluster probabilities
#' @param n number of varaints to visualize
#' @return ggplot plot object
stacked_barplot <- function(assignments, n, do_save=F, fname="plot.png",orderby="assignment") {
  if (missing(n)) {
    n <- nrow(assignments)
  }
  melt_assignments <- melt(assignments[1:n,],id.vars = c("ID", "assignment", "POS"))
  if (orderby=="assignment") {
    melt_assignments <- melt_assignments[order(melt_assignments$assignment),]
  }
  else if (orderby=="position") {
    melt_assignments <- melt_assignments[order(melt_assignments$POS),]
  }
  melt_assignments$ID <- factor(melt_assignments$ID, levels = unique(melt_assignments$ID))
  p <- ggplot(melt_assignments) +
    geom_bar(aes(y = value, x = ID, fill = variable),stat="identity",width=1) + ggtitle(paste("Assignment probabilities for variants by", orderby)) +
    ylab("Assignment Probability") + xlab(paste("Variants ordered by", orderby))
  if (do_save) {
    ggsave(filename=paste0("plots/", fname), plot=p)
  }
  return(p)
}

#' @param res gibbs sampler result
plot_gibbs_chain <- function(res, K) {
  plot(res$P[1,1,],type="l", ylim=c(0,1), main="P")
  for (i in 1:K) {
    for (j in 1:35) {
      lines(res$P[i,j,],col=i+j)
    }
  }
  plot(res$pi[,1],type="l", ylim=c(0,1), main="pi")
  for (i in 2:6) {
    lines(res$pi[,i],col=i)
  }
  plot(res$p.z.given.x[1,1,],type="l", ylim=c(0,1), main="p z given x")
  for (i in 1:nrow(res$p.z.given.x)) {
    for (j in 1:20) {
      lines(res$p.z.given.x[i,j,],col=i+j)
    }
  }
}

#' @param trait_prop R by K matrix or df of proportions of variants in cluster for each trait
plot_triangles <- function(trait_prop) {
  triangle_list <- list()
  for (i in 1:nrow(t_prop)) {
    clusters <- which(t_prop[i,] > 0.1)
    if (length(clusters) < 3) {
      clusters <- sort(order(t_prop[i,], decreasing=T)[1:3])
    }
    triangles <- t(combn(clusters,3))
    # add to traits in triangles
    for (j in 1:nrow(triangles)) {
      name <- paste(triangles[j,], collapse='')
      if (name %in% names(triangle_list)) {
        # add cluster values to triangle
        triangle_list[[name]] <- rbind(triangle_list[[name]], t_prop[i,triangles[j,]])
      } else {
        # create new triangle
        triangle_list[[name]] <- t_prop[i,triangles[j,]]
      }
    }
  }
  p_list <- lapply(triangle_list, function(triangle) {
    curr_set <- triangle %>% tibble::rownames_to_column("trait")
    colnames(curr_set) <- c("trait", "x", "y", "z")
    p <- ggtern(curr_set,aes(x,y,z,label=trait)) +
      geom_point(size=1.2, aes(color=rgb(x,y,z)),show.legend=FALSE) +
      geom_text(vjust=1, size=2, aes(color=rgb(x,y,z)),show.legend=F, srt=20) +
      labs(x="1",y="3",z="5",
           xarrow="% variants in cluster 1", yarrow="% variants in cluster 3", zarrow="% variants in cluster 5") +
      theme_showarrows() +
      theme_nomask()
    return(p)
  })
  return(p_list)
}
