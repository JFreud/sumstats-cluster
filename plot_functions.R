
#' @param assignments df of variant assignments and cluster probabilities
#' @param n number of varaints to visualize
#' @param K number of varaints to visualize
#' @return ggplot plot object
stacked_barplot <- function(assignments, n, K,
                            do_save=F,
                            binary = T,
                            fname="plot.png",
                            orderby="assignment") {
  if (missing(n)) {n <- nrow(assignments)}
  melt_assignments <- melt(assignments[1:n,],id.vars = c("ID", "assignment", "POS"))
  if (orderby=="assignment") {
    melt_assignments <- melt_assignments[order(melt_assignments$assignment),]
  }
  else if (orderby=="position") {
    melt_assignments <- melt_assignments[order(melt_assignments$POS),]
  }
  melt_assignments$ID <- factor(melt_assignments$ID, levels = unique(melt_assignments$ID))
  if (binary) {
    model <- "binary model"
  }
  else {
    model <- "continuous model"
  }
  title <- paste("Assignment probabilities by", orderby,
                 ";", model, "K =", K)
  p <- ggplot(melt_assignments) +
    geom_bar(aes(y = value, x = ID, fill = variable),stat="identity",width=1) + 
      ggtitle(title) +
      ylab("Assignment Probability") + 
      xlab(paste("Variants ordered by", orderby))
  if (do_save) {
    ggsave(filename=paste0("plots/", fname), plot=p, width = 6, height = 4)
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

#' @param category_key mapping of trait to pre-defined trait category
#' @param trait_prop R by K matrix or df of proportions of variants in cluster for each trait
#' @param n_top top n cluster memberships to consider when adding to triangle
plot_triangles <- function(trait_prop, n_top=3, category_key) {
  triangle_list <- list()
  for (i in 1:nrow(t_prop)) {
    clusters <- which(t_prop[i,] > 0.1)
    if (length(clusters) < n_top) {
      clusters <- sort(order(t_prop[i,], decreasing=T)[1:n_top])
    }
    # clusters <- sort(order(t_prop[i,], decreasing=T)[1:n_top])
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
  colors <- setNames(c("#4e79a7","#f28e2b","#b07aa1", "#76b7b2", "#59a14f", "#9c755f"), levels(unique(category_key$`Trait category`)))
  print(triangle_list)
  p_list <- lapply(triangle_list, function(triangle) {
    clusters <- colnames(triangle)
    curr_set <- triangle %>% tibble::rownames_to_column("trait")
    curr_set$category <- as.factor(unlist(lapply(curr_set$trait, function(x) {category_key[category_key$Trait_short==x,2]})))
    colnames(curr_set) <- c("trait", "x", "y", "z", "category")
    p <- ggtern(curr_set,aes(x,y,z,label=trait)) +
      # geom_point(size=1.2, aes(color=rgb(x,y,z)),show.legend=FALSE) +
      geom_text(vjust="bottom", size=3, aes(color=category),show.legend=T) +
      scale_colour_manual(values=colors, drop=T) +
      labs(x=clusters[1],y=clusters[2],z=clusters[3],
           xarrow=paste("% variants in cluster", clusters[1]), yarrow=paste("% variants in cluster", clusters[2]), zarrow=paste("% variants in cluster", clusters[3])) +
      theme_showarrows() +
      theme_nomask()
    return(p)
  })
  return(p_list)
}
