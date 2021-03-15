
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
    geom_bar(aes(y = value, x = ID, fill = variable),stat="identity",width=1) + ggtitle("Assignment probabilities for variants") +
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
}
