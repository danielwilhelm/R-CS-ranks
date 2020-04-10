# load and show data
attach(pisa)
print(pisa)

# compute reading ranking
reading_rank <- xrank(reading_score, na.rm=TRUE)


# ------- marginal confidence sets -------

  # compute marginal confidence sets for ranking
  CS_marg <- csranks(reading_score, reading_se, coverage=0.95, simul=FALSE, R=1000, na.rm=TRUE, seed=101)
  reading_rankL_marg <- CS_marg$L
  reading_rankU_marg <- CS_marg$U

  # plot ranking with marginal confidence sets
  grid::current.viewport()
  plotmarg <- plotranking(ranks=reading_rank, L=reading_rankL_marg, U=reading_rankU_marg, popnames=jurisdiction[!is.na(reading_score)], 
    title="Ranking of OECD Countries by 2018 PISA Reading Score", subtitle="(with 95% marginal confidence sets)", 
    caption="Note: Spain's reading score is missing.")
  print(plotmarg)

  # save plot
  ggsave("readingmarg.pdf", plot=plotmarg)



# ------- simultaneous confidence sets -------

  # compute simultaneous confidence sets for ranking
  CS_simul <- csranks(reading_score, reading_se, coverage=0.95, simul=TRUE, R=1000, na.rm=TRUE, seed=101)
  reading_rankL_simul <- CS_simul$L
  reading_rankU_simul <- CS_simul$U

  # plot ranking with simultaneous confidence sets
  grid::current.viewport()
  plotsimul <- plotranking(ranks=reading_rank, L=reading_rankL_simul, U=reading_rankU_simul, popnames=jurisdiction[!is.na(reading_score)], 
    title="Ranking of OECD Countries by 2018 PISA Reading Score", subtitle="(with 95% simultaneous confidence sets)",
    caption="Note: Spain's reading score is missing.")
  print(plotsimul)

  # save plot
  ggsave("readingsimul.pdf", plot=plotsimul)