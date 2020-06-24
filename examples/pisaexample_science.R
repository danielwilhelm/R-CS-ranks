# load and show data
attach(pisa)
print(pisa)

# compute science ranking
science_rank <- xrank(science_score)


# ------- marginal confidence sets -------

  # compute marginal confidence sets for ranking
  CS_marg <- csranks(science_score, science_se, coverage=0.95, simul=FALSE, R=1000, seed=101)
  science_rankL_marg <- CS_marg$L
  science_rankU_marg <- CS_marg$U

  # plot ranking with marginal confidence sets
  grid::current.viewport()
  plotmarg <- plotranking(ranks=science_rank, L=science_rankL_marg, U=science_rankU_marg, popnames=jurisdiction, 
    title="Ranking of OECD Countries by 2018 PISA Science Score", subtitle="(with 95% marginal confidence sets)", colorbins=4)
  print(plotmarg)

  # save plot
  ggsave("sciencemarg.pdf", plot=plotmarg)



# ------- simultaneous confidence sets -------

  # compute simultaneous confidence sets for ranking
  CS_simul <- csranks(science_score, science_se, coverage=0.95, simul=TRUE, R=1000, seed=101)
  science_rankL_simul <- CS_simul$L
  science_rankU_simul <- CS_simul$U

  # plot ranking with simultaneous confidence sets
  grid::current.viewport()
  plotsimul <- plotranking(ranks=science_rank, L=science_rankL_simul, U=science_rankU_simul, popnames=jurisdiction, 
    title="Ranking of OECD Countries by 2018 PISA Science Score", subtitle="(with 95% simultaneous confidence sets)", colorbins=4)
  print(plotsimul)

  # save plot
  ggsave("sciencesimul.pdf", plot=plotsimul)