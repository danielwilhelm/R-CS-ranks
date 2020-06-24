# load and show data
attach(pisa)
print(pisa)

# compute math ranking
math_rank <- xrank(math_score)


# ------- marginal confidence sets -------

  # compute marginal confidence sets for ranking
  CS_marg <- csranks(math_score, math_se, coverage=0.95, simul=FALSE, R=1000, seed=101)
  math_rankL_marg <- CS_marg$L
  math_rankU_marg <- CS_marg$U

  # plot ranking with marginal confidence sets
  grid::current.viewport()
  plotmarg <- plotranking(ranks=math_rank, L=math_rankL_marg, U=math_rankU_marg, popnames=jurisdiction, 
    title="Ranking of OECD Countries by 2018 PISA Math Score", subtitle="(with 95% marginal confidence sets)", colorbins=4)
  print(plotmarg)

  # save plot
  ggsave("mathmarg.pdf", plot=plotmarg)



# ------- simultaneous confidence sets -------

  # compute simultaneous confidence sets for ranking
  CS_simul <- csranks(math_score, math_se, coverage=0.95, simul=TRUE, R=1000, seed=101)
  math_rankL_simul <- CS_simul$L
  math_rankU_simul <- CS_simul$U

  # plot ranking with simultaneous confidence sets
  grid::current.viewport()
  plotsimul <- plotranking(ranks=math_rank, L=math_rankL_simul, U=math_rankU_simul, popnames=jurisdiction, 
    title="Ranking of OECD Countries by 2018 PISA Math Score", subtitle="(with 95% simultaneous confidence sets)", colorbins=4)
  print(plotsimul)

  # save plot
  ggsave("mathsimul.pdf", plot=plotsimul)