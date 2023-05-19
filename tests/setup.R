compare_for_tests <- function(i,j,v,omega=1){
  omega * (v[i] <= v[j]) + (1-omega) * (v[i] < v[j])
}
