# joint entropy

con_table = matrix(rep(c(1,1,1,1),4), ncol = 4)

jointEntropy <- function(con_table) {
  return(-1 * sum(con_table/sum(con_table) * log2(con_table/sum(con_table)),na.rm = T))
}

jointEntropy(con_table)

# calculate mutual information

x_marginal = c(1,1,1,1)
y_marginal = c(1,1,1,1)

x_y = matrix(rep(c(0,0,0,0),4), ncol = 4)

diag(x_y) = 1

MI <- function(x_marginal, y_marginal, x_y) {
  px = x_marginal / sum(x_marginal)
  py = y_marginal / sum(y_marginal)
  pxy = x_y / sum(x_y)
  mi = 0
  for (i in 1:4) {
    for (j in 1:4) {
      if (px[i] == 0 | py[j] == 0 | pxy[i,j] == 0)
        next
      mi = mi + pxy[i,j] * log2(pxy[i,j]/(px[i]*py[j]))
    }
  }
  mi
}
MI(x_marginal, y_marginal, x_y)
