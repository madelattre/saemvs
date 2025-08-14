#include "maths-utils.h"


// [[Rcpp::export]]
NumericMatrix x_per_indiv_rcpp(List x_list, List index_list) {
  int n = x_list.size();
  
  // Calcul du nombre total de colonnes dans la matrice finale
  int total_cols = 0;
  std::vector<int> widths(n);  // pour stocker la largeur de chaque bloc

  for (int i = 0; i < n; ++i) {
    widths[i] = Rf_length(index_list[i]);
    total_cols += widths[i];
  }

  // La matrice finale (en transposé directement)
  NumericMatrix mat(total_cols, n);

  int col_offset = 0;
  for (int i = 0; i < n; ++i) {
    NumericVector x = x_list[i];
    IntegerVector idx = index_list[i];
    int len = idx.size();

    for (int j = 0; j < len; ++j) {
      // -1 car R est 1-based
      mat(col_offset + j, i) = x[idx[j] - 1];
    }
    col_offset += len;
  }

  return transpose(mat);
}
