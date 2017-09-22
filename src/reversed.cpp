#include "utils.h"

// [[Rcpp::export]]
arma::vec sphericalFromCartesian(const arma::vec& coords) {

    int n = coords.n_elem;
    vec angles(n - 1);

    for (int i = 0; i < coords.n_elem - 1; i++) {
        vec sq = pow(coords.subvec(i, n-1), 2);
        double ssum = sum(sq);
        if (sqrt(ssum) == 0) {
            angles[i] = 0;
        } else {
            angles[i] = acos(coords[i] / sqrt(ssum));
        }
    }

    return(angles);
}



// [[Rcpp::export]]

arma::vec parFromCoefMatrix(const arma::mat& cm, int dims) {
    mat angles(dims - 1, dims);
    mat pre_matrix(dims - 1, dims);
    mat cm_t = cm.t();

    for (int i = 0; i < dims; i++) {
        angles.col(i) = sphericalFromCartesian(cm_t.col(i));
    }
    // Rcout << angles;

    for (int i = 0; i < dims - 1; i++) {
        int k = dims - 1 - i;
        beta_distribution<double> bdist(k / 2.0, k / 2.0);
        for (int j = 0; j < dims; j++) {
            pre_matrix(i, j) = cdf(bdist, cos(angles(i, j)) / 2.0 + 0.5);
            // inverse of
            // pre_matrix(j, i) = acos((quantile(bdist, pre_matrix(j, i)) - 0.5) * 2.0);
        }
    }

    pre_matrix.reshape(dims * (dims - 1), 1);
    return pre_matrix.col(0);

}
