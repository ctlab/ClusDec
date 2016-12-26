#include "utils.h"

// [[Rcpp::export]]
arma::vec cartseianFromSpherical(const arma::vec& angles) {
    vec angles2(angles.n_elem + 2);
    angles2[0] = datum::pi / 2;
    angles2.subvec(1, angles.n_elem) = angles;
    angles2[angles.n_elem + 1] = 0;
    vec angles_cos = cos(angles2);
    vec angles_sin = sin(angles2);
    vec results(angles.n_elem + 1);
    for (int i = 0; i <= angles.n_elem; i++) {
        results[i] = prod(angles_sin.subvec(0, i)) * angles_cos[i + 1];
    }
    return results;
}

// [[Rcpp::export]]
List deconvolve(const arma::mat& mixedData,
                const arma::mat& gg,
                const arma::mat& coefMatrix,
                const arma::vec& coef,
                int dims) {
    //Rcout << coefMatrix;
    colvec planeTransformation = solve(coefMatrix.t(), vec(dims, fill::ones));
    // Rcout << planeTransformation;
    colvec abc = planeTransformation / coef;
    // Rcout << abc;
    arma::mat pseudoBasis(coefMatrix.n_rows, coefMatrix.n_cols);
    for (int i = 0; i < pseudoBasis.n_rows; i++) {
        pseudoBasis.row(i) = coefMatrix.row(i) * abc[i];
    }
    // Rcout << pseudoBasis;
    // Rcout << gg;

    mat props = fcnnls_c(pseudoBasis, gg);
    mat a = props.t();
    mat b = mixedData.t();

    mat basis_results = fcnnls_c(a, b).t();
    mat results = basis_results * props;


    return List::create(
        _["basis"] = basis_results,
        _["proportions"] = props,
        _["resutls"] = results
    );
    // return List::create();
}
// [[Rcpp::export]]
arma::mat getCoefMatrix(const arma::vec& par, int dims) {
    int nrow = par.n_elem / (dims - 1);
    mat pre_matrix;
    pre_matrix.insert_cols(0, par);
    pre_matrix.reshape(dims - 1, nrow);
    pre_matrix = pre_matrix.t();
    for (int i = 0; i < dims - 1; i++) {
        int k = dims - 1 - i;
        beta_distribution<double> bdist(k / 2.0, k / 2.0);
        for (int j = 0; j < nrow; j++) {
            pre_matrix(j, i) = acos((quantile(bdist, pre_matrix(j, i)) - 0.5) * 2.0);
        }
    }
    //Rcout << pre_matrix;
    mat results(nrow, dims);
    for (int i = 0; i < nrow; i++) {
        // Rcout << "subelems:\n" << i * (dims - 1) << " " << (i + 1) * dims - 1 << "--\n";
        results.row(i) = cartseianFromSpherical(pre_matrix.row(i).t()).t();
    }
    return results;
}



double score(const arma::mat& mixedData,
             const arma::mat& gg,
             const arma::vec& coef,
             const arma::vec& par,
             int dims) {

    mat cm = getCoefMatrix(par, dims);
    mat results = deconvolve(mixedData, gg, cm, coef, dims)[2];
    mat error = log2(mixedData + 1) - log2(results + 1);
    return sum(sum(pow(error, 2)));
}

// [[Rcpp::export]]
double score_export(const arma::vec& par,
                    const arma::mat& mixedData,
                    const arma::mat& gg,
                    const arma::vec& coef,
                    int dims) {
    std::ostream nullstream(0);
    set_stream_err2(nullstream);
    return(score(mixedData, gg, coef, par, dims));
}

