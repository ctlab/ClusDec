// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
using namespace boost::math;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double explainedVariance(const arma::vec& x, const arma::vec& y, double slope) {
    mat k(x.n_elem, 2);
    k.col(0) = x;
    k.col(1) = y;
    // Rcout << "input:\n";
    // Rcout << k;
    double p1 = std::sqrt(1.0 / (slope * slope + 1.0));
    double p2 = slope * p1;
    mat22 rotation;
    rotation.col(0) = colvec({p1, p2});
    rotation.col(1) = colvec({-p2, p1});
    // Rcout << "Rotation:\n";
    // Rcout << rotation;
    k = k * rotation;
    // Rcout << "Rotated:\n";
    // Rcout << k;
    mat covar = cov(k);
    // Rcout << "Covar:\n";
    //Rcout << covar;
    double val = covar.at(0, 0) / (covar.at(0, 0)+ abs(covar.at(0, 1)) + covar.at(1, 1));
    // Rcout << "Value: " << val << "\n";
    return val;
}


// [[Rcpp::export]]
List pairwiseLinearFit(const arma::mat& X) {
    int genes_count = X.n_cols;
    int sample_count = X.n_rows;
    mat resid(genes_count, genes_count);
    mat slopes(genes_count, genes_count);
    mat r2(genes_count, genes_count);
    vec tsss(genes_count);
    mat ev(genes_count, genes_count);
    // mat pvals(genes_count, genes_count);

    int df = sample_count - 1;
    //fisher_f dist(1, df);

    for (int i = 0; i < genes_count; i++) {
        mat x_intercepted(sample_count, 1);
        x_intercepted.col(0) = X.col(i);

        mat coef = solve(x_intercepted, X);

        colvec angles = coef.row(0).t();
        slopes.col(i) = angles;

        mat residuals = X - x_intercepted * coef;
        vec rss(genes_count);
        for (int j = 0; j < genes_count; j++) {
            rss[j] = dot(residuals.col(j), residuals.col(j));
        }
        resid.col(i) = rss;
        vec x_int(sample_count);
        double mean_x_intercept = mean(x_intercepted.col(0));
        for (int j = 0; j < sample_count; j++) {
            x_int[j] = x_intercepted.col(0)[j] - mean_x_intercept;
        }
        double tss = dot(x_int, x_int);
        tsss[i] = tss;
    }

    for (int i = 0; i < genes_count; i++) {
        r2.col(i) = ones(genes_count) - resid.col(i) / tsss;
    }

    // for (int i = 0; i < genes_count; i++) {
    //     vec fvals = (tsss / resid.col(i) - ones(genes_count)) * df;
    //     for (int j = 0; j < genes_count; j++) {
    //         pvals.at(j, i) = cdf(complement(dist,
    //                              std::min(std::max(0.0, fvals[j]), 1000000.0)));
    //     }
    // }

    for (int i = 0; i < genes_count; i++) {
        for (int j = 0; j < genes_count; j++) {
            // Rcout << "Genes: " << j << " " << i << "\n";
            ev.at(j, i) = explainedVariance(X.col(i), X.col(j), slopes.at(j, i));
        }
    }


    return List::create
        (
                _["slopes"] = slopes,
                _["R2"] = r2,
                _["ev"] = ev
                //_["tss"] = tsss,
                //_["residuals"] = resid,
                //_["pvals"] = pvals
        );
}


// [[Rcpp::export]]
List pairwiseDemingRegression(const arma::mat& X) {
    int genes_count = X.n_cols;
    int sample_count = X.n_rows;

    mat resid(genes_count, genes_count);
    mat slopes(genes_count, genes_count);
    mat r2(genes_count, genes_count);
    vec tsss(genes_count);

    vec sxx(genes_count);

    for (int i = 0; i < genes_count; i++) {
        sxx[i] = dot(X.col(i), X.col(i));
    }

    for (int i = 0; i < genes_count; i++) {
        rowvec sxy = X.col(i).t() * X;


        for (int j = 0; j < genes_count; j++) {
            double sx = sxx[i];
            double sy = sxx[j];
            double xy = sxy[j];
            slopes.at(j, i) = (sy - sx +
                sqrt( (sy - sx) * (sy - sx) + 4 * xy * xy ) ) / (2 * xy);

        }


        mat residuals = X - X.col(i) * slopes.col(i).t();
        vec rss(genes_count);
        for (int j = 0; j < genes_count; j++) {
            rss[j] = dot(residuals.col(j), residuals.col(j));
        }
        resid.col(i) = rss;
        vec x_int(sample_count);
        double mean_x_intercept = mean(X.col(i));
        for (int j = 0; j < sample_count; j++) {
            x_int[j] = X.at(j, i) - mean_x_intercept;
        }
        double tss = dot(x_int, x_int);
        tsss[i] = tss;

    }

    for (int i = 0; i < genes_count; i++) {
        r2.col(i) = ones(genes_count) - resid.col(i) / tsss;
    }

    return List::create
        (
                _["slopes"] = slopes,
                _["R2"] = r2
        );

}

