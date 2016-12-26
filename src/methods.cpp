// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <omp.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include "fcnnls.h"
#include "utils.h"

using namespace boost::math;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::vec analysis_c(const arma::mat& mixedData,
                     const arma::mat& gg,
                     const arma::vec& coef,
                     int dims,
                     int iterations=1000 * 1000) {

    int len = dims * (dims - 1);
    vec lowerBound(len, fill::zeros);
    vec higherBound(len); higherBound.fill(datum::pi / 2);
    vec start = randu(len) * (datum::pi / 2);

    std::ostream nullstream(0);
    set_stream_err2(nullstream);

    vec current_point = start;
    double current_score = score(mixedData, gg, coef, start, dims);

    for (int i = 0; i < iterations; i++) {
        vec next_point = randu(len) * (datum::pi / 2);
        double next_score = score(mixedData, gg, coef, next_point, dims);
        if (next_score < current_score) {
            Rcout << "Found better value: " << next_score << "\n";
            current_point = next_point;
            current_score = next_score;
        }
    }

    return current_point;
}

// [[Rcpp::export]]
arma::vec local_minimum(const arma::mat& mixedData,
                        const arma::mat& gg,
                        const arma::vec& coef,
                        const arma::vec& initial_solution,
                        int dims,
                        double local_range=0.1,
                        int local_step=100,
                        int iterations=100) {
    int len = dims * (dims - 1);
    std::ostream nullstream(0);
    set_stream_err2(nullstream);

    vec current_point = initial_solution;
    double current_score = score(mixedData, gg, coef, current_point, dims);

    vec local_point = current_point;
    double local_score = current_score;
    Rcout << "Initial solution: " << local_score << "\n";
    for (int i = 0; i < iterations; i++) {
        for (int j = 0; j < local_step; j++) {
            vec next_point = current_point + (randu(len) - 0.5) * local_range;
            next_point.elem( find(next_point < 0) ).zeros();
            next_point.elem( find(next_point > datum::pi / 2) ).fill(datum::pi/2);
            double next_score = score(mixedData, gg, coef, next_point, dims);

            if (next_score < local_score) {
                local_score = next_score;
                local_point = next_point;
            }

        }

        if (local_score < current_score) {
            Rcout << "Found better solution: " << local_score << "\n";
            current_point = local_point;
            current_score = local_score;
        }

    }
    return current_point;

}


// [[Rcpp::export]]
arma::vec particle_swarm_optimisation(const arma::mat& mixedData,
                                      const arma::mat& gg,
                                      const arma::vec& coef,
                                      int dims,
                                      int particles = 50,
                                      int fun_calls = 100000,
                                      double omega = 1,
                                      double phip = 0.05,
                                      double phig = 0.05,
                                      int threads_num=2) {
    // supress warnings
    std::ostream nullstream(0);
    set_stream_err2(nullstream);
    omp_set_num_threads(threads_num);

    int len = dims * (dims - 1);
    vec lowerBound(len, fill::zeros);
    vec higherBound(len); higherBound.fill(datum::pi / 2);

    mat swarm(len, particles);
    mat ps(len, particles);
    vec pscores(particles);
    mat vs(len, particles);
    vec global_best(len);
    double global_best_score = datum::inf;
    for (int i = 0; i < particles; i++) {
        vec sp = randu(len) / 2.0 + 0.5;
        swarm.col(i) = sp;
        ps.col(i) = sp;
        double cur_score = score(mixedData, gg, coef, sp, dims);
        pscores[i] = cur_score;
        if (cur_score < global_best_score) {
            global_best_score = cur_score;
            global_best = sp;
        }
        vs.col(i) = (randu(len) - 0.5);
    }

    Rcout << "Current global solution: " << global_best_score << "\n";

    int iter = 0;
    mat tmp(len, 3);

    for (int iter = 0; iter < fun_calls; iter += particles) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < particles; i++) {
            mat tmp(len, 3);
            tmp.col(0) = vs.col(i); tmp.col(1) = ps.col(i) - swarm.col(i); tmp.col(2) = global_best - swarm.col(i);
            vec multi(3); multi[0] = omega; multi[1] = phip * randu(1)[0]; multi[2] = phig * randu(1)[0];
            vec new_point = swarm.col(i) + tmp * multi;
            new_point.elem( find(new_point > 1) ).ones();
            new_point.elem( find(new_point < 0.5) ).fill(0.5);
            swarm.col(i) = new_point;
            double cur_score = score(mixedData, gg, coef, swarm.col(i), dims);
            if (cur_score < pscores[i]) {
                pscores[i] = cur_score;
                ps.col(i) = new_point;
                #pragma omp critical
                if (cur_score < global_best_score) {
                    global_best_score = cur_score;
                    global_best = new_point;
                    {
                        Rcout << "Curren global solution " << global_best_score << "\n";
                    }
                }
            }

        }
    }

    return global_best;
}


// // [[Rcpp::export]]
// arma::vec simulated_annealing(const arma::mat& mixedData,
//                               const arma::mat& gg,
//                               const arma::vec& coef,
//                               int dims,
//                               //
//                               int threads_num=1) {
//
//     // basic guideline: initial point
//     int len = dims * (dims - 1);
//     vec current = randu(len) * datum::pi / 2;
//     double current_score = score(mixedData, gg, coef, current, dims);
//     return vec(1);
// }

