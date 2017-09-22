// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/beta.hpp>
#include "fcnnls.h"
using namespace boost::math;
using namespace Rcpp;
using namespace arma;

arma::vec cartseianFromSpherical(const arma::vec& angles);
arma::vec sphericalFromCartesian(const arma::vec& coords);
List deconvolve(const arma::mat& mixedData,
                const arma::mat& gg,
                const arma::mat& coefMatrix,
                const arma::vec& coef,
                int dims);
arma::mat getCoefMatrix(const arma::vec& par, int dims);
double score(const arma::mat& mixedData,
             const arma::mat& gg,
             const arma::vec& coef,
             const arma::vec& par,
             int dims);
double score(const arma::vec& par,
             const arma::mat& mixedData,
             const arma::mat& gg,
             const arma::vec& coef,
             int dims);
double sum_to_one(const arma::mat& mixedData,
                  const arma::mat& gg,
                  const arma::vec& coef,
                  const arma::vec& par,
                  int dims);
double sum_to_one_export(const arma::vec& par,
                         const arma::mat& mixedData,
                         const arma::mat& gg,
                         const arma::vec& coef,
                         int dims);
