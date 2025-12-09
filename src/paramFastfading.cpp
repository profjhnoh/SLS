#include <Eigen/Dense>
#include "precision.h"
using namespace Eigen;

MatrixXReal sqrt_corr_matrx_LOS(7, 7);
MatrixXReal sqrt_corr_matrx_NLOS(7, 7);
MatrixXReal sqrt_corr_matrx_OUT2IN(7, 7);

Real corr_dist_DS[3];
Real corr_dist_ASD[3];
Real corr_dist_ASA[3];
Real corr_dist_SF[3];
Real corr_dist_K[3];
Real corr_dist_ZSD[3];
Real corr_dist_ZSA[3];