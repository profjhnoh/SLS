#ifndef _PARAMFASTFADING_H_
#define _PARAMFASTFADING_H_

#include <Eigen/Dense>
using namespace Eigen;

extern MatrixXReal sqrt_corr_matrx_LOS;
extern MatrixXReal sqrt_corr_matrx_NLOS;
extern MatrixXReal sqrt_corr_matrx_OUT2IN;
extern Real corr_dist_DS[3];
extern Real corr_dist_ASD[3];
extern Real corr_dist_ASA[3];
extern Real corr_dist_SF[3];
extern Real corr_dist_K[3];
extern Real corr_dist_ZSD[3];
extern Real corr_dist_ZSA[3];


#endif