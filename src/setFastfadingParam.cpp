#include <iostream>
#include <Eigen/Dense>
#include "common.h"
#include "paramFastFading.h"
#include "const.h"

using namespace Eigen;
using namespace std;

void Set_corr_matrix();
void Set_corr_dist();



void Set_fastfading_param()
{
	Set_corr_dist();
	Set_corr_matrix();
}

void Set_corr_matrix()
{
	MatrixXReal corr_matrx_LOS(7, 7);
	MatrixXReal corr_matrx_NLOS(6, 6);
	MatrixXReal corr_matrx_Out_to_In(6, 6);

	// pathloss_model == 8 for InH_ETRI channel
	if(pathloss_model < 2 || pathloss_model == 8)   // SF	K		DS		ASD		ASA		ZSD		ZSA
	{
		corr_matrx_LOS <<	1.00, 	0.50, - 0.80, - 0.40, - 0.50, 	0.20, 	0.30,
							0.50, 	1.00, - 0.50, 	0.00, 	0.00, 	0.00, 	0.10,
						  - 0.80, - 0.50, 	1.00, 	0.60, 	0.80, 	0.10, 	0.20,
						  - 0.40, 	0.00, 	0.60, 	1.00, 	0.40, 	0.50, 	0.00,
						  - 0.50, 	0.00, 	0.80, 	0.40, 	1.00, 	0.00, 	0.50,
							0.20, 	0.00, 	0.10, 	0.50, 	0.00, 	1.00, 	0.00,
							0.30, 	0.10, 	0.20, 	0.00, 	0.50, 	0.00, 	1.00;

						//	 SF		DS		ASD		ASA		ZSD		ZSA
		corr_matrx_NLOS <<	1.00, - 0.50, 	0.00, - 0.40, 	0.00, 	0.00,
						 -  0.50, 	1.00, 	0.40, 	0.00, - 0.27, - 0.06,
							0.00, 	0.40, 	1.00, 	0.00, 	0.35, 	0.23,
						  - 0.40, 	0.00, 	0.00, 	1.00, - 0.08, 	0.43,
							0.00, - 0.27, 	0.35, - 0.08, 	1.00, 	0.42,
							0.00, - 0.06, 	0.23, 	0.43, 	0.42, 	1.00;

	}
	else if(pathloss_model<4)
	{					
		corr_matrx_LOS <<	1.00, 	0.00, - 0.40, - 0.50, - 0.50, 	0.00, - 0.80,
							0.00, 	1.00, - 0.40, 	0.00, - 0.20, 	0.00, 	0.00,
						  - 0.40, - 0.40, 	1.00, 	0.40, 	0.80, - 0.20, 	0.00,
						  - 0.50, 	0.00, 	0.40, 	1.00, 	0.00, 	0.50, 	0.00,
						  - 0.50, - 0.20, 	0.80, 	0.00, 	1.00, - 0.30, 	0.40,
							0.00, 	0.00, - 0.20, 	0.50, - 0.30, 	1.00, 	0.00,
						  - 0.80, 	0.00, 	0.00, 	0.00, 	0.40, 	0.00,	1.00;
		
		corr_matrx_NLOS <<	1.00, - 0.40, - 0.60, 	0.00, 	0.00, - 0.40,
						  - 0.40,  	1.00, 	0.40, 	0.60, - 0.50, 	0.00,
						  - 0.60,  	0.40, 	1.00, 	0.40, 	0.50, - 0.10,
							0.00,  	0.60, 	0.40,	1.00, 	0.00, 	0.00,
							0.00, - 0.50, 	0.50, 	0.00, 	1.00, 	0.00,
						  - 0.40,  	0.00, - 0.10, 	0.00, 	0.00, 	1.00;

		corr_matrx_Out_to_In << 1.00, - 0.50, 	0.20, 	0.00, 	0.00, 	0.00,
							  - 0.50, 	1.00, 	0.40, 	0.40, - 0.60, - 0.20,
								0.20,  	0.40, 	1.00, 	0.00, - 0.20, 	0.00,
								0.00,  	0.40, 	0.00, 	1.00, 	0.00, 	0.50,
								0.00, - 0.60, - 0.20, 	0.00, 	1.00, 	0.50,
								0.00, - 0.20, 	0.00, 	0.50, 	0.50, 	1.00;


	}
	else if (pathloss_model<6 || pathloss_model == 9)
	{
		corr_matrx_LOS <<	1.00, 	0.50, - 0.40, - 0.50, - 0.40, 	0.00, 	0.00,
							0.50, 	1.00, - 0.70, - 0.20, - 0.30, 	0.00, 	0.00,
						  - 0.40, - 0.70, 	1.00, 	0.50, 	0.80, 	0.00, 	0.20,
						  - 0.50, - 0.20, 	0.50, 	1.00, 	0.40, 	0.50, 	0.30,
						  - 0.40, - 0.30,	0.80, 	0.40, 	1.00, 	0.00, 	0.00,
							0.00, 	0.00, 	0.00, 	0.50, 	0.00,	1.00, 	0.00,
							0.00, 	0.00, 	0.20, 	0.30, 	0.00, 	0.00, 	1.00;
	
		corr_matrx_NLOS <<	1.00, - 0.70, 	0.00, - 0.40, 	0.00, 	0.00,
						  - 0.70,  	1.00, 	0.00, 	0.40, - 0.50, 	0.00,
							0.00,  	0.00, 	1.00, 	0.00, 	0.50, 	0.50,
						  - 0.40,  	0.40,	0.00, 	1.00,	0.00, 	0.20,
							0.00, - 0.50, 	0.50, 	0.00, 	1.00, 	0.00,
							0.00, 	0.00, 	0.50, 	0.20, 	0.00, 	1.00;
	
		corr_matrx_Out_to_In << 1.00, - 0.50,	0.20, 	0.00, 	0.00, 	0.00,
							  - 0.50, 	1.00, 	0.40, 	0.40, - 0.60, - 0.20,
								0.20, 	0.40,	1.00, 	0.00, - 0.20, 	0.00,
								0.00, 	0.40, 	0.00, 	1.00,	0.00, 	0.50,
								0.00, - 0.60, - 0.20, 	0.00, 	1.00, 	0.50,
								0.00, - 0.20, 	0.00, 	0.50, 	0.50, 	1.00;

	}
	else if (pathloss_model<8)
	{
		corr_matrx_LOS <<	1.00, 	0.00, - 0.50, 	0.00, 	0.00, 	0.01, - 0.17,
							0.00, 	1.00, 	0.00, 	0.00, 	0.00, 	0.00, - 0.02,
						  - 0.50, 	0.00, 	1.00, 	0.00, 	0.00, - 0.05, 	0.27,
							0.00, 	0.00, 	0.00, 	1.00, 	0.00, 	0.73, - 0.14,
							0.00, 	0.00, 	0.00, 	0.00, 	1.00, - 0.20, 	0.24,
							0.01, 	0.00, - 0.05, 	0.73, - 0.20, 	1.00, - 0.07,
						  - 0.17, - 0.02, 	0.27, - 0.14, 	0.24, - 0.07, 	1.00;


		corr_matrx_NLOS <<	1.00, - 0.50, 	0.60, 	0.00, - 0.04, - 0.25,
						  - 0.50,  	1.00, - 0.40, 	0.00, - 0.10, - 0.40,
							0.60, - 0.40,	1.00,	0.00, 	0.42, - 0.27,
							0.00, 	0.00, 	0.00,	1.00, - 0.18, 	0.26,
						  - 0.04, - 0.10, 	0.42, - 0.18, 	1.00, - 0.27,
						  - 0.25, - 0.40, - 0.27, 	0.26, - 0.27, 	1.00;
		
		corr_matrx_Out_to_In << 1.00,  	0.00, 	0.00, 	0.00, 	0.00, 	0.00,
								0.00,  	1.00, 	0.00, 	0.00, 	0.00, 	0.00,
								0.00, 	0.00, 	1.00, - 0.70, 	0.66, 	0.47,
								0.00, 	0.00, - 0.70, 	1.00, - 0.55, - 0.22,
								0.00,	0.00, 	0.66, - 0.55, 	1.00, 	0.00,
								0.00, 	0.00, 	0.47, - 0.22, 	0.00, 	1.00;

	}
	else
	{
		cout << " check Set_corr_matrix() pathloss_model" << endl;
		getchar();
	}

	LLT<MatrixXReal> llt_LOS(corr_matrx_LOS);    /// compute the Cholesky decomposition of corr_matrx_LOS
	sqrt_corr_matrx_LOS = llt_LOS.matrixL();  /// retrieve factor L in the decomposition

	LLT<MatrixXReal> llt_NLOS(corr_matrx_NLOS);
	sqrt_corr_matrx_NLOS = llt_NLOS.matrixL();

	if (pathloss_model > 1)
	{
		LLT<MatrixXReal> llt_Out_to_In(corr_matrx_Out_to_In);
		sqrt_corr_matrx_OUT2IN = llt_Out_to_In.matrixL();
	}

}

void Set_corr_dist()
{
	// pathloss_model == 8 for InH_ETRI channel
	if (pathloss_model < 2 || pathloss_model == 8)
	{
		corr_dist_DS [LOS_propagation] = 8;
		corr_dist_ASD[LOS_propagation] = 7;
		corr_dist_ASA[LOS_propagation] = 5;
		corr_dist_SF [LOS_propagation] = 10;
		corr_dist_K  [LOS_propagation] = 4;
		corr_dist_ZSD[LOS_propagation] = 4;
		corr_dist_ZSA[LOS_propagation] = 4;

		corr_dist_DS [NLOS_propagation] = 5;
		corr_dist_ASD[NLOS_propagation] = 3;
		corr_dist_ASA[NLOS_propagation] = 3;
		corr_dist_SF [NLOS_propagation] = 6;
		corr_dist_ZSD[NLOS_propagation] = 4;
		corr_dist_ZSA[NLOS_propagation] = 4;

	}
	else if (pathloss_model < 4)
	{
		corr_dist_DS [LOS_propagation] = 30;
		corr_dist_ASD[LOS_propagation] = 18;
		corr_dist_ASA[LOS_propagation] = 15;
		corr_dist_SF [LOS_propagation] = 37;
		corr_dist_K  [LOS_propagation] = 12;
		corr_dist_ZSD[LOS_propagation] = 15;
		corr_dist_ZSA[LOS_propagation] = 15;

		corr_dist_DS [NLOS_propagation] = 40;
		corr_dist_ASD[NLOS_propagation] = 50;
		corr_dist_ASA[NLOS_propagation] = 50;
		corr_dist_SF [NLOS_propagation] = 50;
		corr_dist_ZSD[NLOS_propagation] = 50;
		corr_dist_ZSA[NLOS_propagation] = 50;

		corr_dist_DS [OUT2IN_propagation] = 10;
		corr_dist_ASD[OUT2IN_propagation] = 11;
		corr_dist_ASA[OUT2IN_propagation] = 17;
		corr_dist_SF [OUT2IN_propagation] = 7;
		corr_dist_ZSD[OUT2IN_propagation] = 25;
		corr_dist_ZSA[OUT2IN_propagation] = 25;
	}
	else if (pathloss_model<6 || pathloss_model == 9)
	{
		corr_dist_DS[LOS_propagation] = 7;
		corr_dist_ASD[LOS_propagation] = 8;
		corr_dist_ASA[LOS_propagation] = 8;
		corr_dist_SF[LOS_propagation] = 10;
		corr_dist_K[LOS_propagation] = 15;
		corr_dist_ZSD[LOS_propagation] = 12;
		corr_dist_ZSA[LOS_propagation] = 12;

		corr_dist_DS[NLOS_propagation] = 10;
		corr_dist_ASD[NLOS_propagation] = 10;
		corr_dist_ASA[NLOS_propagation] = 9;
		corr_dist_SF[NLOS_propagation] = 13;
		corr_dist_ZSD[NLOS_propagation] = 10;
		corr_dist_ZSA[NLOS_propagation] = 10;

		corr_dist_DS[OUT2IN_propagation] = 10;
		corr_dist_ASD[OUT2IN_propagation] = 11;
		corr_dist_ASA[OUT2IN_propagation] = 17;
		corr_dist_SF[OUT2IN_propagation] = 7;
		corr_dist_ZSD[OUT2IN_propagation] = 25;
		corr_dist_ZSA[OUT2IN_propagation] = 25;
	}
	else
	{
		corr_dist_DS[LOS_propagation] = 50;
		corr_dist_ASD[LOS_propagation] = 25;
		corr_dist_ASA[LOS_propagation] = 35;
		corr_dist_SF[LOS_propagation] = 37;
		corr_dist_K[LOS_propagation] = 40;
		corr_dist_ZSD[LOS_propagation] = 15;
		corr_dist_ZSA[LOS_propagation] = 15;

		corr_dist_DS[NLOS_propagation] = 36;
		corr_dist_ASD[NLOS_propagation] = 30;
		corr_dist_ASA[NLOS_propagation] = 40;
		corr_dist_SF[NLOS_propagation] = 120;
		corr_dist_ZSD[NLOS_propagation] = 50;
		corr_dist_ZSA[NLOS_propagation] = 50;

		corr_dist_DS[OUT2IN_propagation] = 36;
		corr_dist_ASD[OUT2IN_propagation] = 30;
		corr_dist_ASA[OUT2IN_propagation] = 40;
		corr_dist_SF[OUT2IN_propagation] = 120;
		corr_dist_ZSD[OUT2IN_propagation] = 50;
		corr_dist_ZSA[OUT2IN_propagation] = 50;
	}

}
