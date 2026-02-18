#include "precision.h"
#include "common.h"
#include "const.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
//#include <fftw3.h>


#include <unsupported/Eigen/FFT>

using namespace std;

#include <cstring>
#include <cstdlib>

Real Get_distance(LOCATION, LOCATION);

/*===================================================================
FUNCTION: UE::Channel_Update_MIMO( void )
===================================================================*/
void MS::Channel_Update_MIMO(int ms_idx)
{
	self_idx = ms_idx;

	/*===================================================================
	FUNCTION : Declare_Ch_matrix()
	OUTPUT   : H_m, HI_m(interference channel)
	===================================================================*/
	if (CH_CAL == 1)
	{
		Declare_ch_matrix();
	}
	/*===================================================================
	FUNCTION : Fourier_Transform_of_Channel(_self_idx)
	OUTPUT   : H_m, HI_m(interference channel)
	===================================================================*/
	//Fourier_Transform_of_Channel(self_idx);
	Fourier_Transform_of_Channel_Optimized(self_idx);

	/*===================================================================
	FUNCTION : Quantization_of_Ch()
	INPUT    : H_m
	OUTPUT   : PMI(cw_idx)
	===================================================================*/
	if ( TDD_mode == 0 ) {
		if (NUM_TX_Port == 2) {
			Quantization_of_Ch_CSIRS_2();
		}
		else {
			Quantization_of_Ch();
		}
	}
	else {
	/*===================================================================
	FUNCTION : Store_CSI_for_TDD()
	INPUT    : H_m
	OUTPUT   : CSI_matrix (for TDD reciprocity-based precoding)
	===================================================================*/
		Store_CSI_for_TDD();
	}

	/*===================================================================
	FUNCTION : CQI_Update()
	INPUT    : H_m, UE_PATH_LOSS, UE_ICI
	OUTPUT   : CQI(esinr)
	===================================================================*/
	CQI_Update();
}

void MS::Quantization_of_Ch_CSIRS_2(void)
{
	if (CH_CAL == 1) //// run time calc
	{
		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		{
			int subband_size =  16;
			int num_subband  = (int)(num_rb/subband_size);

			for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
			{
				int subband_len = subband_size;
				if (subband_idx == num_subband - 1)
					subband_len = subband_size + num_rb%subband_size;

				MatrixXcReal H_sb_eq(2, 2);
				H_sb_eq << 0, 0, 0, 0;
				for(int rbs_idx = subband_idx*subband_size; rbs_idx < subband_idx*subband_size + subband_len; rbs_idx++)
				{
					JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rbs_idx], ComputeThinU | ComputeThinV);

					MatrixXcReal SVD_V;
					MatrixXcReal SVD_V_hermitian;
					MatrixXcReal SVD_S;

					MatrixXcReal P;

					VectorXcReal singular_vector = svd.singularValues();
					SVD_S = singular_vector.asDiagonal();

					SVD_V = svd.matrixV();
					SVD_V_hermitian = SVD_V.adjoint();

					P = SVD_S * SVD_V_hermitian;
					H_sb_eq += P;
				}

				VectorXcReal temp_mat;
				Real sum_of_elements = 0;
				Real max_value = 0;
				int subband_pmi_n = -1;

				for (int n = 0; n < 4; n++)
				{
					temp_mat = H_sb_eq * cb_W_csirs_2_layer_1[n];
					sum_of_elements = 0;
					for (int elem_idx = 0; elem_idx < temp_mat.rows(); elem_idx++)
						sum_of_elements += abs(temp_mat(elem_idx))*abs(temp_mat(elem_idx));

					if (sum_of_elements > max_value)
					{
						subband_pmi_n = n;
						max_value = sum_of_elements;
					}
				}

				for(int rbs_idx = subband_idx*subband_size; rbs_idx < subband_idx*subband_size + subband_len; rbs_idx++)
				{
					PMI[coeff_idx][rbs_idx][t % cqi_history_length]        = {0, 0, subband_pmi_n};
					PMI_vector[coeff_idx][rbs_idx][t % cqi_history_length] = cb_W_csirs_2_layer_1[subband_pmi_n];
				}
			}
		}
	}
}


void MS::Quantization_of_Ch(void)
{
	if (CH_CAL == 1) //// run time calc
	{
		// 220815 jhnoh
		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		{
			//cout << " coeff_idx = " << coeff_idx << endl;
			// Fixed: H_wb_eq size must match H_m size (NUM_RX_Port x NUM_TX_Port)
			MatrixXcReal H_wb_eq(NUM_RX_Port, NUM_TX_Port);
			H_wb_eq.setZero();  // Initialize to zero
			for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++)
			{
				// JacobiSVD <MatrixXcReal> svd(H_m[coeff_idx][rbs_idx], ComputeFullU | ComputeFullV);    // H = RxT, V = TxT, U = RxR
				JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rbs_idx], ComputeThinU | ComputeThinV);
				MatrixXcReal H_temp = H_m[coeff_idx][rbs_idx];

				MatrixXcReal SVD_V;
				MatrixXcReal SVD_V_hermitian;
				MatrixXcReal SVD_U;
				MatrixXcReal SVD_S;

				MatrixXcReal P;

				SVD_U = svd.matrixU();

				VectorXcReal H_eq;
				RowVectorXcReal H_eq_second;

				VectorXcReal singular_vector = svd.singularValues();
				SVD_S = singular_vector.asDiagonal();

				SVD_V = svd.matrixV();
				SVD_V_hermitian = SVD_V.adjoint();

				P = SVD_S * SVD_V_hermitian;
				H_wb_eq += P;

				H_eq_second = H_m[coeff_idx][rbs_idx].row(0);
				H_eq = SVD_V.col(0);
			}

			VectorXcReal temp_mat;
			Real sum_of_elements = 0;
			Real max_value = 0;
			pmi_l = -1;
			pmi_m = -1;
			pmi_n = -1;

			for (int l = 0; l < BS_Np * 4; l++)
			{
				for (int m = 0; m < BS_Mp * 4; m++)
				{
					for (int n = 0; n < 4; n++)
					{
						temp_mat = H_wb_eq * codebook_W[l][m][n];
						sum_of_elements = 0;
						for (int elem_idx = 0; elem_idx < temp_mat.rows(); elem_idx++)
							sum_of_elements += abs(temp_mat(elem_idx))*abs(temp_mat(elem_idx));

						if (sum_of_elements > max_value)
						{
							pmi_l = l;
							pmi_m = m;
							pmi_n = n;

							max_value = sum_of_elements; 
						}
					}
				}
			}

			int subband_size =  16;
			int num_subband  = (int)(num_rb/subband_size);

			for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
			{
				int subband_len = subband_size;
				if (subband_idx == num_subband - 1)
					subband_len = subband_size + num_rb%subband_size;

				// Fixed: H_sb_eq size must match H_m size (NUM_RX_Port x NUM_TX_Port)
				MatrixXcReal H_sb_eq(NUM_RX_Port, NUM_TX_Port);
				H_sb_eq.setZero();  // Initialize to zero

				for(int rbs_idx = subband_idx*subband_size; rbs_idx < subband_idx*subband_size + subband_len; rbs_idx++)
				{
					JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rbs_idx], ComputeThinU | ComputeThinV);

					MatrixXcReal SVD_V;
					MatrixXcReal SVD_V_hermitian;
					MatrixXcReal SVD_S;

					MatrixXcReal P;

					VectorXcReal singular_vector = svd.singularValues();
					SVD_S = singular_vector.asDiagonal();

					SVD_V = svd.matrixV();
					SVD_V_hermitian = SVD_V.adjoint();

					P = SVD_S * SVD_V_hermitian;
					H_sb_eq += P;
				}

				VectorXcReal temp_mat;
				Real sum_of_elements = 0;
				Real max_value = 0;
				int subband_pmi_n = -1;

				for (int n = 0; n < 4; n++)
				{
					temp_mat = H_sb_eq * codebook_W[pmi_l][pmi_m][n];
					sum_of_elements = 0;
					for (int elem_idx = 0; elem_idx < temp_mat.rows(); elem_idx++)
						sum_of_elements += abs(temp_mat(elem_idx))*abs(temp_mat(elem_idx));

					if (sum_of_elements > max_value)
					{
						subband_pmi_n = n;
						max_value = sum_of_elements;
					}
				}

				for(int rbs_idx = subband_idx*subband_size; rbs_idx < subband_idx*subband_size + subband_len; rbs_idx++)
				{
					PMI[coeff_idx][rbs_idx][t % cqi_history_length] = {pmi_l, pmi_m, subband_pmi_n};
					PMI_vector[coeff_idx][rbs_idx][t % cqi_history_length] = codebook_W[pmi_l][pmi_m][subband_pmi_n];
				}
			}
		}
	}
}

/*===================================================================
FUNCTION: MS::Store_CSI_for_TDD()

DESCRIPTION:
  TDD reciprocity-based approach:
  1. Store full channel matrix H_m into CSI_matrix for BS precoding
  2. OPTIMIZATION: Pre-compute dominant eigenvector (for semi-orthogonality check)
     - Calculate H^H * H eigenvector and store in PMI_vector
     - This eliminates redundant eigenvalue decomposition in scheduling

  This maintains the delay structure (using history buffer) consistent
  with FDD mode's PMI_vector feedback.

INPUT:  H_m[coeff_idx][rb_idx] - Current channel matrix (NUM_RX_Port × NUM_TX_Port)
OUTPUT: CSI_matrix[coeff_idx][rb_idx][t % cqi_history_length] - Stored channel matrix
        PMI_vector[coeff_idx][rb_idx][t % cqi_history_length] - Dominant eigenvector (TDD only)
===================================================================*/
void MS::Store_CSI_for_TDD(void)
{
	if (CH_CAL == 1) //// run-time calculation
	{
		const int subband_size = 16;  // Same as scheduling
		int num_subband = (num_rb / subband_size);
		if (num_rb % subband_size != 0) {
			num_subband += 1;  // Add one more subband for remaining RBs
		}

		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		{
			// Step 1: Store CSI_matrix for ALL RBs (per-RB granularity)
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				// Store the full channel matrix H_m (NUM_RX_Port × NUM_TX_Port)
				// This represents uplink channel measured via SRS in TDD
				// Due to reciprocity: H_uplink^T = H_downlink
				CSI_matrix[coeff_idx][rb_idx][t % cqi_history_length] = H_m[coeff_idx][rb_idx];
			}

			// Step 2: Compute PMI_vector per SUBBAND (reduced computation)
			// Use representative RB (middle of subband) to compute eigenvector
			for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
			{
				int rb_start = subband_idx * subband_size;
				int rb_end = rb_start + subband_size;
				if (rb_end > num_rb) {
					rb_end = num_rb;
				}

				// Use middle RB as representative for this subband
				int rep_rb_idx = rb_start + (rb_end - rb_start) / 2;

				// OPTIMIZATION: Pre-compute dominant eigenvector for scheduling
				// This avoids repeated eigenvalue decomposition in semi-orthogonality checks
				MatrixXcReal H = H_m[coeff_idx][rep_rb_idx];
				MatrixXcReal HH = H.adjoint() * H;  // T×T Hermitian matrix

				// Compute eigenvalue decomposition
				SelfAdjointEigenSolver<MatrixXcReal> es(HH);

				// Extract dominant eigenvector (corresponding to largest eigenvalue)
				VectorXcReal eigvec = es.eigenvectors().col(es.eigenvectors().cols() - 1);

				// Normalize
				Real norm = eigvec.norm();
				VectorXcReal normalized_eigvec;
				if (norm > 1e-12) {
					normalized_eigvec = eigvec / norm;
				} else {
					// Degenerate case: use zero vector
					normalized_eigvec = VectorXcReal::Zero(NUM_TX_Port);
				}

				// Store same eigenvector for all RBs in this subband
				for (int rb_idx = rb_start; rb_idx < rb_end; rb_idx++)
				{
					PMI_vector[coeff_idx][rb_idx][t % cqi_history_length] = normalized_eigvec;
				}
			}
		}
	}
}

void MS::CQI_Update(void)
{
	if (CH_CAL == 1) //// run-time calculate
	{
		// if (NUM_TX != 4)   // No Use PMI
		if (Scheduling_Type == 0) // round robin scheduling
		{
			// Optimization 1: Cache dBm2linear conversions (moved outside RB loop)
			Real linear_signal_cached = dBm2linear(links[self_idx].str_signal);
			Real linear_interference_cached = dBm2linear(links[self_idx].interference);

			// Optimization 5: Use cached noise value (computed once in MS::Configuration)
			Real total_interference_noise = linear_interference_cached + noise_linear_cached;

			// Pre-compute SNR ratio (same for all RBs)
			Real snr_ratio = linear_signal_cached / total_interference_noise;

			for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			{
				for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
				{
					// Optimization 2: Use Eigen's squaredNorm() instead of manual loops
					// This enables SIMD vectorization (AVX/SSE) automatically
					Real channel_power = H_m[coeff_idx][rb_idx].squaredNorm();

					CQI[coeff_idx][rb_idx][abs(t) % cqi_history_length] = channel_power * snr_ratio;
				}
			}
		}
		else if (Scheduling_Type == 1) //|| Scheduling_Type == 2)  //MU MIMO scheduling..
		{
			// Subband-based CQI computation (same granularity as scheduling)
			const int subband_size = 16;
			int num_subband = (num_rb / subband_size);
			if (num_rb % subband_size != 0) {
				num_subband += 1;
			}

			// Optimization 1: Cache dBm2linear conversions
			Real linear_signal_cached = dBm2linear(links[self_idx].str_signal);
			Real linear_interference_cached = dBm2linear(links[self_idx].interference);

			// Optimization 5: Use cached noise value
			Real total_interference_noise = linear_interference_cached + noise_linear_cached;

			// Optimization 1: Pre-compute sqrt(linear_signal) once
			Real linear_signal_sqrt = sqrt(linear_signal_cached);

			// Pre-allocate identity matrices (constant across all subbands)
			MatrixXcReal eye_R = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
			MatrixXcReal eye_T = MatrixXcReal::Identity(NUM_TX_Port, NUM_TX_Port);
			MatrixXcReal C = total_interference_noise * eye_R;  // Interference covariance (constant)

			for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			{
				// Compute CQI per subband
				for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
				{
					int rb_start = subband_idx * subband_size;
					int rb_end = rb_start + subband_size;
					if (rb_end > num_rb) {
						rb_end = num_rb;
					}

					// Use middle RB as representative for this subband
					int rep_rb_idx = rb_start + (rb_end - rb_start) / 2;

					// Compute CQI using representative RB
					MatrixXcReal H_eq = linear_signal_sqrt * H_m[coeff_idx][rep_rb_idx]; // R x T

					VectorXcReal w = PMI_vector[coeff_idx][rep_rb_idx][t % cqi_history_length].conjugate(); // Tx1

					// Signal covariance matrix
					MatrixXcReal A = H_eq * w * w.adjoint() * H_eq.adjoint();  // RxR

					// Inter-stream interference covariance matrix
					MatrixXcReal B = (mx_ue_mumimo-1)/((Real)NUM_TX_Port)
					                *(H_eq * (eye_T - w * w.adjoint()) * H_eq.adjoint()); // RxR

					MatrixXcReal BC = B + C;  // Total interference + noise

					// Optimization 4: Use LU decomposition instead of explicit inverse
					MatrixXcReal M = BC.lu().solve(A);

					// Compute dominant eigenvector using selected method
					VectorXcReal dominant_eigvec;
					if (USE_POWER_ITERATION == 1) {
						// Power Iteration method (faster but approximate)
						dominant_eigvec = PowerIteration_DominantEigenvector(M, POWER_ITERATION_MAX_ITER);
					} else {
						// SelfAdjointEigenSolver (more accurate, default)
						SelfAdjointEigenSolver<MatrixXcReal> es(M);
						dominant_eigvec = es.eigenvectors().col(es.eigenvectors().cols() - 1);
					}
					RowVectorXcReal u = dominant_eigvec.adjoint();

					// Compute SINR with numerical stability check
					Real numerator = (u * A * u.adjoint()).norm();
					Real denominator = (u * BC * u.adjoint()).norm();

					Real sinr_estimate = (denominator > 1e-12) ? (numerator / denominator) : 0.0;

					// Store same CQI for all RBs in this subband
					for (int rb_idx = rb_start; rb_idx < rb_end; rb_idx++)
					{
						CQI[coeff_idx][rb_idx][abs(t) % cqi_history_length] = sinr_estimate;
					}
				}
			}
		}
		else if (Scheduling_Type == 2)
		{
			MatrixXcReal w, u, h_u;
			MatrixXcReal H_bar;
			MatrixXcReal H_inf;
			MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
			Real A, B, C;

			for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			{
				Real temp_sinr = 0;
				for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
				{
					if (coeff_idx == 0)
					{
						int sec_idx  = links[self_idx].adj_sector[0];
						int bs_idx   = (int)(sec_idx / 3.);
						Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);

						// jhnoh 220901
						Real linear_signal       = dBm2linear(links[self_idx].str_signal); // dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
						Real linear_interference = dBm2linear(links[self_idx].interference);

						H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length];
						u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
						A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
						C = (linear_interference + noise) * (u * u.adjoint()).norm();

						temp_sinr = A / C;
						if ( comp_mode[self_idx] == 1)
						{
							Real comp_sector_static_gain = dBm2linear(links[self_idx].comp_interf_strength[1]);						
							int comp_sector_idx  = links[self_idx].comp_sector_idx;
							int com_bs_idx = (int) (comp_sector_idx/3);
							if ( comp_sector_idx < num_SECTORS)
								linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);
							else
							{
								linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*com_bs_idx]);
								linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*com_bs_idx+1]);
								linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*com_bs_idx+2]);
							}
							linear_interference += comp_sector_static_gain;
							u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
							A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
							C = (linear_interference + noise) * (u * u.adjoint()).norm();

							CQI_comp[coeff_idx][rb_idx][abs(t) % cqi_history_length] = A/C;
						}
					}
					else if (coeff_idx == 1)
					{
						int sec_idx = links[self_idx].comp_sector_idx;
						int bs_idx = (int)(sec_idx / 3.);
						Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);

						// jhnoh 220901
						Real linear_signal = dBm2linear(links[self_idx].comp_signal_strength[1]);

						Real linear_interference = dBm2linear(links[self_idx].interference);
						if (sec_idx < num_SECTORS)
							linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[sec_idx]);
						else
						{
							linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx]);
							linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx+1]);
							linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx+2]);
						}

						linear_interference += dBm2linear(links[self_idx].comp_interf_strength[0]);

						H_bar = sqrt(linear_signal) * (H_m[coeff_idx][rb_idx]) * PMI_vector[coeff_idx][rb_idx][t % cqi_history_length];

						int main_sector = (int)links[self_idx].adj_sector[0] / 3.;
						u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
						A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
						C = (linear_interference + noise) * (u * u.adjoint()).norm();

						temp_sinr = A / C;
					}
					Real sinr = temp_sinr;
					CQI[coeff_idx][rb_idx][abs(t) % cqi_history_length] = sinr;
				}
			}
		}
		else if (Scheduling_Type == 3)
		{
			CQI_Update_CJT();
		}
	}
}

void MS::CQI_Update_CJT(void)
{
	for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
	{
		MatrixXcReal w, u, h_u;
		MatrixXcReal H_bar1;
		MatrixXcReal H_bar2;
		MatrixXcReal H_bar;
		MatrixXcReal H_inf;

		MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
		Real A, B, C;	

		int sec_idx = links[self_idx].adj_sector[0];
		int bs_idx = (int)(sec_idx / 3.);
		Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
		int comp_sector_idx = links[self_idx].adj_sector[1];
		int comp_bs_idx = (int)comp_sector_idx / 3.;

		Real linear_signal       = dBm2linear(links[self_idx].str_signal); // dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
		Real linear_interference = dBm2linear(links[self_idx].interference);
		linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);

		H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length];
		Real comp_sector_static_gain = dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]); //dBm2linear(links[self_idx].static_gain[1].first); // dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
		H_inf = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length];

		u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + H_inf * H_inf.adjoint() + (linear_interference + noise) * Identity).inverse();
		A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
		B = (u * H_inf * H_inf.adjoint() * u.adjoint()).norm();
		C = (linear_interference + noise) * (u * u.adjoint()).norm();		

		CQI[0][rb_idx][abs(t) % cqi_history_length] = A/(B+C);

		// jhnoh 220901
		if ( comp_mode[self_idx] == 1 && g_comp_mode )
		{
			Real linear_signal       = dBm2linear(links[self_idx].str_signal); 
			Real linear_interference = dBm2linear(links[self_idx].interference);
			linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);

			H_bar1 = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length];

			Real comp_sector_static_gain = dBm2linear(links[self_idx].static_gain[1].first); 
			H_bar2 = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length];
			H_bar = H_bar1 + H_bar2;

			u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
			A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
			B = 0;
			C = (linear_interference + noise) * (u * u.adjoint()).norm();

			CQI[1][rb_idx][abs(t) % cqi_history_length] = A/(B+C);
			CQI_comp[0][rb_idx][abs(t) % cqi_history_length] = A/(B+C);
		}		
	}
}

void MS::Fourier_Transform_of_Channel(int ms_idx)
{
	if (CH_CAL == 1) //// runtime calc
	{
		int self_idx = ms_idx;
		int _bs_idx = 0;

		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		{
			if (TYPE == 11 && num_Indoor_TRxP == 1)
			{
				//_bs_idx = links[self_idx]._sector_in_control;
				_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx]);
			}
			else
			{
				_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx] / 3);
			}

			for (int r_idx = 0; r_idx < NUM_RX_Port; r_idx++) // Received    Channel
			{
				for (int t_idx = 0; t_idx < NUM_TX_Port; t_idx++) // Transmitted Channel
				{
					int pointer = (int(fft_size) - (num_rb * num_freq_per_rbs)) / 2; ////  (1024 - 50*12) / 2 = 212
					Real sum_for_normal = 0;

					//const int nbin = 4096;

					// Real 타입에 맞게 FFT 수행
					//ComplexReal timevec[4096];
					//ComplexReal freqvec[4096];

					//for (int iq = 0; iq < nbin; iq++) {
					//	timevec[iq] = channel[_bs_idx][self_idx].CHIR_vec[(int)(links[self_idx].adj_sector[coeff_idx] % 3)][t_idx][r_idx][iq];
					//}

					//fft.fwd( freqvec, timevec, nbin);

					// Eigen FFT 결과를 H_m에 저장
					//for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++) // Subcarrier idx
					//{
					//	int sbc_idx = pointer + (num_freq_per_rbs / 2) + (rbs_idx * num_freq_per_rbs);
					//	H_m[coeff_idx][rbs_idx](r_idx, t_idx) = freqvec[sbc_idx];
					//}

					/*
					fftw_complex *array_fft_ini;
					fftw_complex *array_fft_fin;

					array_fft_ini = (fftw_complex *)fftw_malloc(nbin * sizeof(fftw_complex));
					array_fft_fin = (fftw_complex *)fftw_malloc(nbin * sizeof(fftw_complex));

					fftw_plan fft_plan; 
					fft_plan = fftw_plan_dft_1d(nbin, array_fft_ini, array_fft_fin,FFTW_FORWARD, FFTW_ESTIMATE);

					for (int iq = 0; iq < nbin; iq++) {
						array_fft_ini[iq][0] = channel[_bs_idx][self_idx].CHIR_vec[(int)(links[self_idx].adj_sector[coeff_idx] % 3)][t_idx][r_idx][iq].real();
						array_fft_ini[iq][1] = channel[_bs_idx][self_idx].CHIR_vec[(int)(links[self_idx].adj_sector[coeff_idx] % 3)][t_idx][r_idx][iq].imag();
					}

					fftw_execute(fft_plan);
					
					for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++) // Subcarrier idx
					{
						int sbc_idx = pointer + (num_freq_per_rbs / 2) + (rbs_idx * num_freq_per_rbs);
						H_m[coeff_idx][rbs_idx](r_idx, t_idx) = ComplexReal (array_fft_fin[sbc_idx][0],array_fft_fin[sbc_idx][1]);
					}
					*/

					
					for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++) // Subcarrier idx
					{
						unsigned long long freq_idx = (pointer + (num_freq_per_rbs / 2) + (rbs_idx * num_freq_per_rbs)) * subcarrier_spacing; //// 15000 = 15kHz subcarrier   ///// (212 + (12/2) + (rbs_idx * 12)) * 15k
						for (int path_idx = 0; path_idx < channel[_bs_idx][self_idx].NUM_PATH_for_channelcoeff; path_idx++)
						{
							complex<Real> dft_exp(0.0, -2 * pi * freq_idx * (channel[_bs_idx][self_idx].delay[path_idx]));
							H_m[coeff_idx][rbs_idx](r_idx, t_idx) += channel[_bs_idx][self_idx].CHIR[(int)(links[self_idx].adj_sector[coeff_idx] % 3)][t_idx][r_idx][path_idx] * exp(dft_exp);
						}
						//sum_for_normal += pow(abs(H_m[coeff_idx][rbs_idx](r_idx, t_idx)), 2.);
					}
					

					//sum_for_normal = sum_for_normal / (Real)num_rb;
					//for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++) // Subcarrier idx
					//{
					//	H_m[coeff_idx][rbs_idx](r_idx, t_idx) = H_m[coeff_idx][rbs_idx](r_idx, t_idx) / sqrt((Real)sum_for_normal);
					//}
				}
			}
		}
	}
}

// ========================================================================================
// OPTIMIZED VERSION: Fourier_Transform_of_Channel_Optimized
// ========================================================================================
// 주요 최적화 사항:
// 1. exp() 호출 횟수 감소: O(num_rb * NUM_PATH) -> O(NUM_PATH)
// 2. 재귀적 phasor 업데이트: 복소곱만으로 다음 RB 계산
// 3. 예상 성능 향상: 10-30배 (num_rb=275, NUM_PATH=20 기준 약 25배)
//
// 원리:
//   주파수 그리드: f_r = f_0 + r * Δf_rb (r = RB index)
//   위상 인자: exp(-j2π f_r τ) = exp(-j2π f_0 τ) * [exp(-j2π Δf_rb τ)]^r
//   재귀 관계: phasor[r+1] = phasor[r] * step (step = exp(-j2π Δf_rb τ))
// ========================================================================================
void MS::Fourier_Transform_of_Channel_Optimized(int ms_idx)
{
	if (CH_CAL == 1) //// runtime calc
	{
		int self_idx = ms_idx;
		int _bs_idx = 0;

		// Debug output for MS 0, t=0 to 10
		static bool debug_output_done = false;
		bool should_output = (ms_idx == 0 && t >= 0 && t <= 10 && !debug_output_done);
		std::ofstream impulse_file, frequency_file;

		if (should_output && t == 10) {
			debug_output_done = true; // Only output once for all time steps
		}

		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		{
			if (TYPE == 11 && num_Indoor_TRxP == 1)
			{
				_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx]);
			}
			else
			{
				_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx] / 3);
			}

			for (int r_idx = 0; r_idx < NUM_RX_Port; r_idx++) // Received Channel
			{
				for (int t_idx = 0; t_idx < NUM_TX_Port; t_idx++) // Transmitted Channel
				{
					// FFT 파라미터 설정
					int pointer = (int(fft_size) - (num_rb * num_freq_per_rbs)) / 2;

					const int NUM_PATH = channel[_bs_idx][self_idx].NUM_PATH_for_channelcoeff;
					const int sec_idx = (int)(links[self_idx].adj_sector[coeff_idx] % 3);

					// RB 중앙 주파수 계산 (첫 번째 RB)
					const unsigned long long f0 = (pointer + (num_freq_per_rbs / 2)) * subcarrier_spacing;
					const Real df   = subcarrier_spacing;          // 30 kHz (numerology 1)
					const Real dfrb = num_freq_per_rbs * df;       // 360 kHz (12 subcarriers per RB)

					// Path별 phasor 초기화 및 스텝 계산 (exp 호출은 여기서만!)
					std::vector<complex<Real>> phasor(NUM_PATH);
					std::vector<complex<Real>> step(NUM_PATH);

					for (int path_idx = 0; path_idx < NUM_PATH; path_idx++)
					{
						const Real tau = channel[_bs_idx][self_idx].delay[path_idx];
						// 기준 RB(r=0) 중앙 주파수에서의 phase factor
						phasor[path_idx] = exp(complex<Real>(0.0, -2.0 * pi * f0 * tau));
						// RB 간 스텝 (다음 RB로 이동시 곱할 값)
						step[path_idx]   = exp(complex<Real>(0.0, -2.0 * pi * dfrb * tau));
					}

					// RB 루프: 재귀적으로 phasor 업데이트하며 채널 계수 계산
					for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++)
					{
						complex<Real> Hc(0.0, 0.0);

						// 모든 path 기여도 합산
						for (int path_idx = 0; path_idx < NUM_PATH; path_idx++)
						{
							const complex<Real>& gain = channel[_bs_idx][self_idx].CHIR[sec_idx][t_idx][r_idx][path_idx];
							Hc += gain * phasor[path_idx];      // 복소곱 + 복소합
							phasor[path_idx] *= step[path_idx]; // 다음 RB를 위한 누적 (복소곱 1회)

							// NaN 검출: gain 또는 phasor에서 NaN 발생 시 즉시 출력
							if (std::isnan(Hc.real()) || std::isnan(Hc.imag())) {
								static int nan_count = 0;
								if (nan_count < 5) {
									nan_count++;
									cout << "\n*** NaN DETECTED in Fourier Transform ***" << endl;
									cout << "  ms_idx=" << ms_idx << " t=" << t
									     << " coeff=" << coeff_idx << " rb=" << rbs_idx
									     << " rx=" << r_idx << " tx=" << t_idx
									     << " path=" << path_idx << endl;
									cout << "  bs_idx=" << _bs_idx << " sec_idx=" << sec_idx
									     << " Prop=" << channel[_bs_idx][self_idx].Propagation
									     << " LOS=" << channel[_bs_idx][self_idx].LOS
									     << " K=" << channel[_bs_idx][self_idx].K_linear << endl;
									cout << "  gain=(" << gain.real() << "," << gain.imag() << ")"
									     << " phasor=(" << phasor[path_idx].real() << "," << phasor[path_idx].imag() << ")"
									     << " step=(" << step[path_idx].real() << "," << step[path_idx].imag() << ")" << endl;
									cout << "  delay[" << path_idx << "]=" << channel[_bs_idx][self_idx].delay[path_idx]
									     << " power[" << path_idx << "]=" << channel[_bs_idx][self_idx].power[path_idx] << endl;
									cout << "  Hc so far=(" << Hc.real() << "," << Hc.imag() << ")" << endl;
									// 전체 CHIR 덤프
									for (int p = 0; p < NUM_PATH; p++) {
										const complex<Real>& g = channel[_bs_idx][self_idx].CHIR[sec_idx][t_idx][r_idx][p];
										cout << "    CHIR[" << p << "]=(" << g.real() << "," << g.imag()
										     << ") delay=" << channel[_bs_idx][self_idx].delay[p]
										     << " power=" << channel[_bs_idx][self_idx].power[p] << endl;
									}
								}
								break; // 이 RB의 나머지 path 스킵
							}
						}

						H_m[coeff_idx][rbs_idx](r_idx, t_idx) = Hc;
					}

					// Output impulse response and frequency response for debugging
					/*
					if (should_output && coeff_idx == 0 && r_idx == 0 && t_idx == 0)
					{
						// Open files on first call
						char impulse_filename[256], frequency_filename[256];
						sprintf(impulse_filename, "channel_impulse_response_t%d.txt", t);
						sprintf(frequency_filename, "channel_frequency_response_t%d.txt", t);

						impulse_file.open(impulse_filename);
						frequency_file.open(frequency_filename);

						// Write impulse response (CHIR)
						impulse_file << "# Time: " << t << ", MS: " << ms_idx << ", RX: " << r_idx << ", TX: " << t_idx << "\n";
						impulse_file << "# path_idx delay(s) real_gain imag_gain abs_gain phase(rad)\n";
						for (int path_idx = 0; path_idx < NUM_PATH; path_idx++)
						{
							const Real tau = channel[_bs_idx][self_idx].delay[path_idx];
							const complex<Real>& gain = channel[_bs_idx][self_idx].CHIR[sec_idx][t_idx][r_idx][path_idx];
							impulse_file << path_idx << " "
								<< std::scientific << std::setprecision(12) << tau << " "
								<< gain.real() << " " << gain.imag() << " "
								<< std::abs(gain) << " " << std::arg(gain) << "\n";
						}
						impulse_file.close();

						// Write frequency response (H_m at RBs)
						frequency_file << "# Time: " << t << ", MS: " << ms_idx << ", RX: " << r_idx << ", TX: " << t_idx << "\n";
						frequency_file << "# rb_idx freq(Hz) real_H imag_H abs_H phase(rad)\n";
						for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++)
						{
							const unsigned long long freq = f0 + rbs_idx * dfrb;
							const complex<Real>& Hf = H_m[coeff_idx][rbs_idx](r_idx, t_idx);
							frequency_file << rbs_idx << " " << freq << " "
								<< std::scientific << std::setprecision(12)
								<< Hf.real() << " " << Hf.imag() << " "
								<< std::abs(Hf) << " " << std::arg(Hf) << "\n";
						}
						frequency_file.close();
					}
					*/
				}
			}
		}
	}
}

void MS::Declare_ch_matrix(void)
{
	// 220815 jhnoh
	if (H_m == NULL)
	{
		H_m = new MatrixXcReal *[num_compute_coef];
		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			H_m[coeff_idx] = new MatrixXcReal[num_rb];
	}

	// 220815 jhnoh
	for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			MatrixXcReal temp(NUM_RX_Port, NUM_TX_Port);

			for (int r_idx = 0; r_idx < NUM_RX_Port; r_idx++)
			{
				for (int t_idx = 0; t_idx < NUM_TX_Port; t_idx++)
				{
					temp(r_idx, t_idx) = 0.;
				}
			}
			H_m[coeff_idx][freq_idx] = temp;
		}
}

/*===================================================================
FUNCTION: MS::PowerIteration_DominantEigenvector
PURPOSE : Compute dominant eigenvector using Power Iteration method
INPUT   : M - Hermitian matrix
          max_iterations - number of iterations (default: 5)
OUTPUT  : Dominant eigenvector (normalized)
NOTES   : Alternative to SelfAdjointEigenSolver for speed comparison
===================================================================*/
VectorXcReal MS::PowerIteration_DominantEigenvector(const MatrixXcReal& M, int max_iterations)
{
	int n = M.rows();

	// Initialize with random vector
	VectorXcReal v = VectorXcReal::Random(n);

	// Normalize initial vector
	Real norm = v.norm();
	if (norm > 1e-12) {
		v = v / norm;
	} else {
		// Fallback to unit vector if random is too small
		v = VectorXcReal::Zero(n);
		v(0) = ComplexReal(1.0, 0.0);
	}

	// Power iteration
	for (int iter = 0; iter < max_iterations; iter++) {
		v = M * v;  // Matrix-vector multiplication

		Real norm = v.norm();
		if (norm > 1e-12) {
			v = v / norm;  // Normalize
		} else {
			// Matrix is nearly zero, return current vector
			break;
		}
	}

	return v;
}


// ====================================================================
// Element-Level Fourier Transform (ns-3 style)
// ====================================================================
// Converts H_usn[cluster](totalRx, totalTx) in delay domain
// to H_m_elem[rb](totalRx, totalTx) in frequency domain.
// Uses the same recursive phasor technique as Fourier_Transform_of_Channel_Optimized.
// ====================================================================

void MS::Fourier_Transform_ElementLevel(int ms_idx)
{
	if (CH_CAL != 1) return;

	int self_idx = ms_idx;
	int _bs_idx = 0;

	// Only process the serving BS (coeff_idx=0)
	int coeff_idx = 0;
	if (TYPE == 11 && num_Indoor_TRxP == 1)
		_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx]);
	else
		_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx] / 3);

	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;
	const int NUM_PATH = channel[_bs_idx][self_idx].NUM_PATH_for_channelcoeff;

	if (NUM_PATH <= 0) return;
	if ((int)channel[_bs_idx][self_idx].H_usn.size() < NUM_PATH) return;

	// Allocate H_m_elem if needed
	if (H_m_elem == NULL) {
		H_m_elem = new MatrixXcReal[num_rb];
		for (int rb = 0; rb < num_rb; rb++) {
			H_m_elem[rb] = MatrixXcReal::Zero(totalRx, totalTx);
		}
	}

	// FFT parameters (same as Fourier_Transform_of_Channel_Optimized)
	int pointer = (int(fft_size) - (num_rb * num_freq_per_rbs)) / 2;
	const unsigned long long f0 = (pointer + (num_freq_per_rbs / 2)) * subcarrier_spacing;
	const Real df   = subcarrier_spacing;
	const Real dfrb = num_freq_per_rbs * df;

	// Initialize phasor and step per path
	std::vector<complex<Real>> phasor(NUM_PATH);
	std::vector<complex<Real>> step_val(NUM_PATH);

	for (int path_idx = 0; path_idx < NUM_PATH; path_idx++) {
		const Real tau = channel[_bs_idx][self_idx].delay[path_idx];
		phasor[path_idx] = exp(complex<Real>(0.0, -2.0 * pi * f0 * tau));
		step_val[path_idx] = exp(complex<Real>(0.0, -2.0 * pi * dfrb * tau));
	}

	// RB loop with recursive phasor updates
	for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++) {
		H_m_elem[rbs_idx].setZero();

		for (int path_idx = 0; path_idx < NUM_PATH; path_idx++) {
			// H_m_elem[rb] += H_usn[path] * phasor[path]
			H_m_elem[rbs_idx] += channel[_bs_idx][self_idx].H_usn[path_idx] * phasor[path_idx];

			// Advance phasor for next RB
			phasor[path_idx] *= step_val[path_idx];
		}
	}

	// Reset phasors for potential re-use
	// (not strictly needed since this is called once per drop)
}

