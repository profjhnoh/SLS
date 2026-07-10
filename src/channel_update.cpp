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
	//Fourier_Transform_of_Channel_Optimized(self_idx);
	Fourier_Transform_WithBF(self_idx);

	/*===================================================================
	FUNCTION : Compute_RI()
	INPUT    : H_m
	OUTPUT   : self_RI, su_capacity
	Only meaningful for eType II with rank_adaptive enabled.
	Otherwise self_RI is fixed to min(g_type2_rank, 1 when Codebook_Type<3).
	===================================================================*/
	Compute_RI();

	/*===================================================================
	FUNCTION : Quantization_of_Ch()
	INPUT    : H_m
	OUTPUT   : PMI(cw_idx)
	===================================================================*/
	if ( TDD_mode == 0 ) {
		if (NUM_TX_Port == 2) {
			Quantization_of_Ch_CSIRS_2();
		}
		else if (g_codebook_type == 3) {
			Quantization_of_Ch_EType2();
		}
		else if (g_codebook_type == 2) {
			Quantization_of_Ch_Type2();
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

/*===================================================================
FUNCTION: MS::Compute_RI()

DESCRIPTION:
  Rank Indicator (RI) selection at UE side.

  When g_rank_adaptive == 1 (and g_codebook_type == 3 = eType II):
    - Sparsely samples RBs and SVD-decomposes H_m[0][rb].
    - For each candidate rank R ∈ {1 .. R_max = g_type2_rank}, computes
      Shannon sum capacity assuming equal-power allocation:
          C_SU(R) = Σ_{ℓ=1..R} log₂(1 + σ_ℓ² · ρ / R)
      where ρ = linear_signal / (linear_interference + noise).
    - Picks argmax_R C_SU(R) as self_RI, stores best capacity in su_capacity.

  Otherwise:
    - self_RI = 1 (for Type I / Type II).
    - self_RI = g_type2_rank (for fixed-rank eType II).

  Rationale:
    - At low SNR, σ_2, σ_3, ... are buried in noise → rank 1 is optimal.
    - At high SNR with rich scattering, rank 2+ gains capacity from
      parallel spatial streams.
    - The SVD Shannon metric captures both effects naturally.

  Called once per time slot, before Quantization_of_Ch_EType2.
  Sparse sampling (every 10th RB) keeps cost ~num_rb/10 extra SVDs per UE.
===================================================================*/
void MS::Compute_RI(void)
{
	if (CH_CAL != 1) return;

	// Default rank
	if (g_codebook_type != 3) {
		self_RI = 1;
		su_capacity = 0.0;
		return;
	}
	if (g_rank_adaptive == 0) {
		self_RI = g_type2_rank;
		su_capacity = 0.0;
		return;
	}

	int R_max = g_type2_rank;
	if (R_max > NUM_RX_Port) R_max = NUM_RX_Port;
	if (R_max < 1) R_max = 1;

	// Total available SNR (linear). Use serving BS pathloss + interference + noise.
	// NOTE: signal = RSRP w/ fixed boresight beam (no analog BF); interference already
	// includes neighbor-cell large-scale gain (path loss + random beam direction RSRP).
	// What is MISSING is the in-cell MU-MIMO inter-stream interference. In a MU group
	// of K UEs co-scheduled, each UE effectively shares the BS spatial resources;
	// per-stream effective SNR is reduced by ~K_mu × (per-stream EPA dilution).
	// Here we model this as snr_per_stream = snr_total / mx_ue_mumimo (the total
	// stream budget), so a higher R UE is penalised by needing more layers in MU.
	const Real noise_lin = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
	const Real signal_lin = dBm2linear(links[self_idx].str_signal);
	const Real intf_lin   = dBm2linear(links[self_idx].interference);
	const Real snr_total  = signal_lin / (intf_lin + noise_lin);

	// MU-aware effective SNR per stream: dilute by total stream budget (mx_ue_mumimo).
	// Rationale: in a saturated MU group K = mx_ue_mumimo / R co-scheduled UEs,
	// each stream gets a 1/mx_ue_mumimo share of the BS spatial DOF (RZF residual model).
	const int total_streams = (mx_ue_mumimo > 0) ? mx_ue_mumimo : 1;
	const Real snr_per_stream_norm = snr_total / (Real)total_streams;

	// Practical per-layer thresholds (NR CQI table style):
	//   MIN_LAYER_SNR_LIN: below this the layer cannot support even low MCS
	//     → contributes ZERO bits (would decode-fail).
	//   MIN_MARGINAL_BITS: extra layer must contribute at least this many bits/s/Hz.
	//     If marginal layer contributes less, RI stops growing — this captures the
	//     "diminishing returns for large σ-spread" intuition that pure Shannon
	//     log(1+x) sub-additivity misses (where any positive SNR keeps growing R).
	//   MAX_BITS_PER_LAYER ≈ 5.55 b/s/Hz: 256-QAM 5/6 max, MCS table cap.
	const Real MIN_LAYER_SNR_LIN  = 1.0;    //  0 dB — layer must clear noise floor handily
	const Real MIN_MARGINAL_BITS  = 1.0;    //  1 b/s/Hz minimum useful contribution
	const Real MAX_BITS_PER_LAYER = 5.55;   // 256-QAM 5/6 cap

	// Sample RBs sparsely (every stride-th RB) to amortize SVD cost.
	const int stride = (num_rb >= 10) ? (num_rb / 5) : 1;  // ~5 samples
	std::vector<Real> cap_per_rank(R_max + 1, 0.0);
	int sample_count = 0;

	for (int rb = 0; rb < num_rb; rb += stride)
	{
		JacobiSVD<MatrixXcReal> svd(H_m[0][rb], Eigen::ComputeThinU | Eigen::ComputeThinV);
		VectorXReal sv = svd.singularValues().array().abs();   // descending

		// Build per-RB "useful capacity vs R" by adding layers greedily and stopping
		// when the next layer's marginal contribution drops below MIN_MARGINAL_BITS.
		// Once a layer fails, all higher R values inherit the last useful R's value.
		Real cumulative_cap = 0.0;
		Real prev_layer_bits = MAX_BITS_PER_LAYER;  // start optimistic (R=1 always allowed)
		bool stopped = false;
		for (int R = 1; R <= R_max; R++) {
			int l = R - 1;
			Real layer_bits = 0.0;
			if (!stopped && l < (int)sv.size()) {
				Real snr_l = sv(l) * sv(l) * snr_per_stream_norm;
				// Diagnostic: log predicted per-layer SE (regardless of threshold)
				if (l < 4) {
					Real bits_raw = std::log2(1.0 + snr_l);
					if (bits_raw > MAX_BITS_PER_LAYER) bits_raw = MAX_BITS_PER_LAYER;
					pred_perlayer_se_sum[l+1] += bits_raw;
				}
				if (snr_l >= MIN_LAYER_SNR_LIN) {
					layer_bits = std::log2(1.0 + snr_l);
					if (layer_bits > MAX_BITS_PER_LAYER) layer_bits = MAX_BITS_PER_LAYER;
				}
				// Stop growing R if marginal layer contribution is below threshold
				// (except R=1 where we always allow at least one layer).
				if (l > 0 && layer_bits < MIN_MARGINAL_BITS) {
					stopped = true;
					layer_bits = 0.0;
				}
			}
			cumulative_cap += layer_bits;
			cap_per_rank[R] += cumulative_cap;
		}
		sample_count++;
	}

	if (sample_count == 0) {
		self_RI = 1;
		su_capacity = 0.0;
		return;
	}

	int best_R = 1;
	Real best_cap = cap_per_rank[1];
	for (int R = 2; R <= R_max; R++) {
		if (cap_per_rank[R] > best_cap) {
			best_cap = cap_per_rank[R];
			best_R = R;
		}
	}

	self_RI = best_R;
	su_capacity = best_cap / (Real)sample_count;  // average per RB
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

/*===================================================================
FUNCTION: MS::Quantization_of_Ch_Type2()

DESCRIPTION:
  Type II CSI codebook (TS 38.214 §5.2.2.2.3, Rel-15, rank 1).
  Reconstructs per-subband precoder as a linear combination of
  2L orthogonal DFT beams × 2 polarizations with quantized
  wideband amplitude (3-bit), subband amplitude (1-bit, optional)
  and subband phase (QPSK or 8-PSK).

  Pipeline (per coeff_idx, the outer-most rank index):
    1. Per-subband k: compute dominant right-singular vector v_k of
       H_sb[k] (sum of H_m over RBs in the subband, following the
       convention used by Quantization_of_Ch for Type I).
    2. Wideband L-beam + rotation selection:
         - For each rotation (q1,q2) ∈ O1×O2, collect N1*N2 orthogonal
           beams b_j = v_(q1+O1 n1, q2+O2 n2) / ||b_j||.
         - Score beam j = Σ_k Σ_pol |b_jᴴ · v_k^pol|².
         - Top-L sum → (q1*, q2*) with largest sum, take those L beams.
    3. SCI = argmax over 2L (L beams × 2 polarizations) of
         Σ_k |b_iᴴ · v_k^pol|².
    4. Wideband amplitude per slot (3-bit, 8 levels incl. 0):
         r_i = sqrt(Σ_k |c_i(k)|² / Σ_k |c_SCI(k)|²), quantized.
    5. Per subband k, per slot i:
         - subband amplitude q^SB_i(k) ∈ {√0.5, 1} (if SA enabled).
         - subband phase φ_i(k) quantized to nearest 2π/M grid.
    6. w_k = Σ_i b_full_i · p^WB_i · p^SB_i(k) · e^{jφ_i(k)}, then
       unit-normalized (consumer code expects ||PMI_vector||=1).

  Output: PMI_vector[coeff][rb][t] (shape NUM_TX_Port×1).
  PMI_FEEDBACK i_11/i_12/i_2 is retained as a debug stub (q1,q2,SCI).

  NOTE: Assumes BS_Mg=BS_Ng=1 (single-panel), matching the Type I
  codebook's memory layout (NUM_TX_Port = BS_P * N1 * N2, port 0..N1N2-1
  is polarization 0, port N1N2..2N1N2-1 is polarization 1).
===================================================================*/
void MS::Quantization_of_Ch_Type2(void)
{
	if (CH_CAL != 1) return;
	if (type2_beam_v == NULL) {
		// Fall back to Type I if beam dictionary was not built
		Quantization_of_Ch();
		return;
	}

	const int N1 = type2_N1;
	const int N2 = type2_N2;
	const int O1 = type2_O1;
	const int O2 = type2_O2;
	const int NN = N1 * N2;
	const int L  = g_type2_L;
	const int M_phase = g_type2_phase_alphabet;
	const bool SA = (g_type2_subband_amplitude != 0);

	// Sanity: this routine expects NUM_TX_Port = BS_P * N1 * N2 (single panel).
	// BS_P is accessed through NUM_TX_Port / NN.
	const int TxPorts = NUM_TX_Port;
	if (TxPorts % NN != 0) {
		// Dimension mismatch — fall back to Type I so simulation still runs.
		Quantization_of_Ch();
		return;
	}
	const int P = TxPorts / NN;   // 1 or 2 polarizations
	const int slotsTotal = L * P; // "2L" when P==2

	// Wideband amplitude grid (3-bit) per TS 38.214
	static const Real WB_AMP_LEVELS[8] = {
		0.0,
		std::sqrt(1.0/64.0),
		std::sqrt(1.0/32.0),
		std::sqrt(1.0/16.0),
		std::sqrt(1.0/8.0),
		std::sqrt(1.0/4.0),
		std::sqrt(1.0/2.0),
		1.0
	};
	// Subband amplitude grid (1-bit) when SA==true
	const Real SB_AMP_LEVELS[2] = { std::sqrt(0.5), 1.0 };

	// Subband configuration — match Quantization_of_Ch and Scheduling
	const int subband_size = 16;
	int num_subband = num_rb / subband_size;
	if (num_rb % subband_size != 0) num_subband += 1;

	const ComplexReal j_img(0.0, 1.0);

	// Pre-build normalized spatial beams for every (l,m):
	//   b_norm[l][m] = type2_beam_v[l][m] / sqrt(N1*N2)
	// This makes |b_norm| = 1, so projections onto v_k are bounded.
	// Done lazily per call to avoid cross-sector race (cheap: O(O1*N1*O2*N2*N1*N2)).
	// Instead of allocating, just divide each projection by sqrt(NN) at use time.

	for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
	{
		// --- Stage 1: per-subband dominant right-singular vector v_k ---
		std::vector<VectorXcReal> v_sb(num_subband);
		for (int sb = 0; sb < num_subband; sb++) {
			int rb_start = sb * subband_size;
			int rb_end   = rb_start + subband_size;
			if (rb_end > num_rb) rb_end = num_rb;

			// Equivalent wideband-over-subband: sum of S * V^H across RBs
			MatrixXcReal H_sb_eq = MatrixXcReal::Zero(NUM_RX_Port, TxPorts);
			for (int rb = rb_start; rb < rb_end; rb++) {
				JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rb], ComputeThinU | ComputeThinV);
				VectorXcReal sv = svd.singularValues();
				MatrixXcReal S  = sv.asDiagonal();
				MatrixXcReal V  = svd.matrixV();
				H_sb_eq += S * V.adjoint();
			}

			// Dominant right-singular vector of H_sb_eq ≈ wideband-eq channel projection
			JacobiSVD<MatrixXcReal> svd_sb(H_sb_eq, ComputeThinU | ComputeThinV);
			VectorXcReal v = svd_sb.matrixV().col(0);   // length TxPorts
			Real nv = v.norm();
			if (nv > 1e-12) v = v / nv;
			else            v = VectorXcReal::Zero(TxPorts);
			v_sb[sb] = v;
		}

		// --- Stage 2: rotation + L-beam selection (wideband) ---
		//
		// For each rotation (q1, q2), the N1*N2 orthogonal beams are
		//   l_n1 = q1 + O1 * n1,  m_n2 = q2 + O2 * n2
		// Score beam (n1,n2) = Σ_sb Σ_pol |b(l,m)ᴴ · v_sb^pol|² / NN
		// Pick rotation with largest sum of top-L beam scores.

		Real    best_rot_score = -1.0;
		int     best_q1 = 0, best_q2 = 0;
		std::vector<int> best_beam_n1(L, 0), best_beam_n2(L, 0);

		for (int q1 = 0; q1 < O1; q1++) {
			for (int q2 = 0; q2 < O2; q2++) {
				// Collect scores for N1*N2 orthogonal beams in this group
				std::vector<Real> scores(NN, 0.0);
				for (int n1 = 0; n1 < N1; n1++) {
					int l = q1 + O1 * n1;
					for (int n2 = 0; n2 < N2; n2++) {
						int m = q2 + O2 * n2;
						VectorXcReal& b = type2_beam_v[l][m];  // length NN
						// Accumulate energy across subbands & polarizations
						Real energy = 0.0;
						for (int sb = 0; sb < num_subband; sb++) {
							for (int p = 0; p < P; p++) {
								ComplexReal proj(0,0);
								// b^H · v_pol where v_pol = v_sb[sb].segment(p*NN, NN)
								for (int t_idx = 0; t_idx < NN; t_idx++) {
									proj += std::conj(b(t_idx)) * v_sb[sb](p * NN + t_idx);
								}
								// b is un-normalized DFT; divide by sqrt(NN) so |b_norm|=1
								energy += std::norm(proj) / (Real)NN;
							}
						}
						scores[n1 * N2 + n2] = energy;
					}
				}
				// Top-L sum
				std::vector<int> idxs(NN);
				for (int k = 0; k < NN; k++) idxs[k] = k;
				int take = (L <= NN) ? L : NN;
				std::partial_sort(idxs.begin(), idxs.begin() + take, idxs.end(),
				                  [&](int a, int b){ return scores[a] > scores[b]; });
				Real topL_sum = 0.0;
				for (int k = 0; k < take; k++) topL_sum += scores[idxs[k]];
				if (topL_sum > best_rot_score) {
					best_rot_score = topL_sum;
					best_q1 = q1;
					best_q2 = q2;
					for (int k = 0; k < take; k++) {
						best_beam_n1[k] = idxs[k] / N2;
						best_beam_n2[k] = idxs[k] % N2;
					}
					for (int k = take; k < L; k++) {  // pad if NN < L
						best_beam_n1[k] = 0;
						best_beam_n2[k] = 0;
					}
				}
			}
		}

		// Selected L beam (l,m) pairs:
		std::vector<int> beam_l(L), beam_m(L);
		for (int k = 0; k < L; k++) {
			beam_l[k] = best_q1 + O1 * best_beam_n1[k];
			beam_m[k] = best_q2 + O2 * best_beam_n2[k];
		}

		// --- Stage 3: Per-slot projections c_i(sb) for i=0..slotsTotal-1 ---
		// Slot index convention: i = pol * L + beam_idx.
		// c_i(sb) = b_spatialᴴ · v_sb^pol / sqrt(NN)
		std::vector<std::vector<ComplexReal>> cproj(slotsTotal,
		                                            std::vector<ComplexReal>(num_subband));
		std::vector<Real> wb_energy(slotsTotal, 0.0);
		for (int k = 0; k < L; k++) {
			VectorXcReal& b = type2_beam_v[beam_l[k]][beam_m[k]];
			for (int p = 0; p < P; p++) {
				int slot = p * L + k;
				for (int sb = 0; sb < num_subband; sb++) {
					ComplexReal proj(0,0);
					for (int t_idx = 0; t_idx < NN; t_idx++) {
						proj += std::conj(b(t_idx)) * v_sb[sb](p * NN + t_idx);
					}
					proj /= std::sqrt((Real)NN);
					cproj[slot][sb] = proj;
					wb_energy[slot] += std::norm(proj);
				}
			}
		}

		// SCI = argmax_slot wb_energy
		int sci_slot = 0;
		Real sci_val = wb_energy[0];
		for (int s = 1; s < slotsTotal; s++) {
			if (wb_energy[s] > sci_val) { sci_val = wb_energy[s]; sci_slot = s; }
		}
		if (sci_val <= 0.0) sci_val = 1e-30;

		// --- Stage 4: Wideband amplitude quantization ---
		std::vector<Real> p_wb(slotsTotal, 0.0);
		for (int s = 0; s < slotsTotal; s++) {
			Real raw = std::sqrt(wb_energy[s] / sci_val);   // ∈ [0, 1]
			if (s == sci_slot) { p_wb[s] = 1.0; continue; }
			// Quantize to nearest level in WB_AMP_LEVELS
			int best_idx = 0;
			Real best_d = std::fabs(raw - WB_AMP_LEVELS[0]);
			for (int q = 1; q < 8; q++) {
				Real d = std::fabs(raw - WB_AMP_LEVELS[q]);
				if (d < best_d) { best_d = d; best_idx = q; }
			}
			p_wb[s] = WB_AMP_LEVELS[best_idx];
		}

		// --- Stage 5: per-subband assembly, amplitude & phase quantization ---
		for (int sb = 0; sb < num_subband; sb++) {
			int rb_start = sb * subband_size;
			int rb_end   = rb_start + subband_size;
			if (rb_end > num_rb) rb_end = num_rb;

			// SCI reference magnitude and phase for this subband
			ComplexReal c_sci = cproj[sci_slot][sb];
			Real mag_sci = std::abs(c_sci);
			Real ph_sci  = (mag_sci > 1e-30) ? std::arg(c_sci) : 0.0;
			if (mag_sci < 1e-30) mag_sci = 1e-30;

			// Build w_sb = Σ_slot b_full_slot · amp · e^{jφ}
			VectorXcReal w_sb = VectorXcReal::Zero(TxPorts);
			for (int s = 0; s < slotsTotal; s++) {
				if (p_wb[s] == 0.0 && s != sci_slot) continue;
				int p = s / L;
				int k = s % L;
				VectorXcReal& b = type2_beam_v[beam_l[k]][beam_m[k]];

				Real coef_mag;
				Real coef_ph;
				if (s == sci_slot) {
					coef_mag = 1.0;  // p_wb=1, p_sb=1, phase=0 (reference)
					coef_ph  = 0.0;
				} else {
					ComplexReal c = cproj[s][sb];
					Real mag = std::abs(c);
					// Subband amplitude: ratio of (this-subband raw amp) / (wb amp) normalized by SCI
					Real raw_sb_ratio = (mag / mag_sci) / (p_wb[s] > 1e-30 ? p_wb[s] : 1.0);
					Real p_sb = 1.0;
					if (SA) {
						int best_idx = 0;
						Real best_d = std::fabs(raw_sb_ratio - SB_AMP_LEVELS[0]);
						for (int q = 1; q < 2; q++) {
							Real d = std::fabs(raw_sb_ratio - SB_AMP_LEVELS[q]);
							if (d < best_d) { best_d = d; best_idx = q; }
						}
						p_sb = SB_AMP_LEVELS[best_idx];
					}
					coef_mag = p_wb[s] * p_sb;
					// Phase: quantize (arg(c) - arg(c_sci)) to M-grid
					Real dph = std::arg(c) - ph_sci;
					Real step = 2.0 * pi / (Real)M_phase;
					int n_phase = (int)std::llround(dph / step);
					coef_ph = step * (Real)n_phase;
				}
				ComplexReal coef = coef_mag * std::exp(j_img * coef_ph);
				// Inject into full-port vector
				for (int t_idx = 0; t_idx < NN; t_idx++) {
					w_sb(p * NN + t_idx) += b(t_idx) * coef;
				}
			}
			// Normalize (consumer expects unit-norm PMI vector)
			Real wn = w_sb.norm();
			if (wn > 1e-12) w_sb = w_sb / wn;
			else            w_sb = VectorXcReal::Zero(TxPorts);

			for (int rb = rb_start; rb < rb_end; rb++) {
				PMI[coeff_idx][rb][t % cqi_history_length] = { best_q1, best_q2, sci_slot };
				PMI_vector[coeff_idx][rb][t % cqi_history_length] = w_sb;
			}
		}
	}
}

/*===================================================================
FUNCTION: MS::Quantization_of_Ch_EType2()

DESCRIPTION:
  Rel-16 Enhanced Type II codebook (TS 38.214 §5.2.2.2.5).
  Generates rank-R precoder per subband using:

      W_ℓ[:, k] = (1/γ_ℓ) · W₁ · W̃₂,ℓ · w_f,k*

  where
    - W₁ : 2N₁N₂ × 2L  — spatial basis (L beams × 2 pols), shared across layers.
    - W_f : N_SB × M_v — FD DFT basis with rotation q₃, shared across layers.
    - W̃₂,ℓ : 2L × M_v — sparse complex coefficients, at most
             K₀ = ⌈β · 2L · M_v⌉ nonzero entries per layer.
    - ℓ ∈ {0, ..., R-1} = layer index.

  Algorithm:

    1. Per-subband dominant R right-singular vectors V_sb[sb] (N_tx × R).

    2. Spatial basis W₁:
         - For each rotation (q₁,q₂), collect N₁N₂ orthogonal DFT beams.
         - Score beam j by Σ_ℓ Σ_sb Σ_pol |b_jᴴ · v_sb^pol,ℓ|².
         - Pick (q₁*, q₂*) and top L beams.

    3. FD basis W_f (Rel-16 novelty):
         - Compute per-(slot, subband, layer) coefficient
             c(i,sb,ℓ) = b_iᴴ · v_sb^pol,ℓ  (i = 0..2L-1).
         - Candidate FD DFT vectors: g_q3,f[sb] = e^{j 2π f (sb+q3/O3)/N3}.
         - O₃ = 4 rotations, pick (q₃*, {f*_0,...,f*_{M_v-1}}) that maximize
             Σ_ℓ Σ_i |C(i,ℓ) · g*|²  summed across selected Mv.

    4. Sparse coefficient quantization (per layer):
         - Project C_ℓ (2L × N_SB) onto W_f: C̃_ℓ = C_ℓ · W_f (2L × M_v).
         - Retain K_NZ ≤ K₀ largest-magnitude entries.
         - SCI per layer (argmax), normalized to 1.
         - Non-SCI entries: 3-bit amplitude + M-PSK phase (M = 8 or 16).

    5. Reconstruct W_ℓ[:, sb] = W₁ · C̃_ℓ · w_f,sb / unit-norm.

  Output: PMI_vector[coeff][rb][t] — MatrixXcReal(N_tx × R).

  Assumptions (same as Rel-15 Type II):
    - BS_Mg = BS_Ng = 1 (single panel).
    - N_tx = BS_P · N₁ · N₂.
===================================================================*/
void MS::Quantization_of_Ch_EType2(void)
{
	if (CH_CAL != 1) return;
	if (type2_beam_v == NULL) {
		Quantization_of_Ch();
		return;
	}

	const int N1 = type2_N1;
	const int N2 = type2_N2;
	const int O1 = type2_O1;
	const int O2 = type2_O2;
	const int NN = N1 * N2;
	const int L  = g_type2_L;
	// Per-UE rank (self_RI). When rank-adaptive mode is off this equals g_type2_rank.
	// Use R_max for memory layout (PMI_vector matrix width), R for actual fill.
	const int R_max = g_type2_rank;
	const int R     = (self_RI > 0 && self_RI <= R_max) ? self_RI : R_max;
	const int M_phase = g_type2_phase_alphabet;

	const int TxPorts = NUM_TX_Port;
	if (TxPorts % NN != 0) {
		Quantization_of_Ch();  // fallback
		return;
	}
	const int P = TxPorts / NN;   // 1 or 2 polarizations
	const int slotsTotal = L * P; // 2L when P==2

	// Subband grid
	const int subband_size = 16;
	int N_SB = num_rb / subband_size;
	if (num_rb % subband_size != 0) N_SB += 1;

	// eType II: M_v = ⌈p_v · N_SB⌉ (rank 1-2), halve for rank 3-4
	int M_v_full = (int)std::ceil(g_etype2_pv * (Real)N_SB);
	if (M_v_full < 1) M_v_full = 1;
	if (M_v_full > N_SB) M_v_full = N_SB;
	int M_v = (R <= 2) ? M_v_full : ((M_v_full + 1) / 2);
	if (M_v < 1) M_v = 1;
	if (M_v > N_SB) M_v = N_SB;

	// Sparsity budget per layer: K₀ = ⌈β · 2L · M_v⌉
	int K0 = (int)std::ceil(g_etype2_beta * (Real)slotsTotal * (Real)M_v);
	if (K0 < 1) K0 = 1;
	if (K0 > slotsTotal * M_v) K0 = slotsTotal * M_v;

	// Wideband amplitude grid (3-bit) and phase step
	static const Real WB_AMP_LEVELS[8] = {
		0.0, std::sqrt(1.0/64.0), std::sqrt(1.0/32.0), std::sqrt(1.0/16.0),
		std::sqrt(1.0/8.0), std::sqrt(1.0/4.0), std::sqrt(1.0/2.0), 1.0
	};

	const ComplexReal j_img(0.0, 1.0);

	for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
	{
		// --- Stage 1: per-subband top-R right singular vectors ---
		// V_sb[sb] has shape (TxPorts × R_eff), where R_eff = min(R, rank(H_sb)).
		std::vector<MatrixXcReal> V_sb(N_SB);
		int R_eff = R;
		for (int sb = 0; sb < N_SB; sb++) {
			int rb_start = sb * subband_size;
			int rb_end   = rb_start + subband_size;
			if (rb_end > num_rb) rb_end = num_rb;

			MatrixXcReal H_sb_eq = MatrixXcReal::Zero(NUM_RX_Port, TxPorts);
			for (int rb = rb_start; rb < rb_end; rb++) {
				JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rb], ComputeThinU | ComputeThinV);
				VectorXcReal sv = svd.singularValues();
				H_sb_eq += sv.asDiagonal() * svd.matrixV().adjoint();
			}

			JacobiSVD<MatrixXcReal> svd_sb(H_sb_eq, ComputeThinU | ComputeThinV);
			MatrixXcReal Vfull = svd_sb.matrixV();   // TxPorts × min(NUM_RX_Port, TxPorts)
			int cols_avail = (int)Vfull.cols();
			int Rs = (R < cols_avail) ? R : cols_avail;
			if (Rs < R_eff) R_eff = Rs;
			MatrixXcReal V = MatrixXcReal::Zero(TxPorts, R);
			for (int r = 0; r < Rs; r++) {
				VectorXcReal v = Vfull.col(r);
				Real nv = v.norm();
				if (nv > 1e-12) v /= nv;
				V.col(r) = v;
			}
			V_sb[sb] = V;
		}

		// --- Stage 2: rotation + L-beam selection (shared across layers) ---
		Real    best_rot_score = -1.0;
		int     best_q1 = 0, best_q2 = 0;
		std::vector<int> best_beam_n1(L, 0), best_beam_n2(L, 0);

		for (int q1 = 0; q1 < O1; q1++) {
			for (int q2 = 0; q2 < O2; q2++) {
				std::vector<Real> scores(NN, 0.0);
				for (int n1 = 0; n1 < N1; n1++) {
					int l = q1 + O1 * n1;
					for (int n2 = 0; n2 < N2; n2++) {
						int m = q2 + O2 * n2;
						VectorXcReal& b = type2_beam_v[l][m];
						Real energy = 0.0;
						for (int sb = 0; sb < N_SB; sb++) {
							for (int r = 0; r < R_eff; r++) {
								for (int p = 0; p < P; p++) {
									ComplexReal proj(0, 0);
									for (int t_idx = 0; t_idx < NN; t_idx++) {
										proj += std::conj(b(t_idx)) * V_sb[sb](p * NN + t_idx, r);
									}
									energy += std::norm(proj) / (Real)NN;
								}
							}
						}
						scores[n1 * N2 + n2] = energy;
					}
				}
				std::vector<int> idxs(NN);
				for (int k = 0; k < NN; k++) idxs[k] = k;
				int take = (L <= NN) ? L : NN;
				std::partial_sort(idxs.begin(), idxs.begin() + take, idxs.end(),
				                  [&](int a, int b){ return scores[a] > scores[b]; });
				Real topL_sum = 0.0;
				for (int k = 0; k < take; k++) topL_sum += scores[idxs[k]];
				if (topL_sum > best_rot_score) {
					best_rot_score = topL_sum;
					best_q1 = q1; best_q2 = q2;
					for (int k = 0; k < take; k++) {
						best_beam_n1[k] = idxs[k] / N2;
						best_beam_n2[k] = idxs[k] % N2;
					}
					for (int k = take; k < L; k++) {
						best_beam_n1[k] = 0; best_beam_n2[k] = 0;
					}
				}
			}
		}

		std::vector<int> beam_l(L), beam_m(L);
		for (int k = 0; k < L; k++) {
			beam_l[k] = best_q1 + O1 * best_beam_n1[k];
			beam_m[k] = best_q2 + O2 * best_beam_n2[k];
		}

		// --- Stage 3: Per-(slot,subband,layer) projection coefficients C ---
		// Slot i = pol * L + beam_k.  Shape: C[layer][slot][sb]
		std::vector<std::vector<std::vector<ComplexReal>>> C(
			R, std::vector<std::vector<ComplexReal>>(
				slotsTotal, std::vector<ComplexReal>(N_SB, ComplexReal(0, 0))));

		for (int k = 0; k < L; k++) {
			VectorXcReal& b = type2_beam_v[beam_l[k]][beam_m[k]];
			for (int p = 0; p < P; p++) {
				int slot = p * L + k;
				for (int sb = 0; sb < N_SB; sb++) {
					for (int r = 0; r < R_eff; r++) {
						ComplexReal proj(0, 0);
						for (int t_idx = 0; t_idx < NN; t_idx++) {
							proj += std::conj(b(t_idx)) * V_sb[sb](p * NN + t_idx, r);
						}
						proj /= std::sqrt((Real)NN);
						C[r][slot][sb] = proj;
					}
				}
			}
		}

		// --- Stage 4: FD basis W_f selection (q₃ rotation + Mv DFT vectors) ---
		// FD DFT vector candidate: g_{q3, f}[sb] = exp(-j 2π sb (f + q3/O3) / N_SB) / sqrt(N_SB)
		// Score(q3, f) = Σ_ℓ Σ_i |Σ_sb C[ℓ][i][sb] · g*[sb]|²
		const int O3 = 4;
		Real    best_fd_score = -1.0;
		int     best_q3 = 0;
		std::vector<int> best_fd_f(M_v, 0);

		std::vector<std::vector<ComplexReal>> g_cand(O3 * N_SB,
		                                             std::vector<ComplexReal>(N_SB));
		for (int q3 = 0; q3 < O3; q3++) {
			for (int f = 0; f < N_SB; f++) {
				for (int sb = 0; sb < N_SB; sb++) {
					Real theta = -2.0 * pi * (Real)sb * ((Real)f + (Real)q3 / (Real)O3) / (Real)N_SB;
					g_cand[q3 * N_SB + f][sb] = std::exp(j_img * theta) / std::sqrt((Real)N_SB);
				}
			}
		}

		for (int q3 = 0; q3 < O3; q3++) {
			// Score each FD vector f within this rotation
			std::vector<Real> fscore(N_SB, 0.0);
			for (int f = 0; f < N_SB; f++) {
				Real sum_sq = 0.0;
				std::vector<ComplexReal>& g = g_cand[q3 * N_SB + f];
				for (int r = 0; r < R_eff; r++) {
					for (int s = 0; s < slotsTotal; s++) {
						ComplexReal inner(0, 0);
						for (int sb = 0; sb < N_SB; sb++) {
							inner += C[r][s][sb] * std::conj(g[sb]);
						}
						sum_sq += std::norm(inner);
					}
				}
				fscore[f] = sum_sq;
			}
			// Top-M_v FD basis within this rotation
			std::vector<int> fidxs(N_SB);
			for (int f = 0; f < N_SB; f++) fidxs[f] = f;
			int takef = (M_v <= N_SB) ? M_v : N_SB;
			std::partial_sort(fidxs.begin(), fidxs.begin() + takef, fidxs.end(),
			                  [&](int a, int b){ return fscore[a] > fscore[b]; });
			Real rot_sum = 0.0;
			for (int t = 0; t < takef; t++) rot_sum += fscore[fidxs[t]];
			if (rot_sum > best_fd_score) {
				best_fd_score = rot_sum;
				best_q3 = q3;
				for (int t = 0; t < takef; t++) best_fd_f[t] = fidxs[t];
				for (int t = takef; t < M_v; t++) best_fd_f[t] = 0;
			}
		}

		// Assemble W_f: N_SB × M_v (each column is a selected FD DFT vector)
		MatrixXcReal W_f = MatrixXcReal::Zero(N_SB, M_v);
		for (int f = 0; f < M_v; f++) {
			int fidx = best_fd_f[f];
			std::vector<ComplexReal>& g = g_cand[best_q3 * N_SB + fidx];
			for (int sb = 0; sb < N_SB; sb++) {
				W_f(sb, f) = g[sb];
			}
		}

		// --- Stage 5: per-layer sparse coefficient quantization ---
		// C̃_ℓ = C_ℓ · W_f  (2L × M_v). Keep K₀ largest-|·| entries with SCI normalization.
		std::vector<MatrixXcReal> C_tilde(R, MatrixXcReal::Zero(slotsTotal, M_v));
		for (int r = 0; r < R_eff; r++) {
			for (int s = 0; s < slotsTotal; s++) {
				for (int f = 0; f < M_v; f++) {
					ComplexReal v(0, 0);
					for (int sb = 0; sb < N_SB; sb++) {
						v += C[r][s][sb] * W_f(sb, f);
					}
					C_tilde[r](s, f) = v;
				}
			}
		}

		std::vector<MatrixXcReal> W2_tilde(R, MatrixXcReal::Zero(slotsTotal, M_v));

		for (int r = 0; r < R_eff; r++) {
			// SCI per layer = argmax_|·|
			int sci_s = 0, sci_f = 0;
			Real sci_val = 0.0;
			for (int s = 0; s < slotsTotal; s++) {
				for (int f = 0; f < M_v; f++) {
					Real a = std::abs(C_tilde[r](s, f));
					if (a > sci_val) { sci_val = a; sci_s = s; sci_f = f; }
				}
			}
			if (sci_val < 1e-30) sci_val = 1e-30;
			Real sci_phase = std::arg(C_tilde[r](sci_s, sci_f));

			// Flatten into (value, index) and sort by magnitude to pick top K₀
			std::vector<std::pair<Real, int>> mag_idx(slotsTotal * M_v);
			for (int s = 0; s < slotsTotal; s++) {
				for (int f = 0; f < M_v; f++) {
					Real ratio = std::abs(C_tilde[r](s, f)) / sci_val;
					mag_idx[s * M_v + f] = std::make_pair(ratio, s * M_v + f);
				}
			}
			int Kkeep = (K0 <= (int)mag_idx.size()) ? K0 : (int)mag_idx.size();
			std::partial_sort(mag_idx.begin(), mag_idx.begin() + Kkeep, mag_idx.end(),
			                  [](const std::pair<Real,int>& a, const std::pair<Real,int>& b){
			                      return a.first > b.first;
			                  });
			std::vector<bool> keep(slotsTotal * M_v, false);
			for (int t_idx = 0; t_idx < Kkeep; t_idx++) keep[mag_idx[t_idx].second] = true;
			keep[sci_s * M_v + sci_f] = true;  // always keep SCI

			// Quantize each kept coefficient
			for (int s = 0; s < slotsTotal; s++) {
				for (int f = 0; f < M_v; f++) {
					if (!keep[s * M_v + f]) continue;
					if (s == sci_s && f == sci_f) {
						W2_tilde[r](s, f) = ComplexReal(1.0, 0.0);
						continue;
					}
					ComplexReal c = C_tilde[r](s, f);
					Real raw_mag = std::abs(c) / sci_val;
					// Amplitude: 3-bit quantization
					int best_i = 0;
					Real best_d = std::fabs(raw_mag - WB_AMP_LEVELS[0]);
					for (int q = 1; q < 8; q++) {
						Real d = std::fabs(raw_mag - WB_AMP_LEVELS[q]);
						if (d < best_d) { best_d = d; best_i = q; }
					}
					Real p_amp = WB_AMP_LEVELS[best_i];
					if (p_amp == 0.0) continue;   // pruned as zero
					// Phase: (M_phase)-PSK quantization, referenced to SCI phase
					Real dph = std::arg(c) - sci_phase;
					Real step = 2.0 * pi / (Real)M_phase;
					int n_ph = (int)std::llround(dph / step);
					Real q_ph = step * (Real)n_ph;
					W2_tilde[r](s, f) = p_amp * std::exp(j_img * q_ph);
				}
			}
		}

		// --- Stage 6: reconstruct W_ℓ[:, sb] per subband, per layer ---
		// w_ℓ,sb = W₁ · W̃₂,ℓ · W_f[sb, :]ᴴ  (length N_tx)
		// W₁ is implicit: slot s = (pol p, beam k) → [0_{p*NN}; b_k; 0_{(1-p)*NN}]
		for (int sb = 0; sb < N_SB; sb++) {
			int rb_start = sb * subband_size;
			int rb_end   = rb_start + subband_size;
			if (rb_end > num_rb) rb_end = num_rb;

			// w_f_sb = W_f[sb, :]^H  (column M_v), i.e. complex conjugate of row
			MatrixXcReal W_rank = MatrixXcReal::Zero(TxPorts, R);
			for (int r = 0; r < R_eff; r++) {
				// coef_per_slot[s] = Σ_f W̃₂(s, f) · conj(W_f(sb, f))
				VectorXcReal coef_s = VectorXcReal::Zero(slotsTotal);
				for (int s = 0; s < slotsTotal; s++) {
					ComplexReal v(0, 0);
					for (int f = 0; f < M_v; f++) {
						v += W2_tilde[r](s, f) * std::conj(W_f(sb, f));
					}
					coef_s(s) = v;
				}
				// Assemble per-port w_r
				VectorXcReal w = VectorXcReal::Zero(TxPorts);
				for (int s = 0; s < slotsTotal; s++) {
					if (std::abs(coef_s(s)) < 1e-30) continue;
					int p = s / L;
					int k = s % L;
					VectorXcReal& b = type2_beam_v[beam_l[k]][beam_m[k]];
					for (int t_idx = 0; t_idx < NN; t_idx++) {
						w(p * NN + t_idx) += b(t_idx) * coef_s(s);
					}
				}
				Real wn = w.norm();
				if (wn > 1e-12) w /= wn;
				W_rank.col(r) = w;
			}

			// (Optional) orthogonalize layers via thin QR for rank > 1
			if (R_eff > 1) {
				Eigen::HouseholderQR<MatrixXcReal> qr(W_rank.leftCols(R_eff));
				MatrixXcReal Q = qr.householderQ() * MatrixXcReal::Identity(TxPorts, R_eff);
				for (int r = 0; r < R_eff; r++) {
					VectorXcReal q = Q.col(r);
					Real qn = q.norm();
					if (qn > 1e-12) q /= qn;
					W_rank.col(r) = q;
				}
			}

			// --- Diagnostic: PMI quality vs true SVD (per subband) ---
			// Compare quantized W_rank columns to the true right singular vectors V_sb[sb].
			// pmi_sigma2_sum[l]  : cumulative σ_l² (rank-spread analysis)
			// pmi_overlap_sum[l] : cumulative |w_l^H v_l|^2 (codebook quantization quality)
			{
				JacobiSVD<MatrixXcReal> sv_chk(V_sb[sb], ComputeThinU | ComputeThinV);
				// V_sb stored is already TxPorts × R unit-norm (from Stage 1).
				// σ_l² is derived from the original H_sb_eq SVD; we re-recover it from
				// the H_m channel for this subband.
				MatrixXcReal H_sb_for_stats = MatrixXcReal::Zero(NUM_RX_Port, TxPorts);
				int rb_start_s = sb * subband_size;
				int rb_end_s   = rb_start_s + subband_size;
				if (rb_end_s > num_rb) rb_end_s = num_rb;
				for (int rb = rb_start_s; rb < rb_end_s; rb++) {
					JacobiSVD<MatrixXcReal> svd(H_m[coeff_idx][rb], ComputeThinU | ComputeThinV);
					H_sb_for_stats += svd.singularValues().asDiagonal() * svd.matrixV().adjoint();
				}
				JacobiSVD<MatrixXcReal> svd_stats(H_sb_for_stats, ComputeThinU | ComputeThinV);
				VectorXReal sv_stats = svd_stats.singularValues().array().abs();
				MatrixXcReal V_true  = svd_stats.matrixV();

				int R_for_stats = (R_eff <= (int)sv_stats.size()) ? R_eff : (int)sv_stats.size();
				for (int l = 0; l < R_for_stats && l < 4; l++) {
					pmi_sigma2_sum[l+1] += sv_stats(l) * sv_stats(l);
					// Per-column overlap (legacy): |w_rank.col(l)^H * v_true.col(l)|^2
					if (l < (int)W_rank.cols() && l < (int)V_true.cols()) {
						ComplexReal inner(0,0);
						for (int t_idx = 0; t_idx < TxPorts; t_idx++) {
							inner += std::conj(W_rank(t_idx, l)) * V_true(t_idx, l);
						}
						pmi_overlap_sum[l+1] += std::norm(inner);
					}
				}

				// Principal angles via SVD(W_rank^H · V_true_R): R×R matrix.
				// cos²(θ_l) = σ_l²(W^H V); R - Σ cos² = subspace chordal distance squared.
				if (R_for_stats > 0 && R_for_stats <= (int)W_rank.cols()) {
					MatrixXcReal V_red = V_true.leftCols(R_for_stats);
					MatrixXcReal W_red = W_rank.leftCols(R_for_stats);
					MatrixXcReal G = W_red.adjoint() * V_red;       // R×R
					JacobiSVD<MatrixXcReal> svd_g(G);
					VectorXReal sv_g = svd_g.singularValues().array().abs();
					Real sum_cos2 = 0.0;
					for (int l = 0; l < R_for_stats && l < 4; l++) {
						Real c2 = (l < (int)sv_g.size()) ? (Real)(sv_g(l) * sv_g(l)) : 0.0;
						if (c2 > 1.0) c2 = 1.0;  // numerical safety
						pmi_cos2_sum[l+1] += c2;
						sum_cos2 += c2;
					}
					// chordal distance² = R - sum_cos²
					Real chord_sq = (Real)R_for_stats - sum_cos2;
					if (chord_sq < 0.0) chord_sq = 0.0;
					pmi_chordal_sum += std::sqrt(chord_sq);
				}

				pmi_quality_count++;
			}

			for (int rb = rb_start; rb < rb_end; rb++) {
				PMI[coeff_idx][rb][t % cqi_history_length] = { best_q1, best_q2, best_q3 };
				PMI_vector[coeff_idx][rb][t % cqi_history_length] = W_rank;
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

					// Use dominant layer (col 0) for CQI estimation. Multi-layer CQI is deferred.
					VectorXcReal w = PMI_vector[coeff_idx][rep_rb_idx][t % cqi_history_length].col(0).conjugate(); // Tx1

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

						H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length].col(0);
						u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
						A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
						C = (linear_interference + noise) * (u * u.adjoint()).norm();

						temp_sinr = A / C;
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

						H_bar = sqrt(linear_signal) * (H_m[coeff_idx][rb_idx]) * PMI_vector[coeff_idx][rb_idx][t % cqi_history_length].col(0);

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

		H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length].col(0);
		Real comp_sector_static_gain = dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]); //dBm2linear(links[self_idx].static_gain[1].first); // dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
		H_inf = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length].col(0);

		u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + H_inf * H_inf.adjoint() + (linear_interference + noise) * Identity).inverse();
		A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
		B = (u * H_inf * H_inf.adjoint() * u.adjoint()).norm();
		C = (linear_interference + noise) * (u * u.adjoint()).norm();		

		CQI[0][rb_idx][abs(t) % cqi_history_length] = A/(B+C);

	}
}

#if 0  // unused: Fourier_Transform_of_Channel — replaced by Fourier_Transform_of_Channel_Optimized
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
#endif  // unused: Fourier_Transform_of_Channel

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

						for (int path_idx = 0; path_idx < NUM_PATH; path_idx++)
						{
							const complex<Real>& gain = channel[_bs_idx][self_idx].CHIR[sec_idx][t_idx][r_idx][path_idx];
							Hc += gain * phasor[path_idx];
							phasor[path_idx] *= step[path_idx];
						}

						H_m[coeff_idx][rbs_idx](r_idx, t_idx) = Hc;
					}

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
// H_usn Fourier Transform (no BF weights)
// Converts H_usn[cluster](totalRx, totalTx) → H_m_elem[rb](totalRx, totalTx)
// ====================================================================

#if 0  // unused: Fourier_Transform_H_usn — never called
void MS::Fourier_Transform_H_usn(int ms_idx)
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
#endif  // unused: Fourier_Transform_H_usn

// ====================================================================
// Element-Level DFT + Beamforming Weights → H_m (port-level)
// ====================================================================
// Data flow:
//   H_usn[cluster](totalRx, totalTx) → DFT per RB → H_elem(totalRx, totalTx)
//   H_m[coeff][rb] = W_rx^H × H_elem × W_tx   (NUM_RX_Port × NUM_TX_Port)
//
// W_tx: (totalTx × NUM_TX_Port) beamforming weight matrix
// W_rx: (totalRx × NUM_RX_Port) beamforming weight matrix
// ====================================================================
void MS::Fourier_Transform_WithBF(int ms_idx)
{
	if (CH_CAL != 1) return;

	int self_idx = ms_idx;
	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;
	int K = BS_M / BS_Mp;
	int L = BS_N / BS_Np;
	int ue_K = MS_M / MS_Mp;
	int ue_L = MS_N / MS_Np;

	for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
	{
		int _bs_idx;
		int adj_sector;
		if (TYPE == 11 && num_Indoor_TRxP == 1) {
			_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx]);
			adj_sector = _bs_idx;
		} else {
			_bs_idx = (int)(links[self_idx].adj_sector[coeff_idx] / 3);
			adj_sector = links[self_idx].adj_sector[coeff_idx];
		}

		CHANNEL& ch = channel[_bs_idx][self_idx];
		const int NUM_PATH = ch.NUM_PATH_for_channelcoeff;
		if (NUM_PATH <= 0) continue;
		if ((int)ch.H_usn.size() < NUM_PATH) continue;

		// --- Build W_tx: (totalTx × NUM_TX_Port) ---
		// TX port index: p*Mp*Np*Mg*Ng + mi*Np*Mg*Ng + ni*Mg*Ng + mg*Ng + ng
		// TX element flat index: s_m*(BS_N*BS_P) + s_n*BS_P + s_p
		// Within port (mi,ni,p,mg,ng): elements are (mi*K+k, ni*L+l, p, mg, ng) for k=0..K-1, l=0..L-1
		MatrixXcReal W_tx = MatrixXcReal::Zero(totalTx, NUM_TX_Port);
		{
			int beam_z = links[self_idx].analog_beam_selection[adj_sector].sector_z;
			int beam_a = links[self_idx].analog_beam_selection[adj_sector].sector_a;

			for (int p = 0; p < BS_P; p++) {
				for (int mi = 0; mi < BS_Mp; mi++) {
					// Row-beam mode: use the sector-wide per-row zenith beam instead of
					// the per-link best. adj_sector covers interfering sectors too, so
					// interference channels automatically carry the interferer's actual
					// assigned beams. Per-port normalization (1/sqrt(K*L)) is baked into
					// the weights and beam-independent — no renormalization needed.
					int bz = row_beam_enable ? sector[adj_sector].row_beam_z[mi] : beam_z;
					for (int ni = 0; ni < BS_Np; ni++) {
						int ba = beam_a;  // az_mode 0: keep per-UE azimuth
						if (row_beam_enable && row_beam_az_mode == 1) ba = row_beam_boresight_a;
						if (row_beam_enable && row_beam_az_mode == 2) ba = sector[adj_sector].col_beam_a[ni];

						// Port flat index (Mg=1, Ng=1 for current configs)
						int port_idx = p * BS_Mp * BS_Np + mi * BS_Np + ni;

						for (int k = 0; k < K; k++) {
							for (int l = 0; l < L; l++) {
								int s_m = mi * K + k;
								int s_n = ni * L + l;
								// Element flat index
								int elem_idx = s_m * (BS_N * BS_P) + s_n * BS_P + p;

								W_tx(elem_idx, port_idx) = virtualization_weight_wv[bz][ba][k][l];
							}
						}
					}
				}
			}
		}

		// --- Build W_rx: (totalRx × NUM_RX_Port) ---
		MatrixXcReal W_rx = MatrixXcReal::Zero(totalRx, NUM_RX_Port);
		{
			if (ue_antenna_element_gain == 0) {
				// Omni UE: identity mapping (1 element per port)
				int min_dim = std::min(totalRx, NUM_RX_Port);
				for (int i = 0; i < min_dim; i++) {
					W_rx(i, i) = ComplexReal(1.0, 0.0);
				}
			} else {
				int beam_z = links[self_idx].analog_beam_selection[adj_sector].z;
				int beam_a = links[self_idx].analog_beam_selection[adj_sector].a;

				for (int p = 0; p < MS_P; p++) {
					for (int mi = 0; mi < MS_Mp; mi++) {
						for (int ni = 0; ni < MS_Np; ni++) {
							int port_idx = p * MS_Mp * MS_Np + mi * MS_Np + ni;

							for (int k = 0; k < ue_K; k++) {
								for (int l = 0; l < ue_L; l++) {
									int u_m = mi * ue_K + k;
									int u_n = ni * ue_L + l;
									int elem_idx = u_m * (MS_N * MS_P) + u_n * MS_P + p;

									W_rx(elem_idx, port_idx) = ue_virtualization_weight_wv[beam_z][beam_a][k][l];
								}
							}
						}
					}
				}
			}
		}

		// --- DFT + BF per RB ---
		// FFT parameters (same as Fourier_Transform_of_Channel_Optimized)
		int pointer = (int(fft_size) - (num_rb * num_freq_per_rbs)) / 2;
		const unsigned long long f0 = (pointer + (num_freq_per_rbs / 2)) * subcarrier_spacing;
		const Real df   = subcarrier_spacing;
		const Real dfrb = num_freq_per_rbs * df;

		// Initialize phasor and step per path
		std::vector<complex<Real>> phasor(NUM_PATH);
		std::vector<complex<Real>> step_val(NUM_PATH);

		for (int path_idx = 0; path_idx < NUM_PATH; path_idx++) {
			const Real tau = ch.delay[path_idx];
			phasor[path_idx] = exp(complex<Real>(0.0, -2.0 * pi * f0 * tau));
			step_val[path_idx] = exp(complex<Real>(0.0, -2.0 * pi * dfrb * tau));
		}

		// Pre-compute W_rx^H (adjoint)
		MatrixXcReal W_rx_H = W_rx.adjoint();  // NUM_RX_Port × totalRx

		for (int rbs_idx = 0; rbs_idx < num_rb; rbs_idx++)
		{
			// DFT: H_elem = Σ_path H_usn[path] * phasor[path]
			MatrixXcReal H_elem = MatrixXcReal::Zero(totalRx, totalTx);
			for (int path_idx = 0; path_idx < NUM_PATH; path_idx++) {
				H_elem += ch.H_usn[path_idx] * phasor[path_idx];
				phasor[path_idx] *= step_val[path_idx];
			}

			// BF: H_m[coeff][rb] = W_rx^H × H_elem × W_tx
			H_m[coeff_idx][rbs_idx] = W_rx_H * H_elem * W_tx;
		}
	}
}

