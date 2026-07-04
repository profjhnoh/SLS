#include "common.h"
#include "const.h"
#include "nr_l2sm.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void Averaging_Sinr(int ue_idx, Real sinr);
void Accumulate_Precoding_Metrics(int ue_idx, Real wideband_sinr_linear, Real coupling_loss_dB);
Real Get_distance(LOCATION, LOCATION);

Real test1 = 12;
void MS::Receive_DL_mTRP( void )
{
	Find_Allocated_Rbs_Mcs_mTRP();
	if (rb_indices.size() > 0)
	{
		if (0) //( self_idx == 1)
		{
			cout << "self_idx = " << self_idx << endl;
		}

		Real bler = 1;
		Compute_RBs_SINR();
		Compute_Transport_Block_Size();
		Compute_Effective_SINR();
		bler = BLER_Mapping();

		// Collect precoding metrics (wideband SINR and coupling loss)
		// 3GPP "Wideband SINR before receiver" definition:
		// - TX precoding (w) applied: YES
		// - RX beamforming (u): NO
		// - MU-MIMO intra-cell interference: NO (only inter-cell)
		// Formula: SINR = ||H * w||^2 * P / (I_intercell + N)
		if (rb_indices.size() > 0 && t > N_pf)
		{
			int serving_sector = links[self_idx]._sector_in_control;
			Real noise_linear = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
			Real linear_signal = dBm2linear(bs_maxpower - channel[(int)(serving_sector/3)][self_idx].pathloss_final);
			Real linear_interference = dBm2linear(links[self_idx].interference);  // Inter-cell only (no MU-MIMO)

			// Compute 3GPP Wideband SINR: average ||H * w||^2 across scheduled RBs
			Real sum_eff_channel_gain = 0.0;
			int valid_rb_count = 0;

			for (int freq_idx = 0; freq_idx < (int)rb_indices.size(); freq_idx++)
			{
				int rb_idx = rb_indices[freq_idx];

				// Find stream number for this UE
				int stream_num = -1;
				for (int i = 0; i < mx_ue_mumimo; i++)
				{
					if (ppSchedulerRead[rb_idx][i].ue_selected == self_idx)
					{
						stream_num = i;
						break;
					}
				}
				if (stream_num == -1) continue;

				// Get channel and precoding matrices
				const MatrixXcReal& H = H_m[0][rb_idx];
				const MatrixXcReal& W = sector[serving_sector].W[rb_idx];

				// Effective channel vector: h_eff = H * w (TX precoding only, no RX beamforming)
				const MatrixXcReal h_eff = H * W.col(stream_num);  // [NUM_RX_Port x 1]

				// Channel power: ||H * w||^2 (Frobenius norm squared)
				Real eff_channel_gain = h_eff.squaredNorm();
				sum_eff_channel_gain += eff_channel_gain;
				valid_rb_count++;
			}

			if (valid_rb_count > 0)
			{
				// Average effective channel gain across RBs
				Real avg_eff_channel_gain = sum_eff_channel_gain / valid_rb_count;

				// 3GPP Wideband SINR (before receiver):
				// SINR = P * ||H * w||^2 / (I_intercell + N)
				Real wideband_sinr_linear = linear_signal * avg_eff_channel_gain / (linear_interference + noise_linear);

				// Coupling Loss: pathloss with TX precoding gain
				// Effective gain (dB) = 10*log10(||H * w||^2)
				Real pathloss = channel[(int)(serving_sector/3)][self_idx].pathloss_final;
				Real eff_gain_dB = 10.0 * log10(avg_eff_channel_gain);
				Real coupling_loss_dB = -pathloss + eff_gain_dB;

				Accumulate_Precoding_Metrics(self_idx, wideband_sinr_linear, coupling_loss_dB);
			}
		}
		if (t > N_pf)
			AddThroughput(bler);
	}
}

// Helper: lookup BLER for a given linear SINR and CQI index (1..15)
static inline Real bler_lookup_dB(Real sinr_linear, int cqi_idx)
{
	Real x_axis = 10 * std::log10(std::max(sinr_linear, (Real)1e-30));
	int mapping_idx = 1425;
	for (int snr_idx = 0; snr_idx < 1426; snr_idx++)
	{
		if (SNR_5G_dB[snr_idx] > x_axis) {
			mapping_idx = snr_idx - 1;
			if (mapping_idx < 0) mapping_idx = 0;
			break;
		}
	}
	int c = cqi_idx - 1;
	if (c < 0) c = 0;
	if (c > 14) c = 14;
	return NEW5GBLER[mapping_idx][c];
}

// Helper: NR transport-block-size (info bits) from MCS index, resource-element count N_RE,
// and a layer multiplier. Mirrors the legacy Compute_Transport_Block_Size math (TS 38.214 5.1.3.2).
// Used by the per-layer TBS path (layer_factor = 1 per layer).
static int tbs_info_bits_from_mcs(int mcs_idx, int N_RE, int layer_factor)
{
	// Math lives in nr_l2sm.cpp (single source shared with the CQI-threshold
	// regeneration). Bit-identity vs. old binaries is NOT a constraint here:
	// the codebase has a layout-sensitive read (adding ANY translation unit
	// perturbs a few UEs' trajectories — verified with a dummy-TU experiment),
	// so regression is judged by same-binary determinism + statistics.
	return nr_tbs_info_bits_from_mcs(mcs_idx, N_RE, layer_factor);
}

Real MS::BLER_Mapping( void )
{
	// Per-layer BLER lookup.
	// TB BLER = 1 - Π_{l=1..R} (1 - BLER_l)   (independent per-layer decode model)
	// Legacy (g_per_layer_mcs==0) keeps the original >1-layer trigger; the single-layer
	// case falls back to the ESINR_linear lookup below.
	// Per-layer mode (==1) MUST populate per_layer_BLER even for a single layer, otherwise
	// AddThroughput reads a stale per_layer_BLER[0]==0 and "always decodes" (throughput blow-up).
	int min_layers_for_perlayer = (g_per_layer_mcs == 1) ? 1 : 2;
	if (num_layers_actual >= min_layers_for_perlayer && (int)rbs_rx_per_layer.size() == num_layers_actual)
	{
		// Per-layer EESM aggregation across RBs
		Real prob_success_TB = 1.0;
		Real fb_err_sum = 0.0;   // realized-ESINR feedback: innovation accumulator (dB)
		int  fb_err_cnt = 0;
		for (int l = 0; l < num_layers_actual && l < 4; l++)
		{
			const std::vector<Real>& v = rbs_rx_per_layer[l];
			if (v.empty()) { per_layer_BLER[l] = 1.0; prob_success_TB = 0.0; continue; }

			// Per-layer MCS: score each layer against ITS OWN cqi, not the dominant layer's.
			int cqi_l = (g_per_layer_mcs == 1) ? _cqi_per_layer[l] : _cqi_idx;
			int mcs_l = (g_per_layer_mcs == 1) ? _mcs_per_layer[l] : _mcs_idx;

			// Effective SINR over this layer's RBs.
			// MATLAB mode uses the RBIR ESM — the AWGN curves are tabulated against an
			// RBIR-effective SINR, and feeding them EESM(legacy Beta) biases the steep
			// LDPC waterfalls (verified: BLER stuck ~0.2, SE collapse). Legacy keeps EESM.
			Real esir_l;
			if (g_matlab_bler && g_matlab_rbir)
			{
				int  qm_l;
				Real rr_l;
				nr_mcs_to_qm_r(mcs_l, qm_l, rr_l);
				esir_l = rbir_effective_sinr_linear(v, qm_l);
			}
			else
			{
				Real beta = Beta[std::max(cqi_l - 1, 0)];
				Real e_sir_l = 0.0;
				for (size_t i = 0; i < v.size(); i++) e_sir_l += std::exp(-v[i] / beta);
				esir_l = -beta * std::log(e_sir_l / (Real)v.size());
			}

			Real bler_l;
			if (g_matlab_bler)
			{
				// (code rate, TBS)-aware lookup. Per-layer mode has its own MCS/TBS;
				// legacy multi-layer shares _mcs_idx and splits _info_bits evenly
				// (approximation: TBS only sets BGN/Zc/C — low sensitivity).
				int tbs_l = (g_per_layer_mcs == 1) ? _info_bits_per_layer[l]
				                                   : (int)(_info_bits / MAX(num_layers_actual, 1));
				bler_l = bler_lookup_matlab(esir_l, mcs_l, tbs_l);
			}
			else
			{
				bler_l = bler_lookup_dB(esir_l, cqi_l);
			}
			per_layer_ESINR[l] = esir_l;
			per_layer_BLER [l] = bler_l;
			prob_success_TB *= (1.0 - bler_l);

			// Realized-ESINR feedback innovation: realized (RBIR) vs scheduler estimate.
			// Uses the raw pre-OLLA estimate so the correction's fixpoint is exactly the
			// estimate-vs-realized bias (OLLA stays the residual BLER-level regulator).
			if (g_matlab_bler && g_matlab_esinr_fb && g_per_layer_mcs == 1 &&
			    _est_sinr_per_layer[l] > 0 && esir_l > 0)
			{
				fb_err_sum += 10.0 * std::log10(esir_l) - 10.0 * std::log10(_est_sinr_per_layer[l]);
				fb_err_cnt++;
			}
		}

		// Integral update of the per-UE correction: corr -> EWMA of (realized - estimate) dB.
		// Gain 0.1/slot converges in tens of scheduled slots; clamp keeps outliers bounded.
		if (fb_err_cnt > 0)
		{
			Real err = fb_err_sum / (Real)fb_err_cnt;
			matlab_sinr_corr += (Real)0.1 * (err - matlab_sinr_corr);
			if (matlab_sinr_corr < (Real)-15.0) matlab_sinr_corr = (Real)-15.0;
			if (matlab_sinr_corr > (Real)5.0)   matlab_sinr_corr = (Real)5.0;
		}

		Real BLER_value = 1.0 - prob_success_TB;
		ue_BLER_Value = BLER_value;
		return BLER_value;
	}

	// Legacy single-SINR fallback (R=1 or per-layer buffer not populated)
	Real esinr_fb = ESINR_linear;   // EESM value from Compute_Effective_SINR
	if (g_matlab_bler && g_matlab_rbir && rbs_rx.size() > 0)
	{
		int  qm_fb;
		Real rr_fb;
		nr_mcs_to_qm_r(_mcs_idx, qm_fb, rr_fb);
		esinr_fb = rbir_effective_sinr_linear(rbs_rx, qm_fb);
	}
	Real BLER_value = g_matlab_bler ? bler_lookup_matlab(esinr_fb, _mcs_idx, (int)_info_bits)
	                                : bler_lookup_dB(ESINR_linear, _cqi_idx);
	ue_BLER_Value = BLER_value;
	return BLER_value;
}

/*===================================================================
FUNCTION: UE::Receive_DL( void )

AUTHOR: SJ, KR, DW

DESCRIPTION: UE check the received resource block and check if it has failed or not.

NOTES:


===================================================================*/
void MS::Receive_DL(int _ms_idx)
{

	receive_self_ms_idx = _ms_idx;

	// initialize capacity
	Find_Allocated_Rbs_Mcs(receive_self_ms_idx);
	if (rb_indices.size() > 0)
	{
		Rate_matching();
		Compute_Num_Bin_Symbol();
		Compute_Effective_SINR();
		////////// BLER mapping //////////////////
		// (the legacy scan below was inlined; bler_lookup_dB is scan-identical,
		//  so matlab_bler=0 stays bit-identical)
		int mapping_cqi_idx = 0;
		Real BLER_value = g_matlab_bler ? bler_lookup_matlab(ESINR_linear, _mcs_idx, (int)_info_bits)
		                                : bler_lookup_dB(ESINR_linear, _cqi_idx);
		ue_BLER_Value  = BLER_value;



		if (t > N_pf)
		{
			AddThroughput(BLER_value);
		}			
	}
}

/*===================================================================
FUNCTION: UE::AddThroughput( void )

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void MS::AddThroughput(Real prob)
{
#ifdef ENABLE_MULTITHREADING
	// Per-thread RNG: AddThroughput runs inside the OpenMP receive region
	// (measure.cpp), where the global `randnum` would be raced by all threads
	// (non-atomic read-modify-write of its p,q,r state = C++ undefined behaviour).
	// Shadowing it with a thread-local instance removes the race. NOTE: this is a
	// CORRECTNESS fix only — it does NOT improve run-to-run reproducibility, because
	// the thread-local RNG is seeded by thread-id and the OpenMP work distribution
	// varies per run (A/B tested: ~17% 1-drop variance unchanged). True reproducibility
	// would require seeding the HARQ coin deterministically per work-item (ue/slot).
	Rand& randnum = get_thread_local_rng();
#endif
	//for check
	if (prob < 0 || prob > 1)
	{
		cout << "something wrong with prob_failure at ue.cpp. prob = " << prob; getchar();
	}

	Real tossing_coin = randnum.u();	

	int my_sec_idx = links[self_idx].adj_sector[0];
	links[self_idx].total_num_scheduled++;
	if (links[self_idx].num_re_tx == 0)
	{
		links[self_idx].total_num_tx++;
		sum_cqi             += _cqi_idx;
		sum_mcs_type        += _mcs_idx;
		total_esinr         += ESINR_linear;
		total_estimate_sinr += _avr_sinr;
	}

	// Helper lambda for OLLA moving window update
	auto update_olla_window = [this](bool is_ack) {
		if (!g_olla_enable || g_olla_window_size <= 0) return;

		// Safety check: ensure olla_history is properly sized
		if (olla_history.size() != (size_t)g_olla_window_size) {
			olla_history.clear();
			olla_history.resize(g_olla_window_size, true);
			olla_history_idx = 0;
			olla_nack_count = 0;
		}

		// Safety check: ensure index is within bounds
		if (olla_history_idx < 0 || olla_history_idx >= (int)olla_history.size()) {
			olla_history_idx = 0;
		}

		// Update NACK count: remove old value, add new value
		bool old_value = olla_history[olla_history_idx];
		if (!old_value) olla_nack_count--;  // Old was NACK, remove from count
		if (!is_ack) olla_nack_count++;     // New is NACK, add to count

		// Store new value and advance index
		olla_history[olla_history_idx] = is_ack;
		olla_history_idx = (olla_history_idx + 1) % g_olla_window_size;

		// Calculate measured BLER from moving window
		Real measured_bler = (Real)olla_nack_count / (Real)g_olla_window_size;
		Real bler_diff = measured_bler - g_olla_target_bler;

		// Adjust offset only if difference exceeds margin (hysteresis)
		if (bler_diff > g_olla_bler_margin)
		{
			// BLER too high: decrease offset (more conservative MCS)
			olla_offset -= g_olla_step_down;
			if (olla_offset < g_olla_min_offset)
				olla_offset = g_olla_min_offset;
		}
		else if (bler_diff < -g_olla_bler_margin)
		{
			// BLER too low: increase offset (more aggressive MCS)
			olla_offset += g_olla_step_up;
			if (olla_offset > g_olla_max_offset)
				olla_offset = g_olla_max_offset;
		}
	};

	// ============================================================================
	// PER-LAYER HARQ (g_per_layer_mcs==1): each layer decodes independently with its
	// own BLER. A weak layer failing does NOT NACK the others; it just retransmits
	// (up to 3) on its own. Throughput credits each decoded layer's sub-TBS. The
	// TB-level ACK is "all active layers resolved", which gates the buffer reset.
	// ============================================================================
	if (g_per_layer_mcs == 1)
	{
		int R = (num_layers_actual > 0) ? num_layers_actual : 1;
		if (R > 4) R = 4;
		bool first_tx = (links[self_idx].num_re_tx == 0);
		bool layer_ok_first[4] = {false,false,false,false};  // decoded on THIS attempt

		for (int l = 0; l < R; l++)
		{
			if (links[self_idx].per_layer_done[l]) continue;
			Real coin   = randnum.u();
			Real bler_l;
			if (g_harq_ir == 1)
			{
				// IR/Chase soft combining: accumulate this (re)transmission's effective SINR
				// and look up BLER at the ACCUMULATED SINR. On the 1st tx the accumulator
				// equals the single-tx ESINR (→ identical to no-combining); each retx adds
				// more SINR so the per-layer BLER drops and decoding becomes progressively
				// more likely. per_layer_ESINR[l] is set by BLER_Mapping for this slot.
				// NOTE: _info_bits_per_layer is frozen during retx (recomputed only on
				// ACK), so the segmentation (BGN/Zc/C) is stable across combining. The
				// scheduler may still refresh _mcs_per_layer mid-HARQ — same semantics
				// as the legacy _cqi_per_layer lookup (future hardening: freeze at 1st tx).
				links[self_idx].per_layer_accum_esinr[l] += per_layer_ESINR[l];
				if (g_matlab_bler)
					bler_l = bler_lookup_matlab(links[self_idx].per_layer_accum_esinr[l],
					                            _mcs_per_layer[l], _info_bits_per_layer[l]);
				else
					bler_l = bler_lookup_dB(links[self_idx].per_layer_accum_esinr[l], _cqi_per_layer[l]);
			}
			else
			{
				bler_l = per_layer_BLER[l];
			}
			if (bler_l < coin || isnan(bler_l))          // layer decoded OK
			{
				links[self_idx]._throughput += _info_bits_per_layer[l];
				links[self_idx].per_layer_done[l] = true;
				layer_ok_first[l] = true;
			}
			else                                          // layer NACK
			{
				links[self_idx].per_layer_num_re_tx[l]++;
				if (links[self_idx].per_layer_num_re_tx[l] > 3)
				{
					links[self_idx].per_layer_done[l] = true;   // drop after 3 retx
					links[self_idx].re_tx_failed++;
					re_tx_failed_cnt++;
				}
			}
		}

		bool all_done = true;
		for (int l = 0; l < R; l++)
			if (!links[self_idx].per_layer_done[l]) { all_done = false; break; }

		if (first_tx)
		{
			// OLLA targets PER-LAYER BLER: feed each layer's first-tx outcome into the
			// window. The single per-UE offset then converges to the per-layer target_bler
			// (NOT TB-level, which for R layers would over-shrink MCS by ~R×).
			bool all_ok_first = true;
			for (int l = 0; l < R; l++)
			{
				update_olla_window(layer_ok_first[l]);
				if (!layer_ok_first[l]) all_ok_first = false;
			}
			if (!all_ok_first) links[self_idx].failed++;   // TB first-tx failure (any layer NACK)
		}

		if (all_done)
		{
			links[self_idx].ACK = true;
			links[self_idx].num_re_tx = 0;
			for (int l = 0; l < 4; l++) { links[self_idx].per_layer_done[l] = false; links[self_idx].per_layer_num_re_tx[l] = 0; }
		}
		else
		{
			links[self_idx].ACK = false;
			links[self_idx].num_re_tx++;
			links[self_idx].total_num_re_tx++;
		}
		return;
	}

	if (prob < tossing_coin || isnan(prob)) // Tossing a biased coin
	{
		links[self_idx].ACK = true;
		links[self_idx]._throughput += _info_bits;

		// OLLA: Track first transmission result (ACK)
		if (links[self_idx].num_re_tx == 0)
			update_olla_window(true);

		links[self_idx].num_re_tx = 0;
	}
	else
	{
		links[self_idx].ACK = false;
		if (links[self_idx].num_re_tx == 0)
		{
			links[self_idx].failed++;

			// OLLA: Track first transmission result (NACK)
			update_olla_window(false);
		}

		// jhnoh 220921
		links[self_idx].num_re_tx++;
		links[self_idx].total_num_re_tx++;
		if (links[self_idx].num_re_tx > 3)
		{
			links[self_idx].num_re_tx = 0;
			links[self_idx].ACK = true;

			// jhnoh 220921
			links[self_idx].re_tx_failed++;
			re_tx_failed_cnt++;
		}
	}
}

Real MS::Return_initial_packet_error(void) const
{
	if (links[self_idx].total_num_tx == 0)
	{
		return 0;
	}
	else
	{
		return ((Real)links[self_idx].failed / (Real)links[self_idx].total_num_tx);
	}
}


void MS::Find_Allocated_Rbs_Mcs_mTRP( void )
{
	SCHEDULE_DECISION decision_of_interest;
	rb_indices.clear(); //vector <int> rb_indices

	int temp_mcs = 0;
	int temp_cqi = 0;
	Real temp_avr_sinr_fb = 0;

	// Search all RBs and streams for this UE's scheduled allocations
	for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
	{
		for (int stream = 0; stream < mx_ue_mumimo; stream++)
		{
			decision_of_interest = pppSchedule_Map[links[self_idx]._sector_in_control][freq_idx][stream];
			if (decision_of_interest.ue_selected == self_idx)
			{
				rb_indices.push_back(freq_idx);
				if (rb_indices.size() == 1)
				{
					temp_mcs  = decision_of_interest.mcs_selected;
					temp_cqi  = decision_of_interest.cqi_selected;
					temp_avr_sinr_fb = decision_of_interest.capacity;

					// PER-LAYER MCS: on the FIRST scheduled RB of a fresh TX (ACK==1), capture
					// every layer-stream's MCS/CQI by layer_idx. Gated by ACK so a retransmitted
					// TB keeps its frozen per-layer MCS (matches legacy _mcs_idx freezing).
					if (g_per_layer_mcs == 1 && links[self_idx].ACK)
					{
						for (int l = 0; l < 4; l++) { _mcs_per_layer[l] = 0; _cqi_per_layer[l] = 0; _est_sinr_per_layer[l] = 0; }
						for (int s2 = 0; s2 < mx_ue_mumimo; s2++)
						{
							const SCHEDULE_DECISION& d2 =
								pppSchedule_Map[links[self_idx]._sector_in_control][freq_idx][s2];
							if (d2.ue_selected != self_idx) continue;
							int L = d2.layer_idx;
							if (L >= 0 && L < 4)
							{
								_mcs_per_layer[L]      = d2.mcs_selected;
								_cqi_per_layer[L]      = d2.cqi_selected;
								_est_sinr_per_layer[L] = d2.capacity;   // scheduler's linear SINR estimate (pre-OLLA)
							}
						}
					}
				}
				break; // Found in this RB, move to next RB
			}
		}
	}

	//if ACK, temp_mcs from schedule map -> _mcs_idx
	if (links[self_idx].ACK)
	{
		if (_mcs_idx == -1)
		{
			cout << "error!!" << endl;
			getchar();
		}

		_mcs_idx = temp_mcs;
		_cqi_idx = temp_cqi;
		_avr_sinr = temp_avr_sinr_fb;
	}
}

/*===================================================================
FUNCTION: UE::Find_allocated_rbs( void )

AUTHOR: DH

DESCRIPTION: UE checks the schedule map, and if allocated rbs is found,
add rbs index to vector ('rb_indices') and check the mcs type.

NOTES: Only used by UE::receiveUL( void )

INPUT: schedule_map,
OUTPUT: rb_indices
===================================================================*/
void MS::Find_Allocated_Rbs_Mcs(int self_ms_idx)
{
	int temp_ms_idx = self_ms_idx; //ms_idx (0~num_MS) from main function
	SCHEDULE_DECISION decision_of_interest;
	rb_indices.clear(); //vector <int> rb_indices

	int temp_mcs = 0;
	int temp_cqi = 0;

	for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
	{
		for (int stream = 0; stream < mx_ue_mumimo; stream++)
		{
			decision_of_interest = pppSchedule_Map[links[temp_ms_idx]._sector_in_control][freq_idx][stream];			
			if (decision_of_interest.ue_selected == temp_ms_idx) 
			{
				rb_indices.push_back(freq_idx);
				if (rb_indices.size() == 1)
				{
					if (temp_mcs == -1)
					{
						cout << "error!!" << endl;
						getchar();
					}

					temp_mcs  = decision_of_interest.mcs_selected;
					temp_cqi  = decision_of_interest.cqi_selected;
				}
				else
				{
					if (temp_mcs != decision_of_interest.mcs_selected)
					{
						cerr << "Something Wrong in the Scheduler" << endl;
						getchar();
					}

					if (temp_cqi != decision_of_interest.cqi_selected)
					{
						cerr << "Something Wrong in the Scheduler CQI" << endl;
						getchar();
					}
				}	
			}
		}
	}

	if (links[temp_ms_idx].ACK)
	{
		if (_mcs_idx == -1)
		{
			cout << "error!!" << endl;
			getchar();
		}

		_mcs_idx = temp_mcs;
		_cqi_idx = temp_cqi;
	}
}


/*===================================================================
FUNCTION: UE::Compute_num_bin_symbol( void )

NOTES: 2009-12-30

===================================================================*/
void MS::Rate_matching(void)
{
	if (links[receive_self_ms_idx].ACK)
	{
		int _num_traffic_ = 0;
		int _info_bits_ = 0;

		int _mod_type_;

		int N_DMRS = 0;
		int N_overhead = 0;

		int N_RE_a;
		int N_RE;
		int N_info;

		N_RE_a = num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms - N_DMRS - N_overhead;

		N_RE = MIN(156, N_RE_a) * rb_indices.size();


		Real coding_rate_x1024 = 0;
		Real coding_rate;
		int TB_Size_bits = 0;

		// MCS -> (Qm, CR) via the shared table (was a duplicated if/else chain)
		nr_mcs_to_qm_r(_mcs_idx, _mod_type_, coding_rate_x1024);

		coding_rate = coding_rate_x1024 / 1024.;

		// BUG FIX 2026-05-25: TBS must scale with the ACTUAL number of layers
		// allocated to this UE (set in compute_tone_SINR()), NOT the fixed cfg
		// constant NUM_UE_Layer. Falls back to NUM_UE_Layer when no layer was
		// allocated (defensive — shouldn't happen if UE is scheduled).
		const int _R_used_ = (num_layers_actual > 0) ? num_layers_actual : NUM_UE_Layer;

		N_info = (int)(N_RE * coding_rate * _mod_type * _R_used_ * (1. - overhead));

		if (N_info <= 3834)
		{
			int n = MAX(3, floor(log2(N_info)) - 6.);

			int N_info_a = MAX(24, pow(2, n) * floor(N_info / pow(2, n)));

			_info_bits_ = N_info_a;
		}
		else
		{
			int n = floor(log2(N_info - 24)) - 5;

			int N_info_a = MAX(3840, pow(2, n) * round(((N_info - 24) / pow(2, n))));

			if (coding_rate <= (1. / 4.))
			{
				int C = ceil((N_info_a + 24) / 3816);

				_info_bits_ = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
			}
			else
			{
				if (N_info_a >= 8424)
				{
					int C = ceil((N_info_a + 24) / 8424);

					_info_bits_ = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
				}
				else
				{
					_info_bits_ = 8 * ceil((N_info_a + 24) / 8) - 24;
				}
			}

		}


		if (((Real)_info_bits_ / (Real)_num_traffic_) > 0.93)
		{
			//_mcs_idx = _mcs_idx - 1;
		}
				
	}
}

void MS::Compute_RBs_SINR( void )
{
	// initial condition
	if (links[self_idx].ACK == 1) // Previously successed or initial transmission case
	{
		rbs_rx.clear();
		comp_mode_rx_flag = 0;
		num_rx_rbs = 0;
		// Also reset per-layer buffers for this packet's reception
		for (auto& v : rbs_rx_per_layer) v.clear();
		rbs_rx_per_layer.clear();
		// Reset per-layer HARQ state for the fresh TB (per_layer_mcs path)
		for (int l = 0; l < 4; l++) { links[self_idx].per_layer_done[l] = false; links[self_idx].per_layer_num_re_tx[l] = 0; links[self_idx].per_layer_accum_esinr[l] = 0.0; }
	}

	for (int freq_idx = 0; freq_idx < (int)rb_indices.size(); freq_idx++)
	{
		Real rx_sinr = 0;
		if (Scheduling_Type == 0)
		{
			rx_sinr = compute_tone_SINR_MRC(rb_indices[freq_idx]);
		}
		else if(Scheduling_Type == 1)
		{
			rx_sinr = compute_tone_SINR(rb_indices[freq_idx]);
		}
		else if(Scheduling_Type == 2)
		{
			rx_sinr = Compute_Tone_SINR_NCJT(links[self_idx]._sector_in_control, freq_idx);
		}
		else if(Scheduling_Type == 3)
		{
			rx_sinr = Compute_Tone_SINR_CJT( links[self_idx]._sector_in_control, freq_idx);
		}
	
		if (rx_sinr)
		{
			rbs_rx.push_back(rx_sinr);
		}

	}

	if (links[self_idx].ACK == 1) // Previously successed or initial transmission case
	{
		num_rx_rbs = rbs_rx.size(); // set the number of resource blocks for TBS calculation
	}

	Real sinr = 0;
	for (int rb_idx = 0; rb_idx < (int)rbs_rx.size(); rb_idx++)
	{
		sinr += rbs_rx[rb_idx];
	}
	received_sinr_avg = sinr / (Real)num_rb;
}

void MS::Compute_Transport_Block_Size( void )
{
	if (links[self_idx].ACK == 1) // Previously successed or initial transmission case
	{
		// PER-LAYER TBS: each layer gets its own sub-TBS from its own MCS (layer factor = 1);
		// total info bits = sum over layers. Sum ~= legacy (x R) when all per-layer MCS equal.
		if (g_per_layer_mcs == 1)
		{
			int N_DMRS = 0, N_overhead = 0;
			int N_RE_a = num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms - N_DMRS - N_overhead;
			int N_RE   = MIN(156, N_RE_a) * num_rx_rbs;

			int R = (num_layers_actual > 0) ? num_layers_actual : 1;
			if (R > 4) R = 4;
			int total = 0;
			for (int l = 0; l < 4; l++) _info_bits_per_layer[l] = 0;
			for (int l = 0; l < R; l++)
			{
				_info_bits_per_layer[l] = tbs_info_bits_from_mcs(_mcs_per_layer[l], N_RE, 1);
				total += _info_bits_per_layer[l];
			}
			_info_bits = total;

			// Keep _mod_type meaningful for legacy stat readers (dominant layer).
			if      (_mcs_per_layer[0] <= 4)  _mod_type = QPSK;
			else if (_mcs_per_layer[0] <= 10) _mod_type = QAM16;
			else if (_mcs_per_layer[0] <= 19) _mod_type = QAM64;
			else                              _mod_type = QAM256;
			return;
		}

		int N_DMRS = 0;
		int N_overhead = 0;

		int N_RE_a;
		int N_RE;
		int N_info;

		N_RE_a = num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms - N_DMRS - N_overhead;
		N_RE = MIN(156, N_RE_a) * num_rx_rbs; //rb_indices.size();

		Real coding_rate_x1024 = 0;
		Real coding_rate;
		int TB_Size_bits = 0;

		// MCS -> (Qm, CR) via the shared table (was a duplicated if/else chain)
		nr_mcs_to_qm_r(_mcs_idx, _mod_type, coding_rate_x1024);

		coding_rate = coding_rate_x1024 / 1024.;
		//N_info = (int)(N_RE * coding_rate * _mod_type * NUM_UE_Layer* (1. - overhead));
		// BUG FIX 2026-05-25: use actual per-UE layer count.
		const int _R_used_tbs = (num_layers_actual > 0) ? num_layers_actual : 1;
		N_info = (int)(N_RE * coding_rate * _mod_type * _R_used_tbs * (1. - overhead));
		int N_info_a = 0;
		//N_info = (int)(N_RE * coding_rate * _mod_type * 1 * (1. - overhead));
		//N_info = (int)(N_RE * coding_rate * _mod_type * mx_ue_mumimo);

		if (N_info <= 3834)
		{
			int n = MAX(3, floor(log2(N_info)) - 6.);
			int N_info_a = MAX(24, pow(2, n) * floor(N_info / pow(2, n)));
			_info_bits = N_info_a;
		}
		else
		{
			int n = floor(log2(N_info - 24)) - 5;
			N_info_a = MAX(3840, pow(2, n) * round(((N_info - 24) / pow(2, n))));
			//int N_info_a = pow(2, n) * round(((N_info - 24) / pow(2, n)));

			if (coding_rate <= (1. / 4.))
			{
				int C = ceil((N_info_a + 24) / 3816);
				_info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
			}
			else
			{
				if (N_info_a >= 8424)
				{
					int C = ceil((N_info_a + 24) / 8424);
					_info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
					
				}
				else
				{
					_info_bits = 8 * ceil((N_info_a + 24) / 8) - 24;
				}
			}

		}		
	}
}

/*===================================================================
FUNCTION: UE::Compute_num_bin_symbol( void )

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void MS::Compute_Num_Bin_Symbol(void)
{
	// initial condition
	num_rx_rb = 0;
	if (links[receive_self_ms_idx].ACK == 0) // Previously failed case
	{
		for (int freq_idx = 0; freq_idx < (int)rb_indices.size(); freq_idx++)
		{
			Real rx_sinr = 0;
			if (Scheduling_Type == 0) // roundrobin 1x2
			{
				rx_sinr = compute_tone_SINR_MRC(rb_indices[freq_idx]); //
			}
			else if(Scheduling_Type == 1)
			{
				rx_sinr = compute_tone_SINR(rb_indices[freq_idx]); //
			}
			else if(Scheduling_Type == 2)
			{
				rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx]._sector_in_control, freq_idx);
			}

			if (rx_sinr)
			{
				rbs_rx.push_back(rx_sinr);
			}

			}
		_num_traffic = (int)(num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms * rbs_rx.size() * (1. - overhead) * _mod_type);
	}
	else // Previously successed or initial transmission case
	{
		rbs_rx.clear();
		comp_mode_rx_flag = 0;

		for (int freq_idx = 0; freq_idx < (int)rb_indices.size(); freq_idx++)
		{
			Real rx_sinr = 0;
			if (Scheduling_Type == 0)
			{
				rx_sinr = compute_tone_SINR_MRC(rb_indices[freq_idx]);
			}
			else if(Scheduling_Type == 1)
			{
				rx_sinr = compute_tone_SINR(rb_indices[freq_idx]);
			}
			else if(Scheduling_Type == 2)
			{
				rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx]._sector_in_control, freq_idx);
			}

			if (rx_sinr)
			{
				rbs_rx.push_back(rx_sinr);
			}

			}

		if (_mcs_idx <= 4)       { _mod_type = QPSK; }
		else if (_mcs_idx <= 10) { _mod_type = QAM16; }
		else if (_mcs_idx <= 19) { _mod_type = QAM64; }
		else if (_mcs_idx <= 28) { _mod_type = QAM256; }
		else                     { cerr << "Something Wrong - rate matching" << endl; }

		int N_DMRS = 0;
		int N_overhead = 0;

		int N_RE_a;
		int N_RE;
		int N_info;

		N_RE_a = num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms - N_DMRS - N_overhead;
		N_RE = MIN(156, N_RE_a) * rb_indices.size();

		Real coding_rate_x1024 = 0;
		Real coding_rate;
		int TB_Size_bits = 0;

		///////
		if (_mcs_idx == 0)        { coding_rate_x1024 = 120; }
		else if (_mcs_idx == 1)   { coding_rate_x1024 = 193; }
		else if (_mcs_idx == 2)   { coding_rate_x1024 = 308; }
		else if (_mcs_idx == 3)   { coding_rate_x1024 = 449; }
		else if (_mcs_idx == 4)   { coding_rate_x1024 = 602; }
		else if (_mcs_idx == 5)   { coding_rate_x1024 = 378; }
		else if (_mcs_idx == 6)   { coding_rate_x1024 = 434; }
		else if (_mcs_idx == 7)   { coding_rate_x1024 = 490; }
		else if (_mcs_idx == 8)   { coding_rate_x1024 = 553; }
		else if (_mcs_idx == 9)   { coding_rate_x1024 = 616; }
		else if (_mcs_idx == 10)  { coding_rate_x1024 = 658; }
		else if (_mcs_idx == 11)  { coding_rate_x1024 = 466; }
		else if (_mcs_idx == 12)  { coding_rate_x1024 = 517; }
		else if (_mcs_idx == 13)  { coding_rate_x1024 = 567; }
		else if (_mcs_idx == 14)  { coding_rate_x1024 = 616; }
		else if (_mcs_idx == 15)  { coding_rate_x1024 = 666; }
		else if (_mcs_idx == 16)   { coding_rate_x1024 = 719; }
		else if (_mcs_idx == 17)   { coding_rate_x1024 = 772; }
		else if (_mcs_idx == 18)   { coding_rate_x1024 = 822; }
		else if (_mcs_idx == 19)   { coding_rate_x1024 = 873; }
		else if (_mcs_idx == 20)   { coding_rate_x1024 = 682.5; }
		else if (_mcs_idx == 21)   { coding_rate_x1024 = 711; }
		else if (_mcs_idx == 22)   { coding_rate_x1024 = 754; }
		else if (_mcs_idx == 23)   { coding_rate_x1024 = 797; }
		else if (_mcs_idx == 24)   { coding_rate_x1024 = 841; }
		else if (_mcs_idx == 25)  { coding_rate_x1024 = 885; }
		else if (_mcs_idx == 26)  { coding_rate_x1024 = 916.5; }
		else if (_mcs_idx == 27)  { coding_rate_x1024 = 948; }
		else
		{
			cout << "Something Wrong rate matching" << endl;
			getchar();
		}

		coding_rate = coding_rate_x1024 / 1024.;
		// BUG FIX 2026-05-25: use actual per-UE layer count instead of cfg constant.
		const int _R_used_cb = (num_layers_actual > 0) ? num_layers_actual : NUM_UE_Layer;
		N_info = (int)(N_RE * coding_rate * _mod_type * _R_used_cb * (1. - overhead));
		//N_info = (int)(N_RE * coding_rate * _mod_type * 1 * (1. - overhead));
		//N_info = (int)(N_RE * coding_rate * _mod_type * mx_ue_mumimo);

		if (N_info <= 3834)
		{
			int n = MAX(3, floor(log2(N_info)) - 6.);
			int N_info_a = MAX(24, pow(2, n) * floor(N_info / pow(2, n)));
			_info_bits = N_info_a;
		}
		else
		{
			int n = floor(log2(N_info - 24)) - 5;
			int N_info_a = MAX(3840, pow(2, n) * round(((N_info - 24) / pow(2, n))));

			if (coding_rate <= (1. / 4.))
			{
				int C = ceil((N_info_a + 24) / 3816);
				_info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
			}
			else
			{
				if (N_info_a >= 8424)
				{
					int C = ceil((N_info_a + 24) / 8424);
					_info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
				}
				else
				{
					_info_bits = 8 * ceil((N_info_a + 24) / 8) - 24;
				}
			}

		}
		_num_traffic = (int)(num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms * rbs_rx.size() * (1. - overhead) * _mod_type);

	}

	Real sinr = 0;
	for (int rb_idx = 0; rb_idx < (int)rbs_rx.size(); rb_idx++)
	{
		sinr += rbs_rx[rb_idx];
	}
	
	received_sinr_avg = 10. * log10(sinr / (Real)num_rb);     ///// dB scale

	//Averaging_Sinr(_self_idx, sinr / (Real)num_rb);
	Averaging_Sinr(receive_self_ms_idx, sinr / (Real)num_rb); 
}

int getIndex(vector<int> v, int K);

int determine_MCS(Real sinr_estimate)
{
	// Throughput-maximizing selection over all 28 MCS (TBS/segmentation aware),
	// precomputed as a dB-grid at init. Bypasses the CQI->MCS map (which reaches
	// only 15 of 28 MCS). sinr_estimate is in dB at every call site.
	if (g_matlab_bler && g_matlab_tput_mcs && TputMCS_Grid_ready())
		return nr_mcs_maxtput_from_dB(sinr_estimate);

	if (false)  //// TS 38.214 , Table 5.1.3.1-2:, use spectral efficiency
	{
		int _MCS_decision = 0;


		if (sinr_estimate < 0.2344)       { _MCS_decision = 0; }
		else if (sinr_estimate < 0.3770)  { _MCS_decision = 1; }
		else if (sinr_estimate < 0.6016)  { _MCS_decision = 2; }
		else if (sinr_estimate < 0.8770)  { _MCS_decision = 3; }
		else if (sinr_estimate < 1.1758)  { _MCS_decision = 4; }
		else if (sinr_estimate < 1.4766)  { _MCS_decision = 5; }
		else if (sinr_estimate < 1.6953)  { _MCS_decision = 6; }
		else if (sinr_estimate < 1.9141)  { _MCS_decision = 7; }
		else if (sinr_estimate < 2.1602)  { _MCS_decision = 8; }
		else if (sinr_estimate < 2.4063)  { _MCS_decision = 9; }
		else if (sinr_estimate < 2.5703)  { _MCS_decision = 10; }
		else if (sinr_estimate < 2.7305)  { _MCS_decision = 11; }
		else if (sinr_estimate < 3.0293)  { _MCS_decision = 12; }
		else if (sinr_estimate < 3.3223)  { _MCS_decision = 13; }
		else if (sinr_estimate < 3.6094)  { _MCS_decision = 14; }
		else if (sinr_estimate < 3.9023)  { _MCS_decision = 15; }
		else if (sinr_estimate < 4.2129)  { _MCS_decision = 16; }
		else if (sinr_estimate < 4.5234)  { _MCS_decision = 17; }
		else if (sinr_estimate < 4.8164)  { _MCS_decision = 18; }
		else if (sinr_estimate < 5.1152)  { _MCS_decision = 19; }
		else if (sinr_estimate < 5.3320)  { _MCS_decision = 20; }
		else if (sinr_estimate < 5.5547)  { _MCS_decision = 21; }
		else if (sinr_estimate < 5.8906)  { _MCS_decision = 22; }
		else if (sinr_estimate < 6.2266)  { _MCS_decision = 23; }
		else if (sinr_estimate < 6.5703)  { _MCS_decision = 24; }
		else if (sinr_estimate < 6.9141)  { _MCS_decision = 25; }
		else if (sinr_estimate < 7.1602)  { _MCS_decision = 26; }
		else if (sinr_estimate < 7.4063)  { _MCS_decision = 27; }
		else                              { _MCS_decision = 27; }


		return _MCS_decision;
	}

	else if (false)  /// CQI
	{
		int _CQI_decision = 0;


		if (sinr_estimate < -9.478)  { _CQI_decision = 1; }
		else if (sinr_estimate < -6.658)  { _CQI_decision = 2; }
		else if (sinr_estimate < -4.098)  { _CQI_decision = 3; }
		else if (sinr_estimate < -1.798)  { _CQI_decision = 4; }
		else if (sinr_estimate < 0.399)   { _CQI_decision = 5; }
		else if (sinr_estimate < 2.424)   { _CQI_decision = 6; }
		else if (sinr_estimate < 4.489)   { _CQI_decision = 7; }
		else if (sinr_estimate < 6.367)   { _CQI_decision = 8; }
		else if (sinr_estimate < 8.456)   { _CQI_decision = 9; }
		else if (sinr_estimate < 10.266)  { _CQI_decision = 10; }
		else if (sinr_estimate < 12.218)  { _CQI_decision = 11; }
		else if (sinr_estimate < 14.122)  { _CQI_decision = 12; }
		else if (sinr_estimate < 15.849)  { _CQI_decision = 13; }
		else if (sinr_estimate < 17.786)  { _CQI_decision = 14; }
		else                              { _CQI_decision = 15; }


		return _CQI_decision;

	}

	else if (true)  
	{
		int _CQI_decision = 0;

		int _MCS_decision = 0;
			/// VIENNA curve

		if (sinr_estimate > SINR_threshold_dB[14])       { _CQI_decision = 15; }
		else if (sinr_estimate > SINR_threshold_dB[13])  { _CQI_decision = 14; }
		else if (sinr_estimate > SINR_threshold_dB[12])  { _CQI_decision = 13; }
		else if (sinr_estimate > SINR_threshold_dB[11])  { _CQI_decision = 12; }
		else if (sinr_estimate > SINR_threshold_dB[10])  { _CQI_decision = 11; }
		else if (sinr_estimate > SINR_threshold_dB[9])  { _CQI_decision = 10; }
		else if (sinr_estimate > SINR_threshold_dB[8])  { _CQI_decision = 9; }
		else if (sinr_estimate > SINR_threshold_dB[7])  { _CQI_decision = 8; }
		else if (sinr_estimate > SINR_threshold_dB[6])  { _CQI_decision = 7; }
		else if (sinr_estimate > SINR_threshold_dB[5])  { _CQI_decision = 6; }
		else if (sinr_estimate > SINR_threshold_dB[4])  { _CQI_decision = 5; }
		else if (sinr_estimate > SINR_threshold_dB[3])  { _CQI_decision = 4; }
		else if (sinr_estimate > SINR_threshold_dB[2])  { _CQI_decision = 3; }
		else if (sinr_estimate > SINR_threshold_dB[1])  { _CQI_decision = 2; }
		else if (sinr_estimate > SINR_threshold_dB[0])  { _CQI_decision = 1; }
		else                                            { _CQI_decision = 1; }

			//_CQI_decision = 8;
			// TS 36.101  TABLE A.4-7 �� MCS Scheme 1
			if      (_CQI_decision == 1)   { _MCS_decision = 0;  }
			else if (_CQI_decision == 2)   { _MCS_decision = 0;  }
			else if (_CQI_decision == 3)   { _MCS_decision = 2;  }
			else if (_CQI_decision == 4)   { _MCS_decision = 4;  }
			else if (_CQI_decision == 5)   { _MCS_decision = 6;  }
			else if (_CQI_decision == 6)   { _MCS_decision = 8;  }
			else if (_CQI_decision == 7)   { _MCS_decision = 11; }
			else if (_CQI_decision == 8)   { _MCS_decision = 13; }
			else if (_CQI_decision == 9)   { _MCS_decision = 16; }
			else if (_CQI_decision == 10)  { _MCS_decision = 18; }
			else if (_CQI_decision == 11)  { _MCS_decision = 21; }
			else if (_CQI_decision == 12)  { _MCS_decision = 23; }
			else if (_CQI_decision == 13)  { _MCS_decision = 25; }
			else if (_CQI_decision == 14)  { _MCS_decision = 27; }
			else if (_CQI_decision == 15)  { _MCS_decision = 27; }



			return _MCS_decision;

	
	}
	else  //  old code
	{
	}


	
}





int determine_CQI(Real sinr_estimate)
{

	if (true)  /// CQI
	{
		int _CQI_decision = 0;

		if (sinr_estimate > SINR_threshold_dB[14])       { _CQI_decision = 15; }
		else if (sinr_estimate > SINR_threshold_dB[13])  { _CQI_decision = 14; }
		else if (sinr_estimate > SINR_threshold_dB[12])  { _CQI_decision = 13; }
		else if (sinr_estimate > SINR_threshold_dB[11])  { _CQI_decision = 12; }
		else if (sinr_estimate > SINR_threshold_dB[10])  { _CQI_decision = 11; }
		else if (sinr_estimate > SINR_threshold_dB[9])  { _CQI_decision = 10; }
		else if (sinr_estimate > SINR_threshold_dB[8])  { _CQI_decision = 9; }
		else if (sinr_estimate > SINR_threshold_dB[7])  { _CQI_decision = 8; }
		else if (sinr_estimate > SINR_threshold_dB[6])  { _CQI_decision = 7; }
		else if (sinr_estimate > SINR_threshold_dB[5])  { _CQI_decision = 6; }
		else if (sinr_estimate > SINR_threshold_dB[4])  { _CQI_decision = 5; }
		else if (sinr_estimate > SINR_threshold_dB[3])  { _CQI_decision = 4; }
		else if (sinr_estimate > SINR_threshold_dB[2])  { _CQI_decision = 3; }
		else if (sinr_estimate > SINR_threshold_dB[1])  { _CQI_decision = 2; }
		else if (sinr_estimate > SINR_threshold_dB[0])  { _CQI_decision = 1; }
		else                                            { _CQI_decision = 1; }
		return _CQI_decision;
	}
}


/*===================================================================
FUNCTION: UE::Compute_Code_Rate( void )

AUTHOR: SJ, KR, DW

DESCRIPTION: Compute the code rate

NOTES:

===================================================================*/
void MS::Compute_Effective_SINR(void)
{
	// compute ESINR
	// EESM method
	Real e_sir = 0.;
	for (int i = 0; i < rbs_rx.size(); i++) {
		e_sir += exp(-rbs_rx[i] / Beta[_cqi_idx-1]);
	}

	e_sir = -Beta[_cqi_idx-1] * log(e_sir / (Real)(rbs_rx.size()));	

	int determined_mcs_test;
    int cqi_test_result;
	int determined_cqi;

	if (0) // ( self_idx == 1)  // Disabled to avoid disrupting progress bar
	{
		cout << " e_sir = " << 10*log10(e_sir) << endl;
	}

	// Compute esinr again with 15 cqi_index
	if (isinf(e_sir) && _cqi_idx < 14)
	{
		e_sir = 0.;
		for (int i = 0; i < rbs_rx.size(); i++) {
			e_sir += exp(-rbs_rx[i] / Beta[14]);
		}
		e_sir = -Beta[14] * log(e_sir / (Real)(rbs_rx.size()));	
	}

	// Compute esinr with linear average
	if (isinf(e_sir))
	{
		e_sir = 0.;
		for (int i = 0; i < rbs_rx.size(); i++) {
			e_sir += rbs_rx[i];		
		}
		e_sir = e_sir/rbs_rx.size();
		ESINR_linear = e_sir;
	} 
	else
	{
		ESINR_linear = e_sir;
	}
}

/*===================================================================
FUNCTION: UE::return_throughput( void )

AUTHOR: SJ

DESCRIPTION: return the total throughput over the frequency

NOTES:

===================================================================*/
INPUT_TO_LLS_MAPPING MS::Return_mapping_info(void)
{
	INPUT_TO_LLS_MAPPING input;
	input.mod_type = _mod_type;
	input.code_rate = _code_rate;
	input.snr_scale_factor = _snr_scale_factor;

	return (input);
}

Real MS::Compute_Tone_SINR_CJT(int sec_idx, int rb_idx)
{
	MatrixXcReal w, u, h_u;
	MatrixXcReal H_bar1;
	MatrixXcReal H_bar2;
	MatrixXcReal H_bar;
	MatrixXcReal H_inf;
	MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
	Real A, B, C;	
	Real Received_SINR;

	//int sec_idx = links[self_idx].adj_sector[0];
	int bs_idx = (int)(sec_idx / 3.);
	Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);

	int comp_sector_idx = links[self_idx].adj_sector[1];

	Real linear_signal       = dBm2linear(links[self_idx].str_signal); // dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
	Real linear_interference = dBm2linear(links[self_idx].interference);
	linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);

	// Use dominant layer (col 0) for mTRP CJT SINR (Phase I: single-layer)
	H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length].col(0);

	Real comp_sector_static_gain = dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]); //dBm2linear(links[self_idx].static_gain[1].first); // dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
	H_inf = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length].col(0);

	u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + H_inf * H_inf.adjoint() + (linear_interference + noise) * Identity).inverse();
	A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
	B = (u * H_inf * H_inf.adjoint() * u.adjoint()).norm();
	C = (linear_interference + noise) * (u * u.adjoint()).norm();			
	Received_SINR = A/(B+C);
	return Received_SINR;
}

Real test111;
Real MS::Compute_Tone_SINR_NCJT(int sec_idx, int rb_idx)
{
	Real Received_SINR = 0;
	MatrixXcReal w, u, h_u;
	MatrixXcReal H_bar;
	MatrixXcReal H_inf;
	MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
	Real A,B,C;

	if ( sec_idx == links[self_idx]._sector_in_control )
	{
		if ( ppSchedulerRead[rb_idx][0].ue_selected == self_idx )
		{
				num_rx_rb++;

			int    bs_idx              = (int) sec_idx/3.;
			Real noise               = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
			// jhnoh 220901
			Real linear_signal       = dBm2linear(links[self_idx].str_signal);//dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
			Real linear_interference = dBm2linear(links[self_idx].interference);

			H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * sector[links[self_idx]._sector_in_control].W[rb_idx];

			//int comp_sector_idx = links[self_idx].adj_sector[1];
			int comp_sector_idx = links[self_idx].comp_sector_idx;

			int comp_bs_idx = (int) comp_sector_idx/3.;
			SCHEDULE_DECISION comp_sector_decision;
			comp_sector_decision = pppSchedule_Map[comp_sector_idx][rb_idx][0];

			if ( comp_sector_decision.ue_selected == self_idx ) 
			{
				// jhnoh 220901
				//Real comp_sector_static_gain = dBm2linear(links[self_idx].static_gain[1].first);//dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
				Real comp_sector_static_gain = dBm2linear(links[self_idx].comp_interf_strength[1]);
				//H_inf = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * sector[comp_sector_idx].W[rb_idx];

				if ( comp_bs_idx < num_BS )
				{
					if (links[self_idx].intf_w_rnd_RSRP[comp_sector_idx] != 0)
						linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);
				}
				else
				{
					if ((3*comp_bs_idx)   != links[self_idx]._sector_in_control )
						linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*comp_bs_idx]);
					if ((3*comp_bs_idx+1) != links[self_idx]._sector_in_control )
						linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*comp_bs_idx+1]);
					if ((3*comp_bs_idx+2) != links[self_idx]._sector_in_control )
						linear_interference -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*comp_bs_idx+2]);
				}
				linear_interference += comp_sector_static_gain;
				u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise)* Identity).inverse();
				A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
				C = (linear_interference + noise) * (u * u.adjoint()).norm();

				if (C < 0)
					cout << " something wrong with C 0" << endl;

			}	
			else
			{
				u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise)* Identity).inverse();
				A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
				C = (linear_interference + noise) * (u * u.adjoint()).norm();
			}	
			Received_SINR = A/C;

			if (Received_SINR < 0)
				cout << " something wrong with SINR 00" << endl;
		}
	}
	else if ( sec_idx == links[self_idx].comp_sector_idx)
	{
		SCHEDULE_DECISION comp_sector_decision;
		comp_sector_decision = pppSchedule_Map[sec_idx][rb_idx][0];
		if ( comp_sector_decision.ue_selected == self_idx )
		{
			num_rx_rb++;
			
			int bs_idx                 = (int)sec_idx/3.;
			Real noise               = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);

			// jhnoh 220901
			Real linear_signal       = dBm2linear(links[self_idx].comp_signal_strength[1]);
			Real linear_interference = dBm2linear(links[self_idx].interference);
			
			if ( bs_idx < num_BS )
			{
				if (links[self_idx].intf_w_rnd_RSRP[sec_idx] != 0)
					linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[sec_idx]);
			}
			else
			{
				if (links[self_idx].intf_w_rnd_RSRP[3*bs_idx] != 0)
					linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx  ]);
				
				if (links[self_idx].intf_w_rnd_RSRP[3*bs_idx+1] != 0)
					linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx+1]);
				
				if (links[self_idx].intf_w_rnd_RSRP[3*bs_idx+2] != 0)
					linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[3*bs_idx+2]);
			}
			linear_interference       += dBm2linear(links[self_idx].comp_interf_strength[0]);

			H_bar = sqrt(linear_signal) * (H_m[1][rb_idx]) * sector[sec_idx].W[rb_idx];

			SCHEDULE_DECISION comp_sector_decision;
			int comp_sector_idx  = links[self_idx]._sector_in_control;
			int comp_bs_idx      = (int)comp_sector_idx/3.;
			comp_sector_decision = pppSchedule_Map[comp_sector_idx][rb_idx][0];

			if ( comp_sector_decision.ue_selected == self_idx ) 
			{
				u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise)* Identity).inverse();
				A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
				C = (linear_interference + noise) * (u * u.adjoint()).norm();
			}	
			else
			{
				cerr << "Compute_Tone_SINR_NCJT : UE is only scheduled in an adjacent sector?? " << endl;
			}	
			Received_SINR = A/C;

			if (Received_SINR < 0)
				cout << " something wrong with SINR 1" << endl;

			Real temp_sinr = 0;
		}
	}

	if (Received_SINR < 0)
		cout << " something wrong with SINR 2" << endl;

	return Received_SINR;
}


/*===================================================================
FUNCTION: UE::compute_tone_SINR( )

DESCRIPTION:
===================================================================*/
Real MS::compute_tone_SINR(int rb_idx)
{
	// Collect ALL streams for this UE (multi-layer support: eType II rank > 1)
	std::vector<int> my_streams;
	for (int i = 0; i < mx_ue_mumimo; i++)
	{
		if (ppSchedulerRead[rb_idx][i].ue_selected == self_idx)
		{
			my_streams.push_back(i);
		}
	}

	// If UE not scheduled on this RB, return 0
	if (my_streams.empty())
		return 0.0;

	// Count received RBs (once per RB, not per layer)
	num_rx_rb++;

	// Record actual layers allocated to this UE on this RB (used by Rate_matching for TBS)
	num_layers_actual = (int)my_streams.size();

	// Calculate power levels
	const Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
	const Real linear_signal = dBm2linear(bs_maxpower - channel[(int)(links[self_idx]._sector_in_control/3)][self_idx].pathloss_final);
	const Real linear_interference = dBm2linear(links[self_idx].interference);

	// Effective channel matrix H_bar = sqrt(P) * H * W  (N_rx × mx_ue_mumimo)
	const MatrixXcReal& H = H_m[0][rb_idx];
	const MatrixXcReal& W = sector[links[self_idx]._sector_in_control].W[rb_idx];
	const MatrixXcReal H_bar = sqrt(linear_signal) * H * W;

	const MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);

	// --- Per-layer SINR computation -----------------------------------------------
	// For each of own UE's R layers, compute its post-equalizer SINR independently.
	// We do NOT sum into a single SINR (which over-counted in the previous version);
	// instead we store per-layer SINRs so that BLER lookup can be done per-layer.
	//
	// Use_SIC modes:
	//   0: per-stream MMSE — Ryy includes ALL streams. Each layer treats every other
	//      stream (own & other UE) as interference. Most pessimistic.
	//   1: legacy "hybrid" — Ryy w/ own layers but own excluded from intf sum.
	//   2: parallel ideal SIC — Ryy excludes own UE's layers entirely. Each layer's
	//      filter sees only other UE+noise covariance. Information-theoretic bound.
	//   3: sequential MMSE-SIC — order layers by SINR, decode strongest first,
	//      cancel it from Ryy if decode succeeded (BLER trial), then proceed.
	//      Most realistic.

	const int R_act = (int)my_streams.size();
	const int Rcap = (R_act < 4) ? R_act : 4;

	// Resize per-layer buffer at first RB of the alloc; appending across RBs.
	if (rbs_rx_per_layer.size() != (size_t)R_act)
		rbs_rx_per_layer.assign(R_act, std::vector<Real>());

	// Initial Ryy:
	//   SIC=2 (parallel ideal):  exclude own UE's layers from the start
	//   SIC=3 (sequential SIC):  start with ALL streams in Ryy. We then
	//                            sequentially subtract own UE's layers as
	//                            they are "decoded" (BLER-trial based).
	//   SIC=0/1 (legacy):        all streams in Ryy
	MatrixXcReal Ryy_base;
	if (g_use_sic == 2) {
		Ryy_base = linear_interference * Identity;
		for (int s = 0; s < mx_ue_mumimo; s++) {
			int sched_ue = ppSchedulerRead[rb_idx][s].ue_selected;
			if (sched_ue == -1 || sched_ue == self_idx) continue;
			Ryy_base += H_bar.col(s) * H_bar.col(s).adjoint();
		}
	} else {
		// SIC=0, 1, 3: include all streams initially
		Ryy_base = H_bar * H_bar.adjoint() + linear_interference * Identity;
	}

	// SIC=3: sort layers by initial SINR (descending) for sequential decoding
	std::vector<int> layer_order(R_act);
	for (int i = 0; i < R_act; i++) layer_order[i] = i;
	if (g_use_sic == 3 && R_act > 1) {
		std::vector<Real> sinr_init(R_act);
		MatrixXcReal Ryy_inv_init = Ryy_base.inverse();
		for (int sidx = 0; sidx < R_act; sidx++) {
			const MatrixXcReal h_u = H_bar.col(my_streams[sidx]);
			const MatrixXcReal u = h_u.adjoint() * Ryy_inv_init;
			Real sig = std::norm((u * h_u)(0,0));
			Real intf = 0.0;
			for (int s = 0; s < mx_ue_mumimo; s++) {
				int other_ue = ppSchedulerRead[rb_idx][s].ue_selected;
				if (other_ue == -1 || other_ue == self_idx) continue;
				intf += std::norm((u * H_bar.col(s))(0,0));
			}
			intf += (linear_interference + noise) * (u * u.adjoint()).norm();
			sinr_init[sidx] = (intf > 1e-30) ? sig / intf : 0.0;
		}
		std::sort(layer_order.begin(), layer_order.end(),
			[&](int a, int b) { return sinr_init[a] > sinr_init[b]; });
	}

	// Process layers (in SIC=3 order if applicable)
	MatrixXcReal Ryy_cur = Ryy_base;
	MatrixXcReal Ryy_cur_inv = Ryy_cur.inverse();
	std::vector<Real> per_layer_sinr_thisRB(R_act, 0.0);
	// Track which own UE layers have been cancelled (SIC=3 only).
	// cancelled[sidx] = true once layer sidx has been "decoded" and subtracted.
	std::vector<bool> cancelled(R_act, false);

	for (int order_idx = 0; order_idx < R_act; order_idx++)
	{
		int sidx = layer_order[order_idx];
		int stream_num = my_streams[sidx];
		const MatrixXcReal h_u = H_bar.col(stream_num);
		const MatrixXcReal u   = h_u.adjoint() * Ryy_cur_inv;

		const ComplexReal A_scalar = (u * h_u)(0, 0);
		Real sig = std::norm(A_scalar);

		Real intf = 0.0;
		for (int s = 0; s < mx_ue_mumimo; s++)
		{
			int other_ue = ppSchedulerRead[rb_idx][s].ue_selected;
			if (other_ue == -1) continue;
			if (s == stream_num) continue;            // exclude this layer itself

			if (other_ue == self_idx) {
				// Own UE's sibling layer: handling depends on SIC mode
				if (g_use_sic == 2) continue;            // perfectly cancelled (ideal)
				if (g_use_sic == 3) {
					// SIC=3: only cancelled (already processed) layers are excluded.
					// Find sibling layer index among my_streams to check cancelled[].
					int sib_idx = -1;
					for (int k = 0; k < R_act; k++) if (my_streams[k] == s) { sib_idx = k; break; }
					if (sib_idx >= 0 && cancelled[sib_idx]) continue;
					// else: still uncancelled, include as interference
				}
				// SIC=1 (legacy): exclude own UE entirely (ideal-but-inconsistent)
				if (g_use_sic == 1) continue;
				// SIC=0: include as interference (fall through)
			}
			intf += std::norm((u * H_bar.col(s))(0, 0));
		}
		Real u_norm_sq = (u * u.adjoint()).norm();
		intf += (linear_interference + noise) * u_norm_sq;

		Real sinr_l = (intf > 1e-30) ? sig / intf : 0.0;
		per_layer_sinr_thisRB[sidx] = sinr_l;

		// SIC=3: mark this layer as cancelled and remove from Ryy for downstream layers.
		// (Optimistic: we assume successful decode at SINR computation. The decode-fail
		//  probability is captured via BLER lookup at the resulting per-layer SINR.)
		if (g_use_sic == 3 && order_idx < R_act - 1) {
			cancelled[sidx] = true;
			Ryy_cur -= h_u * h_u.adjoint();
			Real reg = 1e-12 * (Ryy_cur.trace().real() / (Real)NUM_RX_Port);
			Ryy_cur += reg * Identity;
			Ryy_cur_inv = Ryy_cur.inverse();
		}
	}

	// Save per-layer SINR per-RB into accumulator and global counters
	for (int sidx = 0; sidx < R_act && sidx < 4; sidx++) {
		Real s = per_layer_sinr_thisRB[sidx];
		if (rbs_rx_per_layer[sidx].size() < (size_t)num_rb)
			rbs_rx_per_layer[sidx].push_back(s);
		if (s > 0) {
			recv_perlayer_sinr_sum[sidx+1] += s;
			recv_perlayer_count[sidx+1]++;
		}
	}

	// For legacy compatibility, return average per-layer SINR (NOT sum).
	Real avg_sinr = 0.0;
	int avg_n = 0;
	for (int sidx = 0; sidx < R_act; sidx++) {
		if (per_layer_sinr_thisRB[sidx] > 0) {
			avg_sinr += per_layer_sinr_thisRB[sidx];
			avg_n++;
		}
	}
	const Real Received_SINR = (avg_n > 0) ? avg_sinr / avg_n : 0.0;

	if (Received_SINR > 0.0)
		compute_CQI_offset(Received_SINR, rb_idx);

	return Received_SINR;
}

void MS::compute_CQI_offset(Real Received_SINR, int rb_idx)
{
	num_added++;

	Real difference = linear2dB(Received_SINR) - linear2dB(ppppCQI_Map[receive_self_ms_idx][0][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length]);
	CQI_offset = (CQI_offset * (num_added - 1) + difference) / (Real)num_added;
}


/*===================================================================
FUNCTION: UE::MESC_receiver_filter( )

DESCRIPTION:
===================================================================*/
MatrixXcReal MS::MMSE_Receiver(int rb_idx, int stream_num)
{
	MatrixXcReal H_bar, w, h;

	MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);


	Real geometry = links[receive_self_ms_idx].str_signal - links[receive_self_ms_idx].interference; // �̰Ŵ� dB // �̰� (link[_self_idx].UE_PATH_LOSS / link[_self_idx].UE_ICI) �̰Ŷ� ������ (linear)

	Real linear_signal = pow(10., links[receive_self_ms_idx].str_signal / 10.);
	Real linear_interference = pow(10., links[receive_self_ms_idx].interference / 10.);

	H_bar = sqrt(linear_signal) * H_m[0][rb_idx] * sector[links[receive_self_ms_idx]._sector_in_control].W[rb_idx];
/// Rx2                              RxT                    Tx2

	h = (H_bar.col(stream_num));
//  Rx1

	w = h.adjoint() * (H_bar * H_bar.adjoint() + linear_interference * Identity).inverse();
///1xR   1xR            Rx2       2xR                                     RxR 



	return w;

}


/*===================================================================
FUNCTION: UE::compute_tone_SINR( )

DESCRIPTION:
===================================================================*/
Real MS::compute_tone_SINR_MRC(int rb_idx)
{

	if (CH_CAL == 1) //// use run-time calculate
	{
		Real Received_SINR = 0.;

		Real  sum_of_ch_coef = 0.;

		for (int r_idx = 0; r_idx < NUM_RX_Port; r_idx++) // Received    Channel
		{
			for (int t_idx = 0; t_idx < NUM_TX_Port; t_idx++) // Transmitted Channel
			{
				sum_of_ch_coef += pow(abs(H_m[0][rb_idx](r_idx, t_idx)), 2.);

			}
		}
		Real geometry = links[receive_self_ms_idx].str_signal - links[receive_self_ms_idx].interference; // �̰Ŵ� dB // �̰� (links[_self_idx].UE_PATH_LOSS / links[_self_idx].UE_ICI) �̰Ŷ� ������ (linear)
		Real linear_signal = dBm2linear(links[receive_self_ms_idx].str_signal);
		Real linear_interference = dBm2linear(links[receive_self_ms_idx].interference);

		linear_interference = linear_interference + noise;
		Received_SINR = sum_of_ch_coef * (linear_signal / linear_interference);   // 
		return Received_SINR;
	}
	else
	{

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// For Load Pre-calculate Parameter ///////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int los_idx;
		int indoor_idx;
		int v_idx;

		int SELF_BS_idx = (int)(links[receive_self_ms_idx]._sector_in_control / 3);

		if (channel[SELF_BS_idx][receive_self_ms_idx].LOS == 1)
		{
			los_idx = 1;
		}
		else
		{
			los_idx = 0;
		}

		if (Indoor == 1) //////// InH������ ���� ����
		{
			indoor_idx = 1;
		}
		else
		{
			indoor_idx = 0;
		}

		Real distance = Get_distance(ms[receive_self_ms_idx].loc, bs[SELF_BS_idx].loc);

		if (TYPE == 1 || TYPE == 2 || TYPE == 3 || TYPE == 4 || TYPE == 5)
		{
			Real angle_v = atan2((bs_height - ms_height), distance) * 180. / pi;

			//// 
			if (angle_v > MAX_v_angle - (MAX_v_angle / anglev_mapping_value))
			{
				v_idx = 2;
			}
			else if (angle_v < (MAX_v_angle / anglev_mapping_value))
			{
				v_idx = 0;
			}
			else
			{
				v_idx = 1;
			}
		}
		else
		{
			v_idx = 0;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Real Received_SINR = 0.;
		Real  sum_of_ch_coef = 0.;

		sum_of_ch_coef = sum_of_ch_coef_precalc[t % precalculate_time][los_idx][indoor_idx][v_idx][rb_idx];

		Real geometry = links[receive_self_ms_idx].str_signal - links[receive_self_ms_idx].interference; // �̰Ŵ� dB // �̰� (links[_self_idx].UE_PATH_LOSS / links[_self_idx].UE_ICI) �̰Ŷ� ������ (linear)
		Real linear_signal = dBm2linear(links[receive_self_ms_idx].str_signal);
		Real linear_interference = dBm2linear(links[receive_self_ms_idx].interference);

		linear_interference = linear_interference + noise;
		Received_SINR = sum_of_ch_coef * (linear_signal / linear_interference);

		return Received_SINR;
	}
}


