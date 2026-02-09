#include "common.h"
#include "const.h"
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

Real MS::BLER_Mapping( void )
{
	int mapping_snr_dB_idx = 0;
	int mapping_cqi_idx    = 0;
	Real BLER_value      = 0.;

	Real x_axis = 10 * log10(ESINR_linear);		
	for (int snr_idx = 0; snr_idx < 1426; snr_idx++)
	{
		if (SNR_5G_dB[snr_idx] > x_axis)
		{
			mapping_snr_dB_idx = snr_idx - 1;           
			if (mapping_snr_dB_idx < 0)
			{
				mapping_snr_dB_idx = 0;
			}
			break;
		}
		else  
		{
			mapping_snr_dB_idx = 1425;
		}
	}
	/////////////////////////////////////////////////////////////////////////
	BLER_value     = (NEW5GBLER[mapping_snr_dB_idx][_cqi_idx - 1]);
	ue_BLER_Value  = BLER_value;

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
		////////// Vienna BLER mapping //////////////////
		int mapping_snr_dB_idx = 0;
		int mapping_cqi_idx    = 0;
		Real BLER_value      = 0.;

		Real x_axis = 10 * log10(ESINR_linear);		
		for (int snr_idx = 0; snr_idx < 1426; snr_idx++)
		{
			if (SNR_5G_dB[snr_idx] > x_axis)
			{
				mapping_snr_dB_idx = snr_idx - 1;           

				if (mapping_snr_dB_idx < 0)
				{
					mapping_snr_dB_idx = 0;
				}

				break;
			}
			else  
			{
				mapping_snr_dB_idx = 1425;
			}
		}
		/////////////////////////////////////////////////////////////////////////
		BLER_value     = (NEW5GBLER[mapping_snr_dB_idx][_cqi_idx - 1]);
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

		if (_mcs_idx <= 4)       { _mod_type_ = QPSK; }
		else if (_mcs_idx <= 10) { _mod_type_ = QAM16; }
		else if (_mcs_idx <= 19) { _mod_type_ = QAM64; }
		else if (_mcs_idx <= 28) { _mod_type_ = QAM256; }
		else                      { cerr << "Something Wrong - rate matching" << endl; }


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


		N_info = (int)(N_RE * coding_rate * _mod_type * NUM_UE_Layer* (1. - overhead));

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

		if(Scheduling_Type == 2 && comp_mode[self_idx] == 1)
		{
			//rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].adj_sector[1], freq_idx);
			rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].comp_sector_idx, freq_idx);
			if (rx_sinr)
			{
				rbs_rx.push_back(rx_sinr);
				if (links[self_idx].ACK == 1)
					comp_mode_rx_flag =  1;
			}
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
		N_RE = MIN(156, N_RE_a) * num_rx_rbs; //rb_indices.size();

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
		else if (_mcs_idx == 16)  { coding_rate_x1024 = 719; }
		else if (_mcs_idx == 17)  { coding_rate_x1024 = 772; }
		else if (_mcs_idx == 18)  { coding_rate_x1024 = 822; }
		else if (_mcs_idx == 19)  { coding_rate_x1024 = 873; }
		else if (_mcs_idx == 20)  { coding_rate_x1024 = 682.5;}
		else if (_mcs_idx == 21)  { coding_rate_x1024 = 711; }
		else if (_mcs_idx == 22)  { coding_rate_x1024 = 754; }
		else if (_mcs_idx == 23)  { coding_rate_x1024 = 797; }
		else if (_mcs_idx == 24)  { coding_rate_x1024 = 841; }
		else if (_mcs_idx == 25)  { coding_rate_x1024 = 885; }
		else if (_mcs_idx == 26)  { coding_rate_x1024 = 916.5; }
		else if (_mcs_idx == 27)  { coding_rate_x1024 = 948; }
		else
		{
			cout << "Something Wrong rate matching" << endl;
			getchar();
		}

		coding_rate = coding_rate_x1024 / 1024.;
		//N_info = (int)(N_RE * coding_rate * _mod_type * NUM_UE_Layer* (1. - overhead));
		N_info = (int)(N_RE * coding_rate * _mod_type * (1. - overhead));
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

			if(Scheduling_Type == 2 && comp_mode[self_idx] == 1)
			{
				//rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].adj_sector[1], freq_idx);
				rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].comp_sector_idx, freq_idx);
				if (rx_sinr)
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

			if(Scheduling_Type == 2 && comp_mode[self_idx] == 1)
			{
				//rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].adj_sector[1], freq_idx);
				rx_sinr = Compute_Tone_SINR_NCJT( links[self_idx].comp_sector_idx, freq_idx);
				if (rx_sinr)
				{
					rbs_rx.push_back(rx_sinr);
					comp_mode_rx_flag =  1;
				}
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
		N_info = (int)(N_RE * coding_rate * _mod_type * NUM_UE_Layer* (1. - overhead));
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
	int comp_bs_idx = (int)comp_sector_idx / 3.;

		// jhnoh 220901
	if ( comp_mode[self_idx] == 1 && 
		 g_comp_mode && 
		 sector[comp_sector_idx].ppschedulewrite[0][0].ue_selected == self_idx )
	{
		Real linear_signal       = dBm2linear(links[self_idx].str_signal); // dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
		Real linear_interference = dBm2linear(links[self_idx].interference);
		linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);

		H_bar1 = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length];

		Real comp_sector_static_gain = dBm2linear(links[self_idx].static_gain[1].first); // dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
		H_bar2 = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length];
		H_bar = H_bar1 + H_bar2;

		u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + (linear_interference + noise) * Identity).inverse();
		A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
		B = 0;
		C = (linear_interference + noise) * (u * u.adjoint()).norm();
	}
	else
	{
		Real linear_signal       = dBm2linear(links[self_idx].str_signal); // dBm2linear(bs_maxpower - channel[bs_idx][self_idx].pathloss_final);
		Real linear_interference = dBm2linear(links[self_idx].interference);
		linear_interference       -= dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]);

		H_bar = sqrt(linear_signal) * (H_m[0][rb_idx]) * PMI_vector[0][rb_idx][t % cqi_history_length];

		int comp_sector_idx = links[self_idx].adj_sector[1];
		int comp_bs_idx = (int)comp_sector_idx / 3.;

		Real comp_sector_static_gain = dBm2linear(links[self_idx].intf_w_rnd_RSRP[comp_sector_idx]); //dBm2linear(links[self_idx].static_gain[1].first); // dBm2linear(bs_maxpower - channel[comp_bs_idx][self_idx].pathloss_final);
		H_inf = sqrt(comp_sector_static_gain) * H_m[1][rb_idx] * PMI_vector[1][rb_idx][t % cqi_history_length];

		u = H_bar.adjoint() * (H_bar * H_bar.adjoint() + H_inf * H_inf.adjoint() + (linear_interference + noise) * Identity).inverse();
		A = (u * H_bar * H_bar.adjoint() * u.adjoint()).norm();
		B = (u * H_inf * H_inf.adjoint() * u.adjoint()).norm();
		C = (linear_interference + noise) * (u * u.adjoint()).norm();			
	}
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
			if ( comp_mode[self_idx] && t > 1000 )
				test111 = 12;

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
	// Find stream number for this UE (early exit optimization)
	int stream_num = -1;
	for (int i = 0; i < mx_ue_mumimo; i++)
	{
		if (ppSchedulerRead[rb_idx][i].ue_selected == self_idx)
		{
			stream_num = i;
			break;  // Early exit once found
		}
	}

	// If UE not scheduled on this RB, return 0
	if (stream_num == -1)
		return 0.0;

	// Count received RBs
	num_rx_rb++;

	// Calculate power levels (reuse across function)
	const Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
	const Real linear_signal = dBm2linear(bs_maxpower - channel[(int)(links[self_idx]._sector_in_control/3)][self_idx].pathloss_final);
	const Real linear_interference = dBm2linear(links[self_idx].interference);

	// Compute effective channel matrix H_bar = sqrt(P) * H * W
	// Dimensions: [NUM_RX_Port x mx_ue_mumimo]
	const MatrixXcReal& H = H_m[0][rb_idx];
	const MatrixXcReal& W = sector[links[self_idx]._sector_in_control].W[rb_idx];
	const MatrixXcReal H_bar = sqrt(linear_signal) * H * W;

	// Extract desired user's channel vector h_u
	// Dimensions: [NUM_RX_Port x 1]
	const MatrixXcReal h_u = H_bar.col(stream_num);

	// Compute MMSE receiver filter (inline to avoid redundant H_bar calculation)
	// w = h_u^H * (H_bar * H_bar^H + σ²I)^(-1)
	// Dimensions: [1 x NUM_RX_Port]
	const MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
	const MatrixXcReal Ryy = H_bar * H_bar.adjoint() + linear_interference * Identity;
	const MatrixXcReal u = h_u.adjoint() * Ryy.inverse();

	// Desired signal after receiver combining
	// A = u * h_u (scalar)
	const ComplexReal A_scalar = (u * h_u)(0, 0);
	const Real signal_power = std::norm(A_scalar);  // |A|²

	// Compute MU-MIMO interference power from other streams
	Real interference_power = 0.0;
	for (int stream = 0; stream < mx_ue_mumimo; stream++)
	{
		// Skip if it's the desired stream or no UE scheduled
		if (stream == stream_num || ppSchedulerRead[rb_idx][stream].ue_selected == -1)
			continue;

		// Interference from stream i: B = u * h_i
		const MatrixXcReal h_i = H_bar.col(stream);
		const ComplexReal B_scalar = (u * h_i)(0, 0);
		interference_power += std::norm(B_scalar);  // |B|²
	}

	// Add external interference and noise (with MMSE enhancement)
	// Total noise power = (I_ext + N) * ||u||²
	const Real u_norm_sq = (u * u.adjoint()).norm();
	interference_power += (linear_interference + noise) * u_norm_sq;

	// Compute SINR
	const Real Received_SINR = signal_power / interference_power;

	// Update CQI offset if SINR is valid
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


