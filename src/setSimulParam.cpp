#include <stdio.h>
#include <fstream>
#include <iostream>
//#include <direct.h>
//#include <unistd.h>
#include <sys/stat.h> 
#include <sys/types.h> 

//#ifndef WIN
//int mkdir(const char *pathname, mode_t mode);
//#endif

#include <time.h>
#include "const.h"
#include "common.h"
//#include "rv_generate.h"





using namespace std;

bool Open_configure_file(fstream &file, char *name);
double Get_parameter(fstream &file, const char *name, double _default);
char* Get_character(fstream &file, const char *name, const char* _default);
char* Get_character_time(fstream &file, const char *name, const char* _default);
void eatline(fstream &file);
void Decide_test_environment(char * test_environment_char);
void Decide_pathloss_model(char * pathloss_model_char);
void make_directory(void);


void Set_simul_param(int argc, char *argv[])
{
	cout << "= DownLink Simulation= " << endl;
	cout << endl;
	cout << "Simulation_parameters" << endl;
	cout << endl;
	if (argc == 2)
	{
		fstream infile;
		if (!Open_configure_file(infile, argv[1]))
		{
			cout << "error : Can't find the file " << "'" << argv[1] << "'" << endl;
			exit(1);
		}

		cfg_name = argv[1];

		_seed                    = int(Get_parameter(infile, "_seed", 100000));
		num_drop                 = int(Get_parameter(infile, "num_drops", 10000));
		run_times                = int(Get_parameter(infile, "run_times", 1000));
		num_MS_persector         = int(Get_parameter(infile, "num_user_per_sector", 10));
		num_Indoor_TRxP          = 3;   // hardcoded: 3 sectors per BS
		simple_num_BS            = int(Get_parameter(infile, "simple_num_BS", 19));
		num_rb                   = int(Get_parameter(infile, "num_rb", 50));
		fft_size                 = int(Get_parameter(infile, "fft_size", 1024));
		subcarrier_spacing       = double(Get_parameter(infile, "subcarrier_spacing", 15000.));
		bandwidth                = double(Get_parameter(infile, "bandwidth", 10000000.));
		numerology               = int(Get_parameter(infile, "numerology", 0));
		carrier_freq             = double(Get_parameter(infile, "carrier_freq", 4000000000.));
		Scheduling_Type          = int(Get_parameter(infile, "Scheduling_Type", 1));         //// 1 = MU MIMO scheduling, 0 = roundrobin scheduling
		Calibration_mode         = int(Get_parameter(infile, "Calibration_mode", 0));        //// 1= Calibration mode, 0 = Full simulation
		TDD_mode                 = int(Get_parameter(infile, "TDD_mode", 0));                //// 1 = TDD reciprocity-based precoding, 0 = FDD codebook-based precoding
		USE_POWER_ITERATION      = int(Get_parameter(infile, "USE_POWER_ITERATION", 0));     //// 1 = Power Iteration (fast), 0 = SelfAdjointEigenSolver (accurate)
		POWER_ITERATION_MAX_ITER = int(Get_parameter(infile, "POWER_ITERATION_MAX_ITER", 5)); //// Number of iterations for Power Iteration (3-10 recommended)
		USE_RAY_LEVEL_DOPPLER    = int(Get_parameter(infile, "USE_RAY_LEVEL_DOPPLER", 0));   //// 1 = Ray-level Doppler (accurate), 0 = Cluster-average Doppler (fast)
		USE_PRECODING_BASED_SINR = int(Get_parameter(infile, "USE_PRECODING_BASED_SINR", 0)); //// 1 = Calculate SINR from precoding matrix, 0 = Use feedback CQI (default)
		INTERCELL_INTERFERENCE_MARGIN_DB = Real(Get_parameter(infile, "INTERCELL_INTERFERENCE_MARGIN_DB", 0.0)); //// Inter-cell interference margin in dB (for future use)
		single_cell_mode         = int(Get_parameter(infile, "single_cell_mode", 0));        //// 1= Center cell only (BS 0), 0 = All cells 
		BS_M                     = int(Get_parameter(infile, "BS_num_of_ant_elements_with_same_pol_each_col(M)", 1));
		BS_N                     = int(Get_parameter(infile, "BS_num_of_ant_elements_columns(N)", 1));				
		BS_P                     = int(Get_parameter(infile, "BS_antenna_polarization(P)", 1));
		BS_Mg                    = int(Get_parameter(infile, "BS_num_of_panels_in_a_column(Mg)", 1));
		BS_Ng                    = int(Get_parameter(infile, "BS_num_of_panels_in_a_row(Ng)", 1));
		BS_Mp                    = int(Get_parameter(infile, "BS_Mp", 1));
		BS_Np                    = int(Get_parameter(infile, "BS_Np", 1));
		BS_dH                    = double(Get_parameter(infile, "BS_dH", 0.5));
		BS_dV                    = double(Get_parameter(infile, "BS_dV", 0.5));
		BS_dgH                   = double(Get_parameter(infile, "BS_dgH", 0.5));
		BS_dgV                   = double(Get_parameter(infile, "BS_dgV", 0.5));
		MS_M                     = int(Get_parameter(infile, "MS_num_of_ant_elements_with_same_pol_each_col(M)", 1));
		MS_N                     = int(Get_parameter(infile, "MS_num_of_ant_elements_columns(N)", 1));
		MS_P                     = int(Get_parameter(infile, "MS_antenna_polarization(P)", 1));
		MS_Mg                    = int(Get_parameter(infile, "MS_num_of_panels_in_a_column(Mg)", 1));
		MS_Ng                    = int(Get_parameter(infile, "MS_num_of_panels_in_a_row(Ng)", 1));
		MS_Mp                    = int(Get_parameter(infile, "MS_Mp", 1));
		MS_Np                    = int(Get_parameter(infile, "MS_Np", 1));
		MS_dH                    = double(Get_parameter(infile, "MS_dH", 0.5));
		MS_dV                    = double(Get_parameter(infile, "MS_dV", 0.5));
		MS_dgH                   = double(Get_parameter(infile, "MS_dgH", 0.5));
		MS_dgV                   = double(Get_parameter(infile, "MS_dgV", 0.5));
		Mechanic_tilt            = double(Get_parameter(infile, "Mechanic_tilt", 90));
		num_compute_coef         = int(Get_parameter(infile, "num_compute_coef", 1));
		num_neighbor             = int(Get_parameter(infile, "num_neighbor", 4));
		grid_interval            = int(Get_parameter(infile, "grid_interval", 10));
		cqi_history_length       = int(Get_parameter(infile, "cqi_history_length", 5));
		mx_ue_mumimo             = int(Get_parameter(infile, "mx_ue_mumimo", 2));
		NUM_UE_Layer             = int(Get_parameter(infile, "NUM_UE_Layer", 2));
		N_pf                     = double(Get_parameter(infile, "N_pf", 50.));
		SCHEDULE_DELAY           = int(Get_parameter(infile, "SCHEDULE_DELAY", 6));
		num_of_threads           = int(Get_parameter(infile, "num_threads", 32));
		g_use_min_sinr_for_mcs   = int(Get_parameter(infile, "use_min_sinr_for_mcs", 0));

		// OLLA (Outer Loop Link Adaptation) parameters
		g_olla_enable            = int(Get_parameter(infile, "olla_enable", 0));
		g_olla_step_down         = Get_parameter(infile, "olla_step_down", 0.1);
		g_olla_step_up           = Get_parameter(infile, "olla_step_up", 0.1);
		g_olla_min_offset        = Get_parameter(infile, "olla_min_offset", -10.0);
		g_olla_max_offset        = Get_parameter(infile, "olla_max_offset", 3.0);
		g_olla_init_offset       = Get_parameter(infile, "olla_init_offset", 0.0);
		g_olla_target_bler       = Get_parameter(infile, "olla_target_bler", 0.1);
		g_olla_bler_margin       = Get_parameter(infile, "olla_bler_margin", 0.02);
		g_olla_window_size       = int(Get_parameter(infile, "olla_window_size", 50));

		// MU-MIMO Scheduling Algorithm Selection
		g_mumimo_scheduling_algorithm = int(Get_parameter(infile, "mumimo_scheduling_algorithm", 0));
		g_chordal_alpha               = Get_parameter(infile, "chordal_alpha", 0.22);

		// Type II Codebook (TS 38.214 §5.2.2.2.3 / §5.2.2.2.5)
		g_codebook_type           = int(Get_parameter(infile, "Codebook_Type", 1));            // 1=Type I, 2=Type II, 3=eType II
		g_type2_L                 = int(Get_parameter(infile, "Type2_L", 4));                  // 2/3/4 (overridden by ParamComb if eType II)
		g_type2_phase_alphabet    = int(Get_parameter(infile, "Type2_PhaseAlphabet", 8));      // 4, 8, or 16 (for eType II)
		g_type2_subband_amplitude = int(Get_parameter(infile, "Type2_SubbandAmplitude", 1));
		g_type2_rank              = int(Get_parameter(infile, "Type2_Rank", 1));               // MAX rank cap 1..4
		g_etype2_param_comb       = int(Get_parameter(infile, "eType2_ParamComb", 3));          // 1..8
		g_rank_adaptive           = int(Get_parameter(infile, "Rank_Adaptive", 0));             // 0 = fixed, 1 = per-UE RI
		g_su_fallback             = int(Get_parameter(infile, "SU_Fallback", 0));               // 0 = always MU, 1 = SU vs MU
		g_use_sic                 = int(Get_parameter(infile, "Use_SIC", 1));                    // 0 = per-stream MMSE, 1 = ideal SIC (historical default)
		g_per_layer_mcs           = int(Get_parameter(infile, "per_layer_mcs", 0));              // 0 = legacy single-MCS/TB-HARQ, 1 = per-layer MCS + HARQ
		g_harq_ir                 = int(Get_parameter(infile, "harq_ir", 0));                    // 0 = HARQ Type I, 1 = IR/Chase SINR combining
		g_matlab_bler             = int(Get_parameter(infile, "matlab_bler", 0));                // 0 = legacy BLER_5G.dat, 1 = MATLAB L2SM (code rate + TBS aware)
		g_matlab_cqi_thresholds   = int(Get_parameter(infile, "matlab_cqi_thresholds", 0));      // 1 = regenerate SINR_threshold_dB from MATLAB tables
		g_matlab_bler_selftest    = int(Get_parameter(infile, "matlab_bler_selftest", 0));       // 1 = validate lookup vs oracle CSV at startup
		g_matlab_rbir             = int(Get_parameter(infile, "matlab_rbir", 1));                // 1 = RBIR eff-SINR (curves' axis) when matlab_bler=1; 0 = EESM (diagnostic)
		g_matlab_tput_mcs         = int(Get_parameter(infile, "matlab_tput_mcs", 0));            // 1 = throughput-max MCS over all 28 MCS (needs matlab_bler=1)
		g_tput_mcs_ref_rbs        = int(Get_parameter(infile, "tput_mcs_ref_rbs", 0));           // reference per-UE RBs for the MCS grid; 0 = auto num_rb/3
		g_matlab_esinr_fb         = int(Get_parameter(infile, "matlab_esinr_fb", 0));            // 1 = realized-ESINR feedback correction (needs matlab_bler=1, per_layer_mcs=1)

		if (g_type2_rank < 1) g_type2_rank = 1;
		if (g_type2_rank > 4) g_type2_rank = 4;

		// ParameterCombination (TR 38.214 Table 5.2.2.2.5-1 + simplification): (L, p_v, β)
		// Combinations 1-8, with different L values. L overrides Type2_L when eType II is on.
		if (g_codebook_type == 3) {
			switch (g_etype2_param_comb) {
				case 1: g_type2_L = 2; g_etype2_pv = 0.25; g_etype2_beta = 0.25; break;
				case 2: g_type2_L = 2; g_etype2_pv = 0.25; g_etype2_beta = 0.50; break;
				case 3: g_type2_L = 4; g_etype2_pv = 0.25; g_etype2_beta = 0.25; break;
				case 4: g_type2_L = 4; g_etype2_pv = 0.25; g_etype2_beta = 0.50; break;
				case 5: g_type2_L = 4; g_etype2_pv = 0.25; g_etype2_beta = 0.75; break;
				case 6: g_type2_L = 4; g_etype2_pv = 0.50; g_etype2_beta = 0.50; break;
				case 7: g_type2_L = 6; g_etype2_pv = 0.25; g_etype2_beta = 0.50; break;
				case 8: g_type2_L = 6; g_etype2_pv = 0.25; g_etype2_beta = 0.75; break;
				default: g_type2_L = 4; g_etype2_pv = 0.25; g_etype2_beta = 0.25; break;
			}
		}

		if (g_type2_L < 2) g_type2_L = 2;
		if (g_type2_L > 6) g_type2_L = 6;
		if (g_type2_phase_alphabet != 4 && g_type2_phase_alphabet != 8 && g_type2_phase_alphabet != 16) {
			cout << "[WARNING] Type2_PhaseAlphabet must be 4, 8, or 16, got " << g_type2_phase_alphabet << ", using 16" << endl;
			g_type2_phase_alphabet = 16;
		}

		// Row/Column analog beam allocation (per-drop population vote)
		row_beam_enable         = int(Get_parameter(infile, "row_beam_enable", 0));
		row_beam_az_mode        = int(Get_parameter(infile, "row_beam_az_mode", 0));
		row_beam_force_uniform  = int(Get_parameter(infile, "row_beam_force_uniform", 0));
		row_beam_max_cand       = int(Get_parameter(infile, "row_beam_max_cand", 3));
		row_beam_x_db           = Real(Get_parameter(infile, "row_beam_x_db", 3.0));
		row_beam_zenith_min_deg = Real(Get_parameter(infile, "row_beam_zenith_min_deg", 0.0));
		row_beam_zenith_max_deg = Real(Get_parameter(infile, "row_beam_zenith_max_deg", 0.0));
		row_beam_num_zenith     = int(Get_parameter(infile, "row_beam_num_zenith", 0));

		// Singular Value CDF Collection
		g_collect_singular_values = int(Get_parameter(infile, "collect_singular_values", 0));

		// Cluster parameter CSV dump
		g_dump_cluster_params = int(Get_parameter(infile, "dump_cluster_params", 0));

		// Channel parameter version: 0=V19 (default), 1=old TR 38.901 (pre-V19)
		channel_param_legacy     = int(Get_parameter(infile, "channel_param_legacy", 0));

		// Handheld UT antenna model
		handheld_mode            = int(Get_parameter(infile, "handheld_mode", 0));
		handheld_beta_deg        = Real(Get_parameter(infile, "handheld_beta", 45.0));
		handheld_num_ports       = int(Get_parameter(infile, "handheld_num_ports", 4));
		handheld_port_indices[0] = int(Get_parameter(infile, "handheld_port_1", 1));
		handheld_port_indices[1] = int(Get_parameter(infile, "handheld_port_2", 7));
		handheld_port_indices[2] = int(Get_parameter(infile, "handheld_port_3", 3));
		handheld_port_indices[3] = int(Get_parameter(infile, "handheld_port_4", 5));
		// Ports 5-8 (for 8-port handheld): default to the 4 edge-center positions (2,4,6,8),
		// complementing the 4 corner positions (1,3,5,7) used by the 4-port config.
		handheld_port_indices[4] = int(Get_parameter(infile, "handheld_port_5", 2));
		handheld_port_indices[5] = int(Get_parameter(infile, "handheld_port_6", 4));
		handheld_port_indices[6] = int(Get_parameter(infile, "handheld_port_7", 6));
		handheld_port_indices[7] = int(Get_parameter(infile, "handheld_port_8", 8));

		// BS-side Spatial Non-Stationarity (SNS)
		// Default values: UMi (Table 7.6.14.1.2-1/2/3)
		g_sns_bs_enabled   = int(Get_parameter(infile, "sns_bs_enabled", 0));
		g_sns_mu_P_vis     = Get_parameter(infile, "sns_mu_P_vis", 0.49);
		g_sns_sigma_P_vis  = Get_parameter(infile, "sns_sigma_P_vis", 0.18);
		g_sns_vr_A         = Get_parameter(infile, "sns_vr_A", 0.12);
		g_sns_vr_B         = Get_parameter(infile, "sns_vr_B", 0.48);
		g_sns_vr_R         = Get_parameter(infile, "sns_vr_R", 50.0);
		g_sns_vr_delta     = Get_parameter(infile, "sns_vr_delta", 0.0316);
		g_sns_rolloff_C    = Get_parameter(infile, "sns_rolloff_C", 13.0);

		cfg_BS_Tx_Power      = Real(Get_parameter(infile, "BS_Tx_Power", -9999));
		cfg_UT_Noise_Figure  = Real(Get_parameter(infile, "UT_Noise_Figure", -9999));
		cfg_UE_antenna_element_gain = Real(Get_parameter(infile, "ue_antenna_element_gain", -9999));
		cfg_inter_site_distance = Real(Get_parameter(infile, "inter_site_distance", -1));

		num_mTRP = simple_num_BS*3;

		// ====================================================================
		// Validate antenna configuration parameters against array size limits
		// ====================================================================
		bool config_valid = true;
		if (BS_M > MAX_BS_M) {
			cout << "ERROR: BS_M (" << BS_M << ") exceeds maximum allowed (" << MAX_BS_M << ")" << endl;
			config_valid = false;
		}
		if (BS_N > MAX_BS_N) {
			cout << "ERROR: BS_N (" << BS_N << ") exceeds maximum allowed (" << MAX_BS_N << ")" << endl;
			config_valid = false;
		}
		if (BS_P > MAX_BS_P) {
			cout << "ERROR: BS_P (" << BS_P << ") exceeds maximum allowed (" << MAX_BS_P << ")" << endl;
			config_valid = false;
		}
		if (BS_Mg > MAX_BS_Mg) {
			cout << "ERROR: BS_Mg (" << BS_Mg << ") exceeds maximum allowed (" << MAX_BS_Mg << ")" << endl;
			config_valid = false;
		}
		if (BS_Ng > MAX_BS_Ng) {
			cout << "ERROR: BS_Ng (" << BS_Ng << ") exceeds maximum allowed (" << MAX_BS_Ng << ")" << endl;
			config_valid = false;
		}
		if (MS_M > MAX_MS_M) {
			cout << "ERROR: MS_M (" << MS_M << ") exceeds maximum allowed (" << MAX_MS_M << ")" << endl;
			config_valid = false;
		}
		if (MS_N > MAX_MS_N) {
			cout << "ERROR: MS_N (" << MS_N << ") exceeds maximum allowed (" << MAX_MS_N << ")" << endl;
			config_valid = false;
		}
		if (MS_P > MAX_MS_P) {
			cout << "ERROR: MS_P (" << MS_P << ") exceeds maximum allowed (" << MAX_MS_P << ")" << endl;
			config_valid = false;
		}
		if (MS_Mg > MAX_MS_Mg) {
			cout << "ERROR: MS_Mg (" << MS_Mg << ") exceeds maximum allowed (" << MAX_MS_Mg << ")" << endl;
			config_valid = false;
		}
		if (MS_Ng > MAX_MS_Ng) {
			cout << "ERROR: MS_Ng (" << MS_Ng << ") exceeds maximum allowed (" << MAX_MS_Ng << ")" << endl;
			config_valid = false;
		}
		// Non-handheld DIRECTIONAL UEs: the beam-search subarray buffers are 8x8
		// (ue_pos local in find_best_tx_beam; ue_w/ue_v/ue_virtualization_weight_wv
		// globals), indexed by K_ue = MS_M/MS_Mp and L_ue = MS_N/MS_Np. MAX_MS_M(16)
		// admits configs that would overflow them — reject here when the cfg makes the
		// UE explicitly directional. (ue_antenna_element_gain may also come from the
		// scenario default, resolved AFTER parsing — a runtime guard at the buffer
		// fill in find_best_tx_beam covers that path.)
		if (handheld_mode == 0 &&
		    cfg_UE_antenna_element_gain != -9999 && cfg_UE_antenna_element_gain != 0.0) {
			int K_ue = (MS_Mp > 0) ? MS_M / MS_Mp : MS_M;
			int L_ue = (MS_Np > 0) ? MS_N / MS_Np : MS_N;
			if (K_ue > 8 || L_ue > 8) {
				cout << "ERROR: directional UE subarray " << K_ue << "x" << L_ue
				     << " exceeds the 8x8 beam-search buffers (MS_M/MS_Mp, MS_N/MS_Np)" << endl;
				config_valid = false;
			}
		}
		// Handheld arrays (positions, port indices, per-port ray buffers) are 8-deep.
		if (handheld_mode == 1 && (handheld_num_ports < 1 || handheld_num_ports > 8)) {
			cout << "ERROR: handheld_num_ports (" << handheld_num_ports
			     << ") must be in [1, 8]" << endl;
			config_valid = false;
		}

		if (!config_valid) {
			cout << "Configuration validation failed. Please adjust parameters or increase array limits in const.h" << endl;
			exit(1);
		}

		char * test_environment_char = Get_character(infile, "test_environment", "Dense_Urban_eMBB_A");
		char * pathloss_char = Get_character(infile, "pathloss_model", "UMa_A");

		int max_num_thread = omp_get_max_threads();
		if ( num_of_threads == 100)
			num_of_threads = max_num_thread;


		simple_wrap_mat = new int * [simple_num_BS];
		for(int bs_idx1 = 0; bs_idx1 < simple_num_BS; bs_idx1++)
		{
			simple_wrap_mat[bs_idx1] = new int [simple_num_BS];
			for(int bs_idx2 = 0; bs_idx2 < simple_num_BS; bs_idx2++)
			{
				if ( simple_num_BS == 1) {
					// 1-SITE configuration: single BS with 3 sectors
					simple_wrap_mat[bs_idx1][bs_idx2] = simple_wrap_mat_1site[bs_idx1][bs_idx2];
				}
				else if ( simple_num_BS == 7) {
					simple_wrap_mat[bs_idx1][bs_idx2] = simple_wrap_mat_1tier[bs_idx1][bs_idx2];
				}
				else if ( simple_num_BS == 19) {
					simple_wrap_mat[bs_idx1][bs_idx2] = simple_wrap_mat_2tier[bs_idx1][bs_idx2];
				}
			}
		}

		Decide_test_environment(test_environment_char);
		Decide_pathloss_model(pathloss_char);

		if (row_beam_enable) {
			// The row-beam path relies on the full-simulation flow (Assign_Row_Beams,
			// W_tx build, Get_CouplingLoss DU branch) and on the 3-sector macro beam
			// grid; calibration mode and InH are out of scope.
			if (Calibration_mode != 0) {  // any nonzero value runs calibration branches
				cout << "ERROR: row-beam allocation is not supported in Calibration_mode" << endl;
				exit(1);
			}
			if (scenario == 11) {
				cout << "ERROR: row-beam allocation is not supported for InH scenarios" << endl;
				exit(1);
			}
			if (row_beam_max_cand < 1) row_beam_max_cand = 1;
		}
		// Zenith grid override validation runs regardless of row_beam_enable: the
		// override itself is gated only on row_beam_num_zenith so that a baseline
		// A/B arm (enable 0) can use the same grid as the allocation arms.
		if (row_beam_num_zenith > 0) {
			if (Calibration_mode != 0) {  // any nonzero value runs calibration branches
				// 3GPP calibration relies on the scenario's fixed beam grids; a silent
				// override would invalidate calibration results.
				cout << "ERROR: zenith grid override is not supported in Calibration_mode" << endl;
				exit(1);
			}
			if (scenario == 11) {
				cout << "ERROR: zenith grid override is not supported for InH scenarios" << endl;
				exit(1);
			}
			if (row_beam_num_zenith > MAX_RSRP_SEC_ZENITH) {
				cout << "ERROR: row_beam_num_zenith exceeds MAX_RSRP_SEC_ZENITH ("
				     << MAX_RSRP_SEC_ZENITH << ")" << endl;
				exit(1);
			}
			if (!(row_beam_zenith_min_deg < row_beam_zenith_max_deg)) {
				cout << "ERROR: zenith grid override needs min_deg < max_deg" << endl;
				exit(1);
			}
		}
		file_name = Get_character_time(infile, "file_name", "output.txt");

		Rand new_randum((unsigned long long int)_seed);
		randnum = new_randum;
		default_random_engine new_e(_seed);
		e       = new_e;
	}

	else if (argc == 1)
	{
		_seed = 1000;
		num_drop = 100;
		run_times = 2000;
		num_MS_persector = 10;
		num_Indoor_TRxP = 3;
		simple_num_BS = 19;

		Scheduling_Type = 1;         //// 1 = MU MIMO scheduling, 0 = roundrobin scheduling
		Calibration_mode = 0;        //// 1= Calibration mode, 0 = Full simulation 

		//NUM_RX = 1;
		//NUM_TX = 1;
		scenario = 12;
		Configuration_Type = 0;
		Channel_Model_Type = 0;

		BS_M = 1;
		BS_N = 1;
		BS_P = 2;
		BS_Mg = 1;
		BS_Ng = 1;
		MS_M = 1;
		MS_N = 1;
		MS_P = 2;
		MS_Mg = 1;
		MS_Ng = 1;


		BS_Mp = 1;
		BS_Np = 1;
		MS_Mp = 1;
		MS_Np = 1;

		num_neighbor = 4;
		grid_interval = 10;

		cqi_history_length = 5;
		mx_ue_mumimo = 2;
		NUM_UE_Layer = 2;
		subcarrier_spacing = 15000.;
		N_pf = 50.;
		SCHEDULE_DELAY = 6;

		file_name = (char *)"output.txt";

	}
	else if (argc > 2)
	{
		cout << " Usage : run.exe filename " << endl;
		exit(1);
	}

	if ( test_environment < 3 ) // IMT2020 eval. 8.3.1
	{
		num_BS      = 12;
		num_SECTORS = num_BS * num_Indoor_TRxP;
		num_mTRP         = 0;
		num_mTRP_SECTORS = 0;
	}
	else if (test_environment == 6) //Dense_Urban_eMBB_C_Micro
	{
		num_BS = simple_num_BS;
		//num_BS = 19;
		num_SECTORS = num_BS * num_Indoor_TRxP;
		num_mTRP         = 0;
		num_mTRP_SECTORS = 0;
	}
	else
	{
		num_BS = simple_num_BS;
		//num_BS = 19;
		num_SECTORS = num_BS * 3;
		num_mTRP         = 0;
		num_mTRP_SECTORS = 0;
	}


	// In single_cell_mode, only create users for center BS (3 sectors)
	if (single_cell_mode == 1) {
		num_MS = 3 * num_MS_persector;  // Center BS only: 3 sectors
		cout << "*** Single Cell Mode: Creating users for center BS only (3 sectors) ***" << endl;
	}
	else {
		num_MS = num_SECTORS * num_MS_persector;  // All sectors
	}

	
	NUM_RX = MS_P * MS_M*MS_N*MS_Mg*MS_Ng;
	NUM_TX = BS_P * BS_M*BS_N*BS_Mg*BS_Ng;

	NUM_RX_Port = MS_P * MS_Mp*MS_Np*MS_Mg*MS_Ng;
	NUM_TX_Port = BS_P * BS_Mp*BS_Np*BS_Mg*BS_Ng;
	
	make_directory();


	
	cout << "num_drops               : " << num_drop << endl;
	cout << "run_times               : " << run_times << endl;
	cout << "num_user_per_sector     : " << num_MS_persector << endl;
	cout << "num_rx_antenna          : " << NUM_RX << endl;
	cout << "num_tx_antenna          : " << NUM_TX << endl;
	cout << "num_rx_antenna_port     : " << NUM_RX_Port << endl;
	cout << "num_tx_antenna_port     : " << NUM_TX_Port << endl;
	cout << "num_sectors             : " << num_SECTORS << endl;
	cout << "num_neighbor            : " << num_neighbor << endl;
	cout << "grid_interval           : " << grid_interval << endl;
	cout << "num_of_threads          : " << num_of_threads << endl;
	cout << "use_min_sinr_for_mcs    : " << g_use_min_sinr_for_mcs << " (0=avg, 1=min)" << endl;
	cout << "olla_enable             : " << g_olla_enable << " (0=off, 1=on)" << endl;
	if (g_olla_enable) {
		cout << "  ├─ window_size        : " << g_olla_window_size << " (moving window)" << endl;
		cout << "  ├─ target_bler        : " << g_olla_target_bler * 100 << " %" << endl;
		cout << "  ├─ bler_margin        : " << g_olla_bler_margin * 100 << " % (hysteresis)" << endl;
		cout << "  ├─ step_down          : " << g_olla_step_down << " dB" << endl;
		cout << "  ├─ step_up            : " << g_olla_step_up << " dB" << endl;
		cout << "  ├─ init_offset        : " << g_olla_init_offset << " dB" << endl;
		cout << "  └─ offset_range       : [" << g_olla_min_offset << ", " << g_olla_max_offset << "] dB" << endl;
	}
	cout << "Scheduling_Type         : " << Scheduling_Type << endl;
	if (Scheduling_Type == 1) {
		cout << "  ├─ algorithm          : " << g_mumimo_scheduling_algorithm << " (0=SUS, 1=Chordal)" << endl;
		if (g_mumimo_scheduling_algorithm == 1) {
			cout << "  └─ chordal_alpha      : " << g_chordal_alpha << endl;
		}
	}
	cout << "Codebook_Type           : " << g_codebook_type << " (1=Type I, 2=Type II, 3=eType II)" << endl;
	if (g_codebook_type == 2) {
		cout << "  ├─ L                  : " << g_type2_L << endl;
		cout << "  ├─ PhaseAlphabet      : " << g_type2_phase_alphabet << " (" << (g_type2_phase_alphabet==4?"QPSK":"8-PSK") << ")" << endl;
		cout << "  └─ SubbandAmplitude   : " << g_type2_subband_amplitude << endl;
	}
	if (g_codebook_type == 3) {
		const char* phase_name = (g_type2_phase_alphabet==4) ? "QPSK" :
		                         (g_type2_phase_alphabet==8) ? "8-PSK" : "16-PSK";
		cout << "  ├─ Rank_max (R)       : " << g_type2_rank << endl;
		cout << "  ├─ ParamComb          : " << g_etype2_param_comb << " (L=" << g_type2_L
		     << ", p_v=" << g_etype2_pv << ", β=" << g_etype2_beta << ")" << endl;
		cout << "  ├─ PhaseAlphabet      : " << g_type2_phase_alphabet << " (" << phase_name << ")" << endl;
		cout << "  ├─ Rank_Adaptive      : " << g_rank_adaptive << " (0=fixed, 1=per-UE RI)" << endl;
		cout << "  ├─ SU_Fallback        : " << g_su_fallback << " (0=always MU, 1=SU vs MU)" << endl;
		cout << "  ├─ Use_SIC            : " << g_use_sic << " (0=per-stream MMSE, 1=ideal SIC)" << endl;
		cout << "  ├─ per_layer_mcs      : " << g_per_layer_mcs << " (0=legacy single-MCS, 1=per-layer MCS+HARQ)" << endl;
		cout << "  ├─ harq_ir            : " << g_harq_ir << " (0=Type I no-combining, 1=IR/Chase SINR accumulation)" << endl;
		cout << "  └─ matlab_bler        : " << g_matlab_bler << " (0=legacy SINRxCQI table, 1=MATLAB L2SM CR+TBS)"
		     << "  [rbir=" << g_matlab_rbir << " tput_mcs=" << g_matlab_tput_mcs << " esinr_fb=" << g_matlab_esinr_fb
		     << " cqi_thresholds=" << g_matlab_cqi_thresholds
		     << " selftest=" << g_matlab_bler_selftest << "]" << endl;
	}
	if (channel_param_legacy)
		cout << "channel_param_legacy    : " << channel_param_legacy << " (old TR 38.901 pre-V19)" << endl;
	if (row_beam_enable) {
		cout << "row_beam_enable         : " << row_beam_enable << " (per-drop population-vote analog beams)" << endl;
		cout << "  ├─ az_mode            : " << row_beam_az_mode << " (0=per-UE az, 1=boresight, 2=per-column vote)" << endl;
		cout << "  ├─ force_uniform      : " << row_beam_force_uniform << " (1=argmax beam on all rows/cols)" << endl;
		cout << "  ├─ x_db / max_cand    : " << row_beam_x_db << " dB / " << row_beam_max_cand << endl;
		if (row_beam_num_zenith > 0)
			cout << "  └─ zenith grid        : " << row_beam_num_zenith << " beams in ["
			     << row_beam_zenith_min_deg << ", " << row_beam_zenith_max_deg << "] deg" << endl;
		else
			cout << "  └─ zenith grid        : legacy (unchanged)" << endl;
	}
	if (g_sns_bs_enabled) {
		cout << "sns_bs_enabled          : " << g_sns_bs_enabled << " (j20 Stochastic VR)" << endl;
		cout << "  ├─ Pr_sns (mu/sigma)  : " << g_sns_mu_P_vis << " / " << g_sns_sigma_P_vis << endl;
		cout << "  ├─ vr_A/B             : " << g_sns_vr_A << " / " << g_sns_vr_B << endl;
		cout << "  ├─ vr_R/sigma_delta   : " << g_sns_vr_R << " / " << g_sns_vr_delta << endl;
		cout << "  └─ rolloff_C          : " << g_sns_rolloff_C << endl;
	}
	cout << "folder_name             : " << folder_name << endl;
	if (handheld_mode) {
		cout << "handheld_mode           : " << handheld_mode << endl;
		cout << "  ├─ beta               : " << handheld_beta_deg << " deg" << endl;
		cout << "  └─ ports              : " << handheld_num_ports << " (";
		for (int i = 0; i < handheld_num_ports; i++) {
			if (i > 0) cout << ",";
			cout << handheld_port_indices[i];
		}
		cout << ")" << endl;
	}

	switch (test_environment)
	{
	case InH_eMBB_A:
		cout << "test_environment        : InH_eMBB_A" << endl;
		break;
	case InH_eMBB_B:
		cout << "test_environment        : InH_eMBB_B" << endl;
		break;
	case InH_eMBB_C:
		cout << "test_environment        : InH_eMBB_C" << endl;
		break;
	case Dense_Urban_eMBB_A:
		cout << "test_environment        : Dense_Urban_eMBB_A" << endl;
		break;
	case Dense_Urban_eMBB_B:
		cout << "test_environment        : Dense_Urban_eMBB_B" << endl;
		break;
	case Dense_Urban_eMBB_C_Macro:
		cout << "test_environment        : Dense_Urban_eMBB_C_Macro" << endl;
		break;
	case Dense_Urban_eMBB_C_Micro:
		cout << "test_environment        : Dense_Urban_eMBB_C_Micro" << endl;
		break;
	case Rural_eMBB_A:
		cout << "test_environment        : Rural_eMBB_A" << endl;
		break;
	case Rural_eMBB_B:
		cout << "test_environment        : Rural_eMBB_B" << endl;
		break;
	case Rural_eMBB_C:
		cout << "test_environment        : Rural_eMBB_C" << endl;
		break;
	default:
		cout << "check test_environment input file" << endl;
		break;
	}

	switch (pathloss_model)
	{
	case InH_A:
		cout << "pathloss_model          : InH_A" << endl;
		break;
	case InH_B:
		cout << "pathloss_model          : InH_B" << endl;
		break;
	case UMa_A:
		cout << "pathloss_model          : UMa_A" << endl;
		break;
	case UMa_B:
		cout << "pathloss_model          : UMa_B" << endl;
		break;
	case UMi_A:
		cout << "pathloss_model          : UMi_A" << endl;
		break;
	case UMi_B:
		cout << "pathloss_model          : UMi_B" << endl;
		break;
	case RMa_A:
		cout << "pathloss_model          : RMa_A" << endl;
		break;
	case RMa_B:
		cout << "pathloss_model          : RMa_B" << endl;
		break;
	case InH_ETRI:
		cout << "pathloss_model          : InH_ETRI" << endl;
		break;
	case UMi_ETRI:
		cout << "pathloss_model          : UMi_ETRI" << endl;
		break;
	default:
		cout << "check pathloss_model input file" << endl;
		break;
	}
	//getchar();
}

bool Open_configure_file(fstream &file, char *name)
{
	bool status;

	file.open(name, ios::in);
	if (file.fail())
		status = false;
	else
		status = true;
	return status;
}

double Get_parameter(fstream &file, const char *name, double _default)
{
	double argument = 0.0;
	bool I_got_it = false;
	char temp[100];
	while (!file.eof() && !I_got_it)
	{

		char comment_notice;
		//file >> comment_notice;
		// cout<<"comment_notice :"<<comment_notice<<endl;
		// if (comment_notice == '#')
		// 	eatline(file);
		//file.seekg(-1, ios::cur);
		file >> temp;

		if (!strcmp(name, temp))
		{
			I_got_it = true;
			file >> argument;
			
		}
	}
	file.clear();
	file.seekg(0, ios::beg);

	if (I_got_it)
		return(argument);
	else
		return(_default);
}

char* Get_character(fstream &file, const char *name, const char* _default)
{
	char* argument = new char[80];;
	bool I_got_it = false;
	char temp[80];

	while (!file.eof() && !I_got_it)
	{
		char comment_notice;
		// file >> comment_notice;
		// if (comment_notice == '#')
		// 	eatline(file);

		// file.seekg(-1, ios::cur);
		file >> temp;
		if (!strcmp(name, temp))
		{
			I_got_it = true;
			file >> argument;
		}
	}
	file.clear();
	file.seekg(0, ios::beg);

	if (I_got_it)
		return(argument);
	else
		return((char*)_default);
}

char* Get_character_time(fstream &file, const char *name, const char* _default)
{
	char* argument = new char[80];;
	bool I_got_it = false;
	char temp[80];

	while (!file.eof() && !I_got_it)
	{
		char comment_notice;
		file >> comment_notice;
		if (comment_notice == '#')
			eatline(file);

		file.seekg(-1, ios::cur);
		file >> temp;
		if (!strcmp(name, temp))
		{
			I_got_it = true;
			file >> argument;

			time_t t = time(NULL);
			struct tm tm = *localtime(&t);

			char c_year[10];
			char c_mon [10];
			char c_day [10];
			char c_hour[10];
			char c_min [10];

			g_year = tm.tm_year + 1900;
			g_mon = tm.tm_mon + 1;
			g_day = tm.tm_mday;
			g_hour = tm.tm_hour;
			g_min = tm.tm_min;

			sprintf(c_year, "%04d", g_year);
			sprintf(c_mon, "%02d",  g_mon);
			sprintf(c_day, "%02d",  g_day);
			sprintf(c_hour, "%02d", g_hour);
			sprintf(c_min, "%02d",  g_min);

			strcat(argument, c_year);
			strcat(argument, "-");
			strcat(argument, c_mon);
			strcat(argument, "-");
			strcat(argument, c_day);
			strcat(argument, "_");
			strcat(argument, c_hour);
			strcat(argument, "-");
			strcat(argument, c_min);
			strcat(argument, ".txt");
		}
	}
	file.clear();
	file.seekg(0, ios::beg);

	if (I_got_it)
	{
		return(argument);
	}

	else
		return((char*)_default);
}

void eatline(fstream &file)
{
	char null_input[1000];
	file.getline(null_input, 1000);
}



void Decide_test_environment(char * test_environment_char)
{
	char *test_environment_type_[10];

	test_environment_type_[0] = (char *)"InH_eMBB_A";
	test_environment_type_[1] = (char *)"InH_eMBB_B";
	test_environment_type_[2] = (char *)"InH_eMBB_C";
	test_environment_type_[3] = (char *)"Dense_Urban_eMBB_A";
	test_environment_type_[4] = (char *)"Dense_Urban_eMBB_B";
	test_environment_type_[5] = (char *)"Dense_Urban_eMBB_C_Macro";
	test_environment_type_[6] = (char *)"Dense_Urban_eMBB_C_Micro";
	test_environment_type_[7] = (char *)"Rural_eMBB_A";
	test_environment_type_[8] = (char *)"Rural_eMBB_B";
	test_environment_type_[9] = (char *)"Rural_eMBB_C";

	for (int idx = 0; idx < 10; idx++)
	{
		if (!strcmp(test_environment_char, test_environment_type_[idx]))
		{
			test_environment = idx;
		}
	}

	if (test_environment == 0) //InH_eMBB_A
	{
		scenario = 11;
		Configuration_Type = 0;
	}
	else if (test_environment == 1) //InH_eMBB_B
	{
		scenario = 11;
		Configuration_Type = 1;
	}
	else if (test_environment == 2) //InH_eMBB_C
	{
		scenario = 11;
		Configuration_Type = 2;
	}
	else if (test_environment == 3) //Dense_Urban_eMBB_A
	{
		scenario = 12;
		Configuration_Type = 0;
	}
	else if (test_environment == 4) //Dense_Urban_eMBB_B
	{
		scenario = 12;
		Configuration_Type = 1;
	}
	else if (test_environment == 5) //Dense_Urban_eMBB_C_Macro
	{
		scenario = 12;
		Configuration_Type = 2;
	}
	else if (test_environment == 6) //Dense_Urban_eMBB_C_Micro
	{
		scenario = 12;
		Configuration_Type = 3;
	}
	else if (test_environment == 7) //Rural_eMBB_A
	{
		scenario = 13;
		Configuration_Type = 0;
	}
	else if (test_environment == 8) //Rural_eMBB_B
	{
		scenario = 13;
		Configuration_Type = 1;
	}
	else if (test_environment == 9) //Rural_eMBB_C
	{
		scenario = 13;
		Configuration_Type = 2;
	}
}


void Decide_pathloss_model(char * pathloss_model_char)
{
	char *pathloss_type_[10];

	pathloss_type_[0] = (char *) "InH_A";
	pathloss_type_[1] = (char *) "InH_B";
	pathloss_type_[2] = (char *) "UMa_A";
	pathloss_type_[3] = (char *) "UMa_B";
	pathloss_type_[4] = (char *) "UMi_A";
	pathloss_type_[5] = (char *) "UMi_B";
	pathloss_type_[6] = (char *) "RMa_A";
	pathloss_type_[7] = (char *) "RMa_B";
	pathloss_type_[8] = (char *) "InH_ETRI";
	pathloss_type_[9] = (char *) "UMi_ETRI";

	for (int idx = 0; idx < 10; idx++)
	{
		if (!strcmp(pathloss_model_char, pathloss_type_[idx]))
		{
			pathloss_model = idx;
		}
	}

	if (pathloss_model == 0 || pathloss_model == 2 || pathloss_model == 4 || pathloss_model == 6)
	{
		Channel_Model_Type = 0;
	}
	else if (pathloss_model == 1 || pathloss_model == 3 || pathloss_model == 5 || pathloss_model == 7)
	{
		Channel_Model_Type = 1;
	}
	else if (pathloss_model == 8 || pathloss_model == 9)
	{
		Channel_Model_Type = 2;
	}
	else
	{
	}
}

void make_directory(void)
{
	for (int idx = 0; idx < 2 * sizeof(file_name); idx++)
	{
		if (file_name[idx] == '.')
		{
			if (file_name[idx + 1] == '/')
			{
				int k = 1;
				while (file_name[idx + 1 + k] != '/')
				{
					folder_name[k - 1] = file_name[idx + 1 + k];
					k++;
				}
				folder_name[k - 1] = '\0';
			}
		}
	}
#ifdef _WIN32
    mkdir(folder_name);
#else
	mkdir(folder_name , 0777);
#endif
	path[0] = '.';   path[1] = '/';   path[2] = '\0';
	strcat(path, folder_name);
}