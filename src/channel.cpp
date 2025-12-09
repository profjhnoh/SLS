#include "common.h"

Real normal(Real mu, Real sigma);
Real Get_LCS_theta(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);
Real Get_LCS_pi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);
Real Get_cos_psi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);
Real Get_sin_psi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);
void Generate_XPR(int ue_idx, int bs_idx, Real mu_XPR, Real sigma_XPR);

Real Get_BS_antenna_pattern(Real theta_GCS, Real pi_GCS, int bs_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2);
Real dot(LOCATION3D a, LOCATION3D b);
Real Transform_angle_minus_180_to_plus_180(Real x);
Real Transform_angle_0_to_plus_180(Real x);

Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2);

void CHANNEL::Reset2Default(void)
{
	self_bs = {0, 0};
	self_ms = {0, 0};
	bs_location = {0, 0};
	ms_location = {0, 0};

	LOS_prob = 0;
	pathloss = 0;
	pathloss_2 = 0;
	pathloss_final = 0;
	LOS = false;
	Propagation = 0;
	Indoor = false;
	RMS_delay_spread = 0;
	circular_angle_spread_AOA = 0;
	circular_angle_spread_AOD = 0;
	self_bs_idx = 0;
	self_ms_idx = 0;
	LOS_AOA_GCS = 0;
	LOS_AOD_GCS = 0;
	LOS_ZOA_GCS = 0;
	LOS_ZOD_GCS = 0;
	strongest_power = 0;
	strongest_power2 = 0;
	strongest_power_idx = 0;
	strongest_power_idx2 = 0;
	d_rx = 0;
	d_tx = 0;
	ms_idx = 0;
	random_phase_vv_LOS = 0;
	random_phase_hh_LOS = 0;
	num_ray = 0;
	num_path = 0;
	NUM_PATH_for_channelcoeff = 0;
	distance = 0;
	D = 0;
	r_tau = 0;
	K_factor = 0;
	DS = 0; // = sigma_tau in M.2135
	mu_K_factor = 0;
	sigma_K_factor = 0;
	mu_DS = 0;
	sigma_DS = 0;
	mu_ASD = 0;
	sigma_ASD = 0;
	mu_ASA = 0;
	sigma_ASA = 0;
	mu_ZSA = 0;
	sigma_ZSA = 0;
	mu_ZSD = 0;
	sigma_ZSD = 0;
	_DS = 0;
	K_dB = 0;
	K_linear = 0;
	ASA = 0;
	ASD = 0;
	ZSA = 0;
	ZSD = 0;
	P_1LOS = 0;
	mu_XPR = 0;
	sigma_XPR = 0;

	XPR_mean[0] = 0;
	XPR_mean[1] = 0;
	XPR_mean[2] = 0;

	XPR_std[0] = 0;
	XPR_std[1] = 0;
	XPR_std[2] = 0;

	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 20; j++)
			kappa[i][j] = 0;

	mu_offset_ZOD = 0;
	cluster_DS = 0;
	cluster_ASD = 0;
	cluster_ASA = 0;
	cluster_ZSA = 0;
	XPR = 0;
	cluster_shadowing = 0;
	sigma_SF = 0;
	BS_d_H = 0;
	BS_d_V = 0;
	BS_d_gH = 0;
	BS_d_gV = 0;
	MS_d_H = 0;
	MS_d_V = 0;
	MS_d_gH = 0;
	MS_d_gV = 0;
	alpha = 0;
	beta = 0;
	gamma = 0;
	correlated_randnum_SF = 0;
	correlated_randnum_DS = 0;
	correlated_randnum_ASD = 0;
	correlated_randnum_ASA = 0;
	correlated_randnum_K = 0;
	final_ant_gain = 0;

	for (int sec_idx = 0; sec_idx < 3; sec_idx++)
		for (int sec_zenith = 0; sec_zenith < 4; sec_zenith++)
			for (int sec_azimuth = 0; sec_azimuth < 8; sec_azimuth++)
				for (int ue_zenith_idx = 0; ue_zenith_idx < 2; ue_zenith_idx++)
					for (int ue_azimuth_idx = 0; ue_azimuth_idx < 4; ue_azimuth_idx++)
						for (int pannel = 0; pannel < 2; pannel++)
							signal_RSRP_gain[sec_idx][sec_zenith][sec_azimuth][ue_zenith_idx][ue_azimuth_idx][pannel] = 0;

	F_tx_v = complex<Real>(0, 0);
	F_tx_h = complex<Real>(0, 0);
	F_rx_v = complex<Real>(0, 0);
	F_rx_h = complex<Real>(0, 0);
	F_vv = complex<Real>(0, 0);
	F_vh = complex<Real>(0, 0);
	F_hv = complex<Real>(0, 0);
	F_hh = complex<Real>(0, 0);

	for (int cluster_idx = 0; cluster_idx < MAX_NUM_CLUSTERS; cluster_idx++)
	{
		delay[cluster_idx] = 0;
		delay_LOS[cluster_idx] = 0;
		power[cluster_idx] = 0;
		power_LOS[cluster_idx] = 0;
		power_NLOS[cluster_idx] = 0;
		AOA[cluster_idx] = 0;
		AOD[cluster_idx] = 0;
		ZOD[cluster_idx] = 0;
		ZOA[cluster_idx] = 0;
		NUM_RAY_per_ClusterNUM[cluster_idx] = 0;

		delay_power[cluster_idx] = 0;
		power_ray_sum[cluster_idx] = 0;

		AOA_power[cluster_idx] = 0;
		AOD_power[cluster_idx] = 0;
		power_ray_sum_circular[cluster_idx] = 0;

		theta_n_m_mu_AOA[cluster_idx] = 0;
		theta_n_m_mu_AOD[cluster_idx] = 0;
		theta_power_AOD[cluster_idx] = 0;
		theta_power_AOA[cluster_idx] = 0;

		for (int ray_idx = 0; ray_idx < MAX_NUM_RAYS; ray_idx++)
		{
			random_phase_vv[cluster_idx][ray_idx] = 0;
			random_phase_vh[cluster_idx][ray_idx] = 0;
			random_phase_hv[cluster_idx][ray_idx] = 0;
			random_phase_hh[cluster_idx][ray_idx] = 0;
			offset_angle[cluster_idx][ray_idx] = 0;
			offset_angle_rand_coupling[cluster_idx][ray_idx] = 0;
		}
	}

	// CHIR 초기화 (할당된 섹터만)
	if (CHIR_allocated)
	{
		for (int sec_idx = 0; sec_idx < 3; sec_idx++)
		{
			if (!sector_allocated[sec_idx]) continue;  // 이 섹터가 할당되지 않았으면 스킵

			for (int u = 0; u < NUM_TX_Port; u++)
			{
				for (int s = 0; s < NUM_RX_Port; s++)
				{
					CHIR_LOS[sec_idx][u][s] = complex<Real>(0, 0);
					for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
					{
						CHIR[sec_idx][u][s][i] = complex<Real>(0, 0);
					}

					// COMMENTED OUT: CHIR_vec not allocated to save memory
					// for (int i = 0; i < fft_size; i++ )
					// {
					// 	CHIR_vec[sec_idx][u][s][i] = complex<Real>(0, 0);
					// }
				}
			}
		}
	}
}

void CHANNEL::Set_channel(Real _distance, bool _Indoor, int _bs_idx, int _ms_idx, LOCATION _bs_location, LOCATION _ms_location)
{
	distance = _distance;
	Indoor = _Indoor;

	self_bs_idx = _bs_idx;
	self_ms_idx = _ms_idx;
	bs_location = _bs_location;
	ms_location = _ms_location;
	
	if ( scenario == 12 || scenario == 13)
	{
		if( _bs_idx < num_BS )
			_bs_height = macro_bs_height;
		else 
			_bs_height = micro_bs_height;
	}

	Set_LOS_Prob();
	Set_PATHLOSS();
	// Set_SHADOWFADING();
	Set_Channel_Parameters();
}

void CHANNEL::Set_SmallScaleParameter(int _bs_idx, int _ue_idx)
{
	K_dB = 0;
	K_linear = 0;
	P_1LOS = 0;

	if (channel[_bs_idx][_ue_idx].Propagation == LOS_propagation)
	{
		_DS = ms[_ue_idx].LSPs[_bs_idx](LOS_DS);
		K_dB = ms[_ue_idx].LSPs[_bs_idx](LOS_K);
		ASA = ms[_ue_idx].LSPs[_bs_idx](LOS_ASA);
		ASD = ms[_ue_idx].LSPs[_bs_idx](LOS_ASD);
		ZSA = ms[_ue_idx].LSPs[_bs_idx](LOS_ZSA);
		ZSD = ms[_ue_idx].LSPs[_bs_idx](LOS_ZSD);
		K_linear = pow(10, K_dB / 10);
		P_1LOS = K_linear / (K_linear + 1);
		mu_XPR = XPR_mean[LOS_propagation];
		sigma_XPR = XPR_std[LOS_propagation];
	}
	else if (channel[_bs_idx][_ue_idx].Propagation == NLOS_propagation)
	{
		_DS = ms[_ue_idx].LSPs[_bs_idx](NLOS_DS);
		ASA = ms[_ue_idx].LSPs[_bs_idx](NLOS_ASA);
		ASD = ms[_ue_idx].LSPs[_bs_idx](NLOS_ASD);
		ZSA = ms[_ue_idx].LSPs[_bs_idx](NLOS_ZSA);
		ZSD = ms[_ue_idx].LSPs[_bs_idx](NLOS_ZSD);
		mu_XPR = XPR_mean[NLOS_propagation];
		sigma_XPR = XPR_std[NLOS_propagation];
	}
	else if (channel[_bs_idx][_ue_idx].Propagation == OUT2IN_propagation)
	{
		_DS = ms[_ue_idx].LSPs[_bs_idx](OUT2IN_DS);
		ASA = ms[_ue_idx].LSPs[_bs_idx](OUT2IN_ASA);
		ASD = ms[_ue_idx].LSPs[_bs_idx](OUT2IN_ASD);
		ZSA = ms[_ue_idx].LSPs[_bs_idx](OUT2IN_ZSA);
		ZSD = ms[_ue_idx].LSPs[_bs_idx](OUT2IN_ZSD);
		mu_XPR = XPR_mean[OUT2IN_propagation];
		sigma_XPR = XPR_std[OUT2IN_propagation];
	}
	else
	{
		cout << "Something wrong with Propagation" << endl;
	}

	Set_DELAY();
	Set_POWER();
	Find_Strong2Clusters();
	Set_AOAAOD(_bs_idx, _ue_idx);
	Set_ZOAZOD(_bs_idx, _ue_idx);

	Generate_XPR(_ue_idx, _bs_idx, mu_XPR, sigma_XPR);

	/*
	if (CH_CAL == 0 || CH_CAL == 2)       // For Pre-Calculate
	{
	Set_AOAAOD_LOS();
	Set_ZOAZOD_LOS();

	Set_AOAAOD_NLOS();
	Set_ZOAZOD_NLOS();

	Set_AOAAOD();
	Set_ZOAZOD();
	}
	else
	{
	Set_AOAAOD();
	Set_ZOAZOD();
	}
	*/

	Set_InitialPhase();
	Set_SUBCLUSTER();

	Sampling_DelaySpread();


	// Set_W_matrix();

	/*
	Set_RMS_delay_spread(); // for calibration
	Set_circular_angle_spread();
	*/
}

void CHANNEL::Allocate_memory()
{
	// jhnoh 240109
	sampled_delay =  new int[fft_size];

	delay = new Real[MAX_NUM_CLUSTERS];
	delay_LOS = new Real[MAX_NUM_CLUSTERS];
	power = new Real[MAX_NUM_CLUSTERS];
	power_LOS = new Real[MAX_NUM_CLUSTERS];
	power_NLOS = new Real[MAX_NUM_CLUSTERS]; /// LOS, NLOS separate for pre-calc

	sync_delay = new Real[MAX_NUM_CLUSTERS];
	sync_delay_LOS = new Real[MAX_NUM_CLUSTERS];
	sync_power = new Real[MAX_NUM_CLUSTERS];
	sync_power_LOS = new Real[MAX_NUM_CLUSTERS];
	sync_power_NLOS = new Real[MAX_NUM_CLUSTERS]; /// LOS, NLOS separate for pre-calc


	AOA = new Real[MAX_NUM_CLUSTERS];
	AOD = new Real[MAX_NUM_CLUSTERS];
	ZOD = new Real[MAX_NUM_CLUSTERS];
	ZOA = new Real[MAX_NUM_CLUSTERS];
	NUM_RAY_per_ClusterNUM = new Real[MAX_NUM_CLUSTERS];

	delay_power = new Real[MAX_NUM_CLUSTERS];
	power_ray_sum = new Real[MAX_NUM_CLUSTERS];
	RMS_delay_power = new Real[MAX_NUM_CLUSTERS];

	AOA_power = new Real[MAX_NUM_CLUSTERS];
	AOD_power = new Real[MAX_NUM_CLUSTERS];
	power_ray_sum_circular = new Real[MAX_NUM_CLUSTERS];

	theta_n_m_mu_AOA = new Real[MAX_NUM_CLUSTERS];
	theta_n_m_mu_AOD = new Real[MAX_NUM_CLUSTERS];
	theta_power_AOD = new Real[MAX_NUM_CLUSTERS];
	theta_power_AOA = new Real[MAX_NUM_CLUSTERS];

	random_phase_vv = new Real *[MAX_NUM_CLUSTERS];
	random_phase_vh = new Real *[MAX_NUM_CLUSTERS];
	random_phase_hv = new Real *[MAX_NUM_CLUSTERS];
	random_phase_hh = new Real *[MAX_NUM_CLUSTERS];
	offset_angle = new Real *[MAX_NUM_CLUSTERS];
	offset_angle_rand_coupling = new Real *[MAX_NUM_CLUSTERS];

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		random_phase_vv[i] = new Real[MAX_NUM_RAYS];
		random_phase_vh[i] = new Real[MAX_NUM_RAYS];
		random_phase_hv[i] = new Real[MAX_NUM_RAYS];
		random_phase_hh[i] = new Real[MAX_NUM_RAYS];
		offset_angle[i] = new Real[MAX_NUM_RAYS];
		offset_angle_rand_coupling[i] = new Real[MAX_NUM_RAYS];
	}

	// CHIR/CHIR_vec/CHIR_LOS는 Allocate_CHIR_memory()에서 할당됨 (지연 할당)
}

void CHANNEL::Allocate_CHIR_memory(int sector_idx)
{
	// ====================================================================
	// MEMORY OPTIMIZATION NOTE:
	// CHIR_vec is currently not used in the optimized Fourier transform.
	// Memory allocation for CHIR_vec has been commented out to reduce
	// memory usage. If FFT-based approach is needed in future, uncomment
	// the corresponding lines below.
	// Memory savings: ~3 * NUM_TX_Port * NUM_RX_Port * fft_size * sizeof(complex<Real>)
	//                 Example: 3 * 8 * 4 * 4096 * 16 bytes = ~6.3 MB per channel
	// ====================================================================

	// 첫 할당 시 포인터 배열 생성
	if (!CHIR_allocated) {
		CHIR = new complex<Real> ***[3];
		// CHIR_vec = new complex<Real> ***[3];  // COMMENTED OUT: Not used in optimized version
		CHIR_LOS = new complex<Real> **[3];

		for (int i = 0; i < 3; i++) {
			CHIR[i] = NULL;
			// CHIR_vec[i] = NULL;  // COMMENTED OUT: Not used in optimized version
			CHIR_LOS[i] = NULL;
		}

		CHIR_allocated = true;
	}

	// 이미 이 섹터가 할당되어 있으면 스킵
	if (sector_allocated[sector_idx]) return;

	// 특정 섹터만 할당
	CHIR[sector_idx] = new complex<Real> **[NUM_TX_Port];
	// CHIR_vec[sector_idx] = new complex<Real> **[NUM_TX_Port];  // COMMENTED OUT: Not used in optimized version
	CHIR_LOS[sector_idx] = new complex<Real> *[NUM_TX_Port];

	for (int u = 0; u < NUM_TX_Port; u++)
	{
		CHIR[sector_idx][u] = new complex<Real> *[NUM_RX_Port];
		// CHIR_vec[sector_idx][u] = new complex<Real> *[NUM_RX_Port];  // COMMENTED OUT: Not used in optimized version
		CHIR_LOS[sector_idx][u] = new complex<Real>[NUM_RX_Port];

		for (int s = 0; s < NUM_RX_Port; s++)
		{
			CHIR[sector_idx][u][s] = new complex<Real>[MAX_NUM_CLUSTERS];
			// CHIR_vec[sector_idx][u][s] = new complex<Real>[fft_size];  // COMMENTED OUT: Not used in optimized version

			for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
			{
				CHIR[sector_idx][u][s][i] = complex<Real>(0, 0);
			}

			// COMMENTED OUT: CHIR_vec initialization not needed
			// for (int i = 0; i < fft_size; i++)
			// {
			// 	CHIR_vec[sector_idx][u][s][i] = complex<Real>(0, 0);
			// }

			CHIR_LOS[sector_idx][u][s] = complex<Real>(0, 0);
		}
	}

	sector_allocated[sector_idx] = true;

	// Allocate ray-level precision arrays if enabled
	if (USE_RAY_LEVEL_DOPPLER && !ray_data_allocated)
	{
		// Allocate ray_gain: [sector][tx][rx][cluster][ray]
		ray_gain = new complex<Real> ****[3];
		for (int sec = 0; sec < 3; sec++)
		{
			ray_gain[sec] = new complex<Real> ***[NUM_TX_Port];
			for (int tx = 0; tx < NUM_TX_Port; tx++)
			{
				ray_gain[sec][tx] = new complex<Real> **[NUM_RX_Port];
				for (int rx = 0; rx < NUM_RX_Port; rx++)
				{
					ray_gain[sec][tx][rx] = new complex<Real> *[MAX_NUM_CLUSTERS];
					for (int cluster = 0; cluster < MAX_NUM_CLUSTERS; cluster++)
					{
						ray_gain[sec][tx][rx][cluster] = new complex<Real>[MAX_NUM_RAYS];
						for (int ray = 0; ray < MAX_NUM_RAYS; ray++)
						{
							ray_gain[sec][tx][rx][cluster][ray] = complex<Real>(0, 0);
						}
					}
				}
			}
		}

		// Allocate ray angle arrays: [cluster][ray][2] - [0]=azimuth, [1]=zenith
		ray_AOA = new Real **[MAX_NUM_CLUSTERS];
		ray_AOD = new Real **[MAX_NUM_CLUSTERS];
		for (int cluster = 0; cluster < MAX_NUM_CLUSTERS; cluster++)
		{
			ray_AOA[cluster] = new Real *[MAX_NUM_RAYS];
			ray_AOD[cluster] = new Real *[MAX_NUM_RAYS];
			for (int ray = 0; ray < MAX_NUM_RAYS; ray++)
			{
				ray_AOA[cluster][ray] = new Real[2];  // [0]=azimuth, [1]=zenith
				ray_AOD[cluster][ray] = new Real[2];
				ray_AOA[cluster][ray][0] = 0.0;
				ray_AOA[cluster][ray][1] = 0.0;
				ray_AOD[cluster][ray][0] = 0.0;
				ray_AOD[cluster][ray][1] = 0.0;
			}
		}

		ray_data_allocated = true;
	}
}

void CHANNEL::Delete_memory()
{
	delete[] delay;
	delete[] delay_LOS;
	delete[] power;
	delete[] power_LOS;
	delete[] power_NLOS; /// LOS, NLOS separate for pre-calc
	delete[] AOA;
	delete[] AOD;
	delete[] ZOD;
	delete[] ZOA;
	delete[] NUM_RAY_per_ClusterNUM;

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		delete[] random_phase_vv[i];
		delete[] random_phase_vh[i];
		delete[] random_phase_hv[i];
		delete[] random_phase_hh[i];

		delete[] offset_angle[i];
	}

	delete[] random_phase_vv;
	delete[] random_phase_vh;
	delete[] random_phase_hv;
	delete[] random_phase_hh;

	delete[] offset_angle;
}

void CHANNEL::Delete_CHIR_memory()
{
	if (!CHIR_allocated) return;  // 할당되지 않았으면 스킵

	// 할당된 섹터만 해제
	for (int sec_idx = 0; sec_idx < 3; sec_idx++)
	{
		if (!sector_allocated[sec_idx]) continue;  // 이 섹터가 할당되지 않았으면 스킵

		for (int u = 0; u < NUM_TX_Port; u++)
		{
			for (int s = 0; s < NUM_RX_Port; s++)
			{
				delete[] CHIR[sec_idx][u][s];
				// delete[] CHIR_vec[sec_idx][u][s];  // COMMENTED OUT: CHIR_vec not allocated
			}
			delete[] CHIR[sec_idx][u];
			// delete[] CHIR_vec[sec_idx][u];  // COMMENTED OUT: CHIR_vec not allocated
			delete[] CHIR_LOS[sec_idx][u];
		}
		delete[] CHIR[sec_idx];
		// delete[] CHIR_vec[sec_idx];  // COMMENTED OUT: CHIR_vec not allocated
		delete[] CHIR_LOS[sec_idx];

		sector_allocated[sec_idx] = false;
	}

	delete[] CHIR;
	// delete[] CHIR_vec;  // COMMENTED OUT: CHIR_vec not allocated
	delete[] CHIR_LOS;

	CHIR = NULL;
	// CHIR_vec = NULL;  // COMMENTED OUT: CHIR_vec not allocated
	CHIR_LOS = NULL;
	CHIR_allocated = false;
}

void CHANNEL::Set_LOS_Prob()
{

	//////////////////////////////
	///
	/// LOS Prob
	///
	//////////////////////////////////

	//////////////////////////////////////////////////////// IMT-2020 /////////////////////////////////////////////////////////////////////////////////////////////////
	////// IMT 2020 EVAL Page 46

	if (TYPE == 11) // InH
	{
		if (distance <= 5)
		{
			LOS_prob = 1;
		}
		else if (distance > 5. && distance < 49.)
		{
			LOS_prob = exp(-(distance - 5.) / 70.8);
		}
		else if (distance >= 49.)
		{
			LOS_prob = exp(-(distance - 49.) / 211.7) * 0.54;
		}

		Real P = randnum.u();

		if (P < LOS_prob)
		{
			LOS = true; ///  LOS
			Propagation = LOS_propagation;
		}
		else
		{
			LOS = false; ///  NLOS
			Propagation = NLOS_propagation;
		}

		ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height;
		ms[self_ms_idx].floor_idx = 0;

		// LOS = true;   ///  LOS
		// Propagation = LOS_propagation;
		// cout << "LOS = true" << endl;
	}
	else if (TYPE == 12) // Dense_Urban
	{
		Real n;
		Real N;

		if (self_bs_idx == 0)
		{
			n = 0;
			N = floor((5. * randnum.u()) + 4.); // 4 ~ 8
			n = floor(((N)*randnum.u()) + 1.);	// 1 ~ N

			ms_height_in = 3 * (n - 1) + 1.5;
		}

		//////////////////////////////////////// UMi
		if (Configuration_Type == 3)
		{
			///////////////////////////////// outdoor
			if (Indoor == false)
			{
				if (distance <= 18.)
				{
					LOS_prob = 1;
				}
				else if (distance > 18.)
				{
					LOS_prob = (18. / distance) + exp(-1 * distance / 36.) * (1. - (18. / distance));
				}

				Real P = randnum.u();
				if (P < LOS_prob)
				{
					LOS = true;
					Propagation = LOS_propagation;
				}
				else
				{
					LOS = false;
					Propagation = NLOS_propagation;
				}

				if (self_bs_idx == 0)
				{
					ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_out;
					ms[self_ms_idx].floor_idx = 0;
				}
			}
			///////////////////////////////////////// Indoor
			else if (Indoor == true)
			{

				Real d_in;
				Real d_out;

				d_in = MIN(25. * randnum.u(), 25. * randnum.u()); //// UMa, UMi = 0~25, RMa = 0~10
				d_out = distance - d_in;

				if (distance <= 18.)
				{
					LOS_prob = 1;
				}
				else if (distance > 18.)
				{
					LOS_prob = (18. / d_out) + exp(-1 * d_out / 36.) * (1. - (18. / d_out));
				}

				Real P = randnum.u();
				if (P < LOS_prob)
				{
					LOS = true;
					Propagation = LOS_propagation;
				}
				else
				{
					LOS = false;
					Propagation = NLOS_propagation;
				}

				if (self_bs_idx == 0)
				{
					ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_in;
					ms[self_ms_idx].floor_idx = 0;
				}
			}
			else
			{
				cout << "Indoor state ERROR! - LOS Prob" << endl;
			}

			////////////////////////////////// Set OUT2IN case /////////////////////////////
			if (Indoor == true)
			{
				Propagation = OUT2IN_propagation; //// OUT2IN
			}

			// LOS = false;
			/// Propagation = NLOS_propagation;
			// cout << "LOS = false" << endl;
		}
		//////////////////////////////////////// UMa
		else
		{
			/////////////////////// outdoor
			if (Indoor == false)
			{
				Real C = 0;

				if (ms_height_out <= 13)
				{
					C = 0;
				}
				else
				{
					C = pow(((ms_height_out - 13.) / 10.), 1.5);
				}

				if (distance <= 18.)
				{
					LOS_prob = 1;
				}
				else if (distance > 18.)
				{
					LOS_prob = ((18. / distance) + exp(-1 * distance / 63.) * (1 - (18. / distance))) * (1 + C * (5 / 4) * pow((distance / 100), 3) * exp(-1 * distance / 150.));
				}

				Real P = randnum.u();

				if (P < LOS_prob)
				{
					LOS = true; /// LOS
					Propagation = LOS_propagation;
				}
				else
				{
					LOS = false; /// NLOS
					Propagation = NLOS_propagation;
				}

				if (self_bs_idx == 0)
				{
					ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_out;
					ms[self_ms_idx].floor_idx = 0;
					ms[self_ms_idx].ms_d_in = 0.;
				}
			}
			///////////////////////////// Indoor
			else if (Indoor == true)
			{

				Real d_in;
				Real d_out;

				d_in = 25. * randnum.u(); //// UMa, UMi = 0~25, RMa = 0~10
				d_out = distance - d_in;

				Real C = 0;

				if (ms_height_in <= 13)
				{
					C = 0;
				}
				else
				{
					C = pow(((ms_height_in - 13.) / 10.), 1.5);
				}

				if (d_out <= 18)
				{
					LOS_prob = 1;
				}
				else if (d_out > 18)
				{
					LOS_prob = ((18. / d_out) + exp(-1 * d_out / 63.) * (1 - (18. / d_out))) * (1 + C * (5 / 4) * pow(d_out / 100, 3) * exp(-1 * d_out / 150.));
				}

				Real P = randnum.u();
				if (P < LOS_prob)
				{
					LOS = true;
					Propagation = LOS_propagation;
				}
				else
				{
					LOS = false;
					Propagation = NLOS_propagation;
				}

				if (self_bs_idx == 0)
				{
					ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_in;
					ms[self_ms_idx].floor_idx = int(n) - 1;
					ms[self_ms_idx].ms_d_in = d_in;
				}
			}
			else
			{
				cout << "Indoor state ERROR! - LOS Prob" << endl;
			}

			/////////////////////// Set OUT2IN case //////////////
			if (Indoor == true)
			{
				Propagation = OUT2IN_propagation; //// OUT2IN
			}

			/// For SSP test

			// LOS = true;
			// Propagation = LOS_propagation;
			// cout << "LOS = true" << endl;

			// LOS = true;
			// Propagation = LOS_propagation;
			// cout << "LOS = true" << endl;
		}
	}
	else if (TYPE == 13) // RMa
	{
		///////////////////////////////// outdoor
		if (Indoor == false)
		{
			if (distance <= 10)
			{
				LOS_prob = 1;
			}
			else if (distance > 10.)
			{
				LOS_prob = exp(-1 * (distance - 10.) / 1000.);
			}

			Real P = randnum.u();
			if (P < LOS_prob)
			{
				LOS = true;
				Propagation = LOS_propagation;
			}
			else
			{
				LOS = false;
				Propagation = NLOS_propagation;
			}

			ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_out;
			ms[self_ms_idx].floor_idx = 0;
		}
		///////////////////////////////////////// Indoor
		else if (Indoor == true)
		{

			Real d_in;
			Real d_out;

			d_in = MIN(10. * randnum.u(), 10. * randnum.u()); //// UMa, UMi = 0~25, RMa = 0~10
			d_out = distance - d_in;

			if (distance <= 10)
			{
				LOS_prob = 1;
			}
			else if (distance > 10.)
			{
				LOS_prob = exp(-1 * (d_out - 10.) / 1000.);
			}

			Real P = randnum.u();
			if (P < LOS_prob)
			{
				LOS = true;
				Propagation = LOS_propagation;
			}
			else
			{
				LOS = false;
				Propagation = NLOS_propagation;
			}

			ms[self_ms_idx].MS_HEIGHT_FINAL = ms_height_in; // 1.5
			ms[self_ms_idx].floor_idx = 0;
		}
		else
		{
			cout << "Indoor state ERROR! - LOS Prob" << endl;
		}

		////////////////////////////////// Set OUT2IN case /////////////////////////////
		if (Indoor == true)
		{
			Propagation = OUT2IN_propagation; //// OUT2IN
		}

		// LOS = false;
		// Propagation = NLOS_propagation;
		// cout << "LOS = false" << endl;
	}
}
Real debugging_temp1 = 0;
Real debugging_temp2 = 0;

void CHANNEL::Set_PATHLOSS()
{
	////////////////////////////// PathLoss ////////////////////////////////
	// PathLoss
	//////////////////////////////////////////////////////// IMT-2020 /////////////////////////////////////////////////////////////////////////////////////////////////
	////// IMT 2020 EVAL Page 39 ~

	/////////////////////////////////////////// InH
	if (TYPE == 11)
	{
		Real distance_3d;
		distance_3d = sqrt(distance * distance + (_bs_height - ms_height) * (_bs_height - ms_height)); // 3D distance
		//////////////////////////////////////// Channel Model A
		if (Channel_Model_Type == 0)
		{
			//////////////////////////////////// LOS
			if (LOS == 1)
			{
				if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
				{
					sigma_SF = 3.;
					pathloss = 16.9 * log10(distance_3d) + 32.8 + 20 * log10(carrier_freq / 1000000000.);
				}
				else // 6GHz < fc <= 100 GHz
				{
					sigma_SF = 3.;
					pathloss = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(carrier_freq / 1000000000.);
				}
			}
			/////////////////////////////////// NLOS
			else if (LOS == 0)
			{
				if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
				{
					sigma_SF = 4.;
					pathloss = 43.3 * log10(distance_3d) + 11.5 + 20 * log10(carrier_freq / 1000000000.);
				}
				else // 6GHz < fc <= 100 GHz
				{
					sigma_SF = 8.03;
					pathloss = 43.3 * log10(distance_3d) + 11.5 + 20 * log10(carrier_freq / 1000000000.);
					pathloss_2 = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(carrier_freq / 1000000000.);
					pathloss = MAX(pathloss, pathloss_2);
				}
			}
		}
		//////////////////////////////////////// Channel Model
		else if (Channel_Model_Type == 1)
		{
			//////////////////////////////////// LOS
			if (LOS == 1)
			{
				sigma_SF = 3.;
				pathloss = 17.3 * log10(distance_3d) + 32.4 + 20 * log10(carrier_freq / 1000000000.);
			}
			//////////////////////////////////// NLOS
			else if (LOS == 0)
			{
				sigma_SF = 8.03;
				pathloss = 17.3 * log10(distance_3d) + 32.4 + 20 * log10(carrier_freq / 1000000000.);
				pathloss_2 = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(carrier_freq / 1000000000.);
				pathloss = MAX(pathloss, pathloss_2);
			}
		}
		// Channel_Model_Type == 2 for InH_ETRI channel
		else if (Channel_Model_Type == 2)
		{
			//////////////////////////////////// LOS
			if (LOS == 1)
			{
				if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
				{
					sigma_SF = 3.;
					pathloss = 16.9 * log10(distance_3d) + 32.8 + 20 * log10(carrier_freq / 1000000000.);
				}
				else // 6GHz < fc <= 100 GHz
				{
					sigma_SF = 2.39; // ETRI_Channel_Model
					pathloss = 32.4 + 17.3 * log10(distance_3d) + 20 * log10(carrier_freq / 1000000000.);
				}
			}
			/////////////////////////////////// NLOS
			else if (LOS == 0)
			{
				if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
				{
					sigma_SF = 4.;
					pathloss = 43.3 * log10(distance_3d) + 11.5 + 20 * log10(carrier_freq / 1000000000.);
				}
				else // 6GHz < fc <= 100 GHz
				{
					sigma_SF = 5.09; // ETRI_Channel_Model
					pathloss = 43.3 * log10(distance_3d) + 11.5 + 20 * log10(carrier_freq / 1000000000.);
					pathloss_2 = 38.3 * log10(distance_3d) + 17.3 + 24.9 * log10(carrier_freq / 1000000000.);
					pathloss = MAX(pathloss, pathloss_2);
				}
			}
		}

	}
	////////////////////////////////////////////////////// Dense_Urban
	else if (TYPE == 12)
	{
		// Real n = 0;
		// Real N = floor((5. * randnum.u()) + 4.); // 4 ~ 8
		// n = floor(((N)* randnum.u()) + 1.); // 1 ~ N

		// ms_height_in = 3 * (n - 1) + 1.5;
		int lc_Urban_env_Config_Type = Configuration_Type;
		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		{
			if ( self_bs_idx < simple_num_BS )
				lc_Urban_env_Config_Type = 1;
			else
				lc_Urban_env_Config_Type = 3;
		}


		Real ms_height_in_channel = 0;

		if (Indoor == true)
		{
			// ms_height_in_channel = ms_height_in;
			ms_height_in_channel = ms[self_ms_idx].MS_HEIGHT_FINAL;
		}
		else if (Indoor == false)
		{
			ms_height_in_channel = ms_height_out;
		}

		////////////////////////////////////////////////////// UMi
		if (lc_Urban_env_Config_Type == 3)
		{
			Real g;
			Real C;
			Real hE;

			if (distance > 18.)
			{
				g = (5 / 4) * pow((distance / 100), 3) * exp(-1 * distance / 150);
			}
			else
			{
				g = 0;
			}

			if (ms_height_in_channel < 13)
			{
				C = 0;
			}
			else
			{
				C = pow((ms_height_in_channel - 13.) / 10., 1.5) * g;
			}

			if (randnum.u() < (1 / (1 + C)))
			{
				hE = 1;
			}
			else
			{
				hE = floor((ms_height_in_channel - 1.5 - 9.) / 3. * randnum.u()) * 3. + 12.;
			}

			Real d_BP = 4 * (ms_height_in_channel - hE) * (_bs_height - hE) * carrier_freq / light_speed;
			Real distance_3d;
			distance_3d = sqrt(distance * distance + (_bs_height - ms_height_in_channel) * (_bs_height - ms_height_in_channel)); // 3D distance

			/////////////////////////////////////////////// Channel Model A
			if (Channel_Model_Type == 0)
			{
				//////////////////////////////// LOS
				if (LOS == 1)
				{
					if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 3.;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 3.;
							pathloss = 40. * log10(distance_3d) + 28. + 20. * log10(carrier_freq / 1000000000.) - 9. * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
					}
					else // 6GHz < fc <= 100 GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 4.;
							pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 4.;
							pathloss = 40. * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
					}
				}
				/////////////////////////////// NLOS
				else if (LOS == 0)
				{
					if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 3.;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 3.;
							pathloss = 40. * log10(distance_3d) + 28. + 20. * log10(carrier_freq / 1000000000.) - 9. * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
						pathloss_2 = 36.7 * log10(distance_3d) + 22.7 + 26. * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);

						pathloss = MAX(pathloss, pathloss_2);
					}
					else // 6GHz < fc <= 100 GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 4.;
							pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 4.;
							pathloss = 40. * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
						pathloss_2 = 35.3 * log10(distance_3d) + 22.4 + 21.3 * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);

						pathloss = MAX(pathloss, pathloss_2);
					}
				}
			}
			////////////////////////////////////////////////////////////// Channel Model B
			else if (Channel_Model_Type == 1)
			{
				/////////////////////////////////////////////////////// LOS
				if ( LOS == 1 )
				{
					if (distance <= d_BP)
					{
						sigma_SF = 4.;
						pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						sigma_SF = 4.;
						pathloss = 40 * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
				}
				//////////////////////////////////////////////// NLOS
				else if ( LOS == 0 )
				{
					if (distance <= d_BP)
					{
						sigma_SF = 4.;
						pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						sigma_SF = 4.;
						pathloss = 40 * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
					pathloss_2 = 22.4 + 35.3 * log10(distance_3d) + 21.3 * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);

					pathloss = MAX(pathloss, pathloss_2);
				}
			}
			else if (Channel_Model_Type == 2)
			{
				//////////////////////////////// LOS
				if (LOS == 1)
				{
					if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 3.;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 3.;
							pathloss = 40. * log10(distance_3d) + 28. + 20. * log10(carrier_freq / 1000000000.) - 9. * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
					}
					else // 6GHz < fc <= 100 GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 1.96; // ETRI_Channel_Model
							pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 1.96; // ETRI_Channel_Model
							pathloss = 40. * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
					}
				}
				/////////////////////////////// NLOS
				else if (LOS == 0)
				{
					if (carrier_freq <= 6 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 3.;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 3.;
							pathloss = 40. * log10(distance_3d) + 28. + 20. * log10(carrier_freq / 1000000000.) - 9. * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
						pathloss_2 = 36.7 * log10(distance_3d) + 22.7 + 26. * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);

						pathloss = MAX(pathloss, pathloss_2);
					}
					else // 6GHz < fc <= 100 GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 5.98; // ETRI_Channel_Model
							pathloss = 32.4 + 21. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 5.98; // ETRI_Channel_Model
							pathloss = 40. * log10(distance_3d) + 32.4 + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}
						pathloss_2 = 35.3 * log10(distance_3d) + 22.4 + 21.3 * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);
						pathloss = MAX(pathloss, pathloss_2);
					}
				}
			}

		}
		////////////////////////////////////////////////////// UMa
		else
		{
			Real g;
			Real C;
			Real hE;

			if (Indoor == false)
			{
				debugging_temp1 = 1;
			}

			if (distance > 18.)
			{
				g = (5 / 4) * pow((distance / 100), 3) * exp(-1 * distance / 150);
			}
			else
			{
				g = 0;
			}

			if (ms_height_in_channel < 13)
			{
				C = 0;
			}
			else
			{
				C = pow(((ms_height_in_channel - 13.) / 10.), 1.5) * g;
			}

			if (randnum.u() < (1 / (1 + C)))
			{
				hE = 1;
			}
			else
			{

				hE = floor((ms_height_in_channel - 1.5 - 9.) / 3. * randnum.u()) * 3. + 12.;
			}

			Real d_BP = 4 * (ms_height_in_channel - hE) * (_bs_height - hE) * carrier_freq / light_speed;
			Real distance_3d;
			distance_3d = sqrt(distance * distance + (_bs_height - ms_height_in_channel) * (_bs_height - ms_height_in_channel)); // 3D distance

			/////////////////////////////////////////////// Channel Model A
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2)
			{
				//////////////////////////////// LOS
				if (LOS == 1)
				{
					if (distance <= d_BP)
					{
						sigma_SF = 4.;
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						sigma_SF = 4.;
						pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
				}
				/////////////////////////////// NLOS
				else if (LOS == 0)
				{
					if (carrier_freq <= 6000000000.) ///// 0.5GHz <= fc <= 6GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 6.;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 6.;
							pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}

						pathloss_2 = 161.04 - 7.1 * log10(20) + 7.5 * log10(20) - (24.37 - 3.7 * (20 / _bs_height) * (20 / _bs_height)) * log10(_bs_height) + (43.42 - 3.1 * log10(_bs_height)) * (log10(distance_3d) - 3) + 20 * log10(carrier_freq / 1000000000.) - (3.2 * (pow(log10(17.625), 2)) - 4.97) - 0.6 * (ms_height_in_channel - 1.5);

						pathloss = MAX(pathloss, pathloss_2);
					}
					else // 6GHz < fc <= 100 GHz
					{
						if (distance <= d_BP)
						{
							sigma_SF = 6;
							pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
						}
						else if (d_BP < distance)
						{
							sigma_SF = 6;
							pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
						}

						pathloss_2 = 13.54 + 39.08 * log10(distance_3d) + 20 * log10(carrier_freq / 1000000000.) - 0.6 * (ms_height_in_channel - 1.5);

						pathloss = MAX(pathloss, pathloss_2);
					}
				}

				if (Propagation == OUT2IN_propagation)
				{
					sigma_SF = 7;
				}
			}
			////////////////////////////////////////////////////////////// Channel Model B
			else if (Channel_Model_Type == 1)
			{
				/////////////////////////////////////////////////////// LOS
				if (LOS == 1)
				{
					if (distance <= d_BP)
					{
						sigma_SF = 4.;
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						sigma_SF = 4.;
						pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
				}
				//////////////////////////////////////////////// NLOS
				else if (LOS == 0)
				{
					if (distance <= d_BP)
					{
						sigma_SF = 6.;
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						sigma_SF = 6.;
						pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}

					pathloss_2 = 13.54 + 39.08 * log10(distance_3d) + 20 * log10(carrier_freq / 1000000000.) - 0.6 * (ms_height_in_channel - 1.5);

					pathloss = MAX(pathloss, pathloss_2);
				}

				if (Propagation == OUT2IN_propagation)
				{
					sigma_SF = 7;
				}
			}
		}
	}
	/////////////////////////////////////////////////////////////////////////////// RMa
	else if (TYPE == 13)
	{

		Real ms_height_in_channel = 0;

		/////////////////////////////////////////////////////////////// Indoor
		if (Indoor == true)
		{
			ms_height_in_channel = ms_height_in;
		}
		///////////////////////////////////////////////////////// Outdoor
		else if (Indoor == false)
		{
			ms_height_in_channel = ms_height_out;
		}

		Real d_BP = 2 * pi * (ms_height_in_channel) * (_bs_height)*carrier_freq / light_speed;
		Real distance_3d;
		distance_3d = sqrt(distance * distance + (_bs_height - ms_height_in_channel) * (_bs_height - ms_height_in_channel)); // 3D distance

		/////////////////////////////////////////////////////////// Channel Model A
		if (Channel_Model_Type == 0)
		{
			Real hh = 5.;
			Real WW = 20.;

			///////////////////////////////////LOS
			if (LOS == 1)
			{
				if (distance <= d_BP)
				{
					sigma_SF = 4.;
					pathloss = 20 * log10(40 * pi * distance_3d * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(distance_3d) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * distance_3d;
				}
				else if (d_BP < distance && distance < 21000.)
				{
					sigma_SF = 6.;
					pathloss_2 = 20 * log10(40 * pi * d_BP * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(d_BP) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * d_BP;

					pathloss = pathloss_2 + 40 * log10(distance_3d / d_BP);
				}
			}
			//////////////////////////////// NLOS
			else if (LOS == 0)
			{
				sigma_SF = 8.;
				pathloss = 161.04 - 7.1 * log10(WW) + 7.5 * log10(hh) - (24.37 - 3.7 * (pow((hh / _bs_height), 2))) * log10(_bs_height) + (43.42 - 3.1 * log10(_bs_height)) * (log10(distance_3d) - 3) + 20 * log10(carrier_freq / 1000000000.) - (3.2 * (log10(11.75 * ms_height_in_channel)) * (log10(11.75 * ms_height_in_channel)) - 4.97);

				if (Configuration_Type == 2) /// Rural config C, LMLC
				{
					if (distance <= d_BP)
					{
						sigma_SF = 4.;
						pathloss_2 = 20 * log10(40 * pi * distance_3d * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(distance_3d) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * distance_3d;
					}
					else if (d_BP < distance && distance < 21000.)
					{
						sigma_SF = 6.;
						pathloss_2 = 20 * log10(40 * pi * d_BP * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(d_BP) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * d_BP;

						pathloss_2 = pathloss_2 + 40 * log10(distance_3d / d_BP);
					}
					pathloss = MAX(pathloss_2, pathloss - 12);
				}

				sigma_SF = 8.;
			}

			if (Propagation == OUT2IN_propagation)
			{
				sigma_SF = 8;
			}
		}
		///////////////////////////////////////////////////// Channel Model B
		else if (Channel_Model_Type == 1)
		{
			Real hh = 5.;
			Real WW = 20.;

			////////////////////////////////////////// LOS
			if (LOS == 1)
			{
				if (distance <= d_BP)
				{
					sigma_SF = 4.;
					pathloss = 20 * log10(40 * pi * distance_3d * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(distance_3d) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * distance_3d;
				}
				else if (d_BP < distance)
				{
					sigma_SF = 6.;
					pathloss_2 = 20 * log10(40 * pi * d_BP * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(d_BP) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * d_BP;

					pathloss = pathloss_2 + 40 * log10(distance_3d / d_BP);
				}
			}
			/////////////////////////////////////////// NLOS
			else if (LOS == 0)
			{

				sigma_SF = 8.;
				pathloss = 161.04 - 7.1 * log10(WW) + 7.5 * log10(hh) - (24.37 - 3.7 * (pow((hh / _bs_height), 2))) * log10(_bs_height) + (43.42 - 3.1 * log10(_bs_height)) * (log10(distance_3d) - 3) + 20 * log10(carrier_freq / 1000000000.) - (3.2 * (log10(11.75 * ms_height_in_channel)) * (log10(11.75 * ms_height_in_channel)) - 4.97);

				if (distance <= d_BP)
				{
					sigma_SF = 4.;
					pathloss_2 = 20 * log10(40 * pi * distance_3d * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(distance_3d) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * distance_3d;
				}
				else if (d_BP < distance)
				{
					sigma_SF = 6.;
					pathloss_2 = 20 * log10(40 * pi * d_BP * (carrier_freq / 1000000000.) / 3) + MIN(0.03 * pow(hh, 1.72), 10) * log10(d_BP) - MIN(0.044 * pow(hh, 1.72), 14.77) + 0.002 * log10(hh) * d_BP;

					pathloss_2 = pathloss_2 + 40 * log10(distance_3d / d_BP);
				}

				if (Configuration_Type == 2) /// Rural config C, LMLC
				{
					pathloss = MAX(pathloss_2, pathloss - 12);
				}
				else
				{
					pathloss = MAX(pathloss, pathloss_2);
				}

				sigma_SF = 8.;
			}

			if (Propagation == OUT2IN_propagation)
			{
				sigma_SF = 8;
			}
		}
	}
}


Real normal(Real mu, Real sigma)
{
	Real U1, U2, W, mult;
	static Real X1, X2;
	static int call = 0;

	if (call == 1)
	{
		call = !call;
		return (mu + sigma * (Real)X2);
	}

	do
	{
		U1 = -1 + randnum.u() * 2;
		U2 = -1 + randnum.u() * 2;
		W = pow(U1, 2) + pow(U2, 2);
	} while (W >= 1 || W == 0);

	mult = sqrt((-2 * log(W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return (mu + sigma * (Real)X1);
}

void CHANNEL::Set_Channel_Parameters()
{

	////////////////////////// Parameter setting /////////////////////////////////////////////////////
	//////////////////////////////////////////////////////// IMT-2020 /////////////////////////////////////////////////////////////////////////////////////////////////
	////// M.2412-0 page 62

	//////////////////////////////////////////////////////////////////// InH
	if (TYPE == 11)
	{
		//////////////////////////////////////////////////// Channel Model A
		if (Channel_Model_Type == 0)
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			{

				XPR_mean[LOS_propagation] = 11;
				XPR_mean[NLOS_propagation] = 10;
				XPR_std[LOS_propagation] = 4;
				XPR_std[NLOS_propagation] = 4;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); // linear

					mu_DS = -7.70;
					sigma_DS = 0.18;
					DS = pow(10, normal(-7.70, 0.18)); // = sigma_tau in M.2135

					mu_ASD = 1.6;
					sigma_ASD = 0.18;
					mu_ASA = 1.62;
					sigma_ASA = 0.22;

					mu_ZSA = 1.22;
					sigma_ZSA = 0.23;
					mu_ZSD = 1.02;
					sigma_ZSD = 0.41;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 8;
					cluster_ZSA = 9;

					XPR = normal(11, 4);
					cluster_shadowing = 6;
				}
				else if (LOS == 0) //// NLOS
				{
					num_path = 19; // Number of clusters
					r_tau = 3;

					K_factor = -1; // N/A
					mu_K_factor = -1;
					sigma_K_factor = -1;

					mu_DS = -7.41;
					sigma_DS = 0.14;
					DS = pow(10, normal(-7.41, 0.14)); // = sigma_tau in M.2135

					mu_ASD = 1.62;
					sigma_ASD = 0.25;
					mu_ASA = 1.77;
					sigma_ASA = 0.16;

					mu_ZSA = 1.26;
					sigma_ZSA = 0.67;
					mu_ZSD = 1.08;
					sigma_ZSD = 0.36;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 11;
					cluster_ZSA = 9;

					XPR = normal(10, 4);
					cluster_shadowing = 3;
				}
			}
			else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
			// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
			{
				XPR_mean[LOS_propagation] = 11;
				XPR_mean[NLOS_propagation] = 10;
				// XPR_mean[OUT2IN_propagation] = 9;

				XPR_std[LOS_propagation] = 4;
				XPR_std[NLOS_propagation] = 4;
				// XPR_std[OUT2IN_propagation] = 5;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

					mu_DS = -0.01 * log10(1 + (carrier_freq / 1000000000)) - 7.692;
					sigma_DS = 0.18;
					DS = pow(10, normal(mu_DS, 0.18)); // = sigma_tau in M.2135

					mu_ASD = 1.6;
					sigma_ASD = 0.18;
					mu_ASA = -0.19 * log10(1 + (carrier_freq / 1000000000)) + 1.781;
					sigma_ASA = 0.12 * log10(1 + (carrier_freq / 1000000000)) + 0.119;

					mu_ZSA = -0.26 * log10(1 + (carrier_freq / 1000000000)) + 1.44;
					sigma_ZSA = -0.04 * log10(1 + (carrier_freq / 1000000000)) + 0.264;
					mu_ZSD = -1.43 * log10(1 + (carrier_freq / 1000000000)) + 2.228;
					sigma_ZSD = 0.13 * log10(1 + (carrier_freq / 1000000000)) + 0.3;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 8;
					cluster_ZSA = 9;

					XPR = normal(11, 4);
					cluster_shadowing = 6;
				}
				else if (LOS == 0) //// NLOS
				{
					num_path = 19; // Number of clusters
					r_tau = 3;

					K_factor = -1; // N/A
					mu_K_factor = -1;
					sigma_K_factor = -1;

					mu_DS = -0.28 * log10(1 + (carrier_freq / 1000000000)) - 7.173;
					sigma_DS = 0.1 * log10(1 + (carrier_freq / 1000000000)) + 0.055;
					DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

					mu_ASD = 1.62;
					sigma_ASD = 0.25;
					mu_ASA = -0.11 * log10(1 + (carrier_freq / 1000000000)) + 1.863;
					sigma_ASA = 0.12 * log10(1 + (carrier_freq / 1000000000)) + 0.059;

					mu_ZSA = -0.15 * log10(1 + (carrier_freq / 1000000000)) + 1.387;
					sigma_ZSA = -0.09 * log10(1 + (carrier_freq / 1000000000)) + 0.746;
					mu_ZSD = 1.08;
					sigma_ZSD = 0.36;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 11;
					cluster_ZSA = 9;

					XPR = normal(10, 4);
					cluster_shadowing = 3;
				}
			}
			else
			{
				cout << "carrier_freq ERROR! - Set channel parameter" << endl;
			}
		}
		//////////////////////////////////////////////////// Channel Model B
		else if (Channel_Model_Type == 1)
		{
			if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
			// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
			{
				Real fc;

				if (carrier_freq < 6000000000) // below 6GHz
				{
					fc = 6;
				}
				else
				{
					fc = (carrier_freq / 1000000000);
				}

				XPR_mean[LOS_propagation] = 11;
				XPR_mean[NLOS_propagation] = 10;
				// XPR_mean[OUT2IN_propagation] = 9;

				XPR_std[LOS_propagation] = 4;
				XPR_std[NLOS_propagation] = 4;
				// XPR_std[OUT2IN_propagation] = 5;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

					mu_DS = -0.01 * log10(1 + fc) - 7.692;
					sigma_DS = 0.18;
					DS = pow(10, normal(mu_DS, 0.18)); // = sigma_tau in M.2135

					mu_ASD = 1.6;
					sigma_ASD = 0.18;
					mu_ASA = -0.19 * log10(1 + fc) + 1.781;
					sigma_ASA = 0.12 * log10(1 + fc) + 0.119;

					mu_ZSA = -0.26 * log10(1 + fc) + 1.44;
					sigma_ZSA = -0.04 * log10(1 + fc) + 0.264;
					mu_ZSD = -1.43 * log10(1 + fc) + 2.228;
					sigma_ZSD = 0.13 * log10(1 + fc) + 0.3;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 8;
					cluster_ZSA = 9;

					XPR = normal(11, 4);
					cluster_shadowing = 6;
				}
				else if (LOS == 0) //// NLOS
				{
					num_path = 19; // Number of clusters
					r_tau = 3;

					K_factor = -1; // N/A
					mu_K_factor = -1;
					sigma_K_factor = -1;

					mu_DS = -0.28 * log10(1 + fc) - 7.173;
					sigma_DS = 0.1 * log10(1 + fc) + 0.055;
					DS = pow(10, normal(-7.41, 0.14)); // = sigma_tau in M.2135

					mu_ASD = 1.62;
					sigma_ASD = 0.25;
					mu_ASA = -0.11 * log10(1 + fc) + 1.863;
					sigma_ASA = 0.12 * log10(1 + fc) + 0.059;

					mu_ZSA = -0.15 * log10(1 + fc) + 1.387;
					sigma_ZSA = -0.09 * log10(1 + fc) + 0.746;
					mu_ZSD = 1.08;
					sigma_ZSD = 0.36;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 11;
					cluster_ZSA = 9;

					XPR = normal(10, 4);
					cluster_shadowing = 3;
				}
			}
			else
			{
				cout << "carrier_freq ERROR! - Set channel parameter" << endl;
			}
		}
		// Channel_Model_Type == 2 for InH_ETRI channel
		else if (Channel_Model_Type == 2)
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			{

				XPR_mean[LOS_propagation] = 11;
				XPR_mean[NLOS_propagation] = 10;
				XPR_std[LOS_propagation] = 4;
				XPR_std[NLOS_propagation] = 4;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); // linear

					mu_DS = -7.70;
					sigma_DS = 0.18;
					DS = pow(10, normal(-7.70, 0.18)); // = sigma_tau in M.2135

					mu_ASD = 1.6;
					sigma_ASD = 0.18;
					mu_ASA = 1.62;
					sigma_ASA = 0.22;

					mu_ZSA = 1.22;
					sigma_ZSA = 0.23;
					mu_ZSD = 1.02;
					sigma_ZSD = 0.41;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 8;
					cluster_ZSA = 9;

					XPR = normal(11, 4);
					cluster_shadowing = 6;
				}
				else if (LOS == 0) //// NLOS
				{
					num_path = 19; // Number of clusters
					r_tau = 3;

					K_factor = -1; // N/A
					mu_K_factor = -1;
					sigma_K_factor = -1;

					mu_DS = -7.41;
					sigma_DS = 0.14;
					DS = pow(10, normal(-7.41, 0.14)); // = sigma_tau in M.2135

					mu_ASD = 1.62;
					sigma_ASD = 0.25;
					mu_ASA = 1.77;
					sigma_ASA = 0.16;

					mu_ZSA = 1.26;
					sigma_ZSA = 0.67;
					mu_ZSD = 1.08;
					sigma_ZSD = 0.36;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 11;
					cluster_ZSA = 9;

					XPR = normal(10, 4);
					cluster_shadowing = 3;
				}
			}
			else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
			// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
			{
				XPR_mean[LOS_propagation] = 11;
				XPR_mean[NLOS_propagation] = 10;
				// XPR_mean[OUT2IN_propagation] = 9;

				XPR_std[LOS_propagation] = 4;
				XPR_std[NLOS_propagation] = 4;
				// XPR_std[OUT2IN_propagation] = 5;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

					/*-------------- modified by jhnoh 230914 -----------------*/
					mu_DS    = -8.41; // ETRI_Channel_Model
					sigma_DS = 0.52;  // ETRI_Channel_Model
					DS = pow(10, normal(mu_DS, 0.18)); // = sigma_tau in M.2135
					/*----------------------------------------------------------*/

					mu_ASD = 1.6;
					sigma_ASD = 0.18;

                    /*-------------- modified by jhnoh 230914 -----------------*/
					mu_ASA    = 1.52;
					sigma_ASA = 0.16;
					/*----------------------------------------------------------*/

					mu_ZSA = -0.26 * log10(1 + (carrier_freq / 1000000000)) + 1.44;
					sigma_ZSA = -0.04 * log10(1 + (carrier_freq / 1000000000)) + 0.264;
					mu_ZSD = -1.43 * log10(1 + (carrier_freq / 1000000000)) + 2.228;
					sigma_ZSD = 0.13 * log10(1 + (carrier_freq / 1000000000)) + 0.3;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 8;
					cluster_ZSA = 9;

					XPR = normal(11, 4);
					cluster_shadowing = 6;
				}
				else if (LOS == 0) //// NLOS
				{
					num_path = 19; // Number of clusters
					r_tau = 3;

					K_factor = -1; // N/A
					mu_K_factor = -1;
					sigma_K_factor = -1;

					/*-------------- modified by jhnoh 230914 -----------------*/
					mu_DS    = -7.9;
					sigma_DS = 0.38;
					DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135
					/*----------------------------------------------------------*/

					mu_ASD = 1.62;
					sigma_ASD = 0.25;

					/*-------------- modified by jhnoh 230914 -----------------*/
					mu_ASA    = 1.80;
					sigma_ASA = 0.15;
					/*----------------------------------------------------------*/

					mu_ZSA = -0.15 * log10(1 + (carrier_freq / 1000000000)) + 1.387;
					sigma_ZSA = -0.09 * log10(1 + (carrier_freq / 1000000000)) + 0.746;
					mu_ZSD = 1.08;
					sigma_ZSD = 0.36;
					mu_offset_ZOD = 0;

					cluster_DS = -1; /// N/A
					cluster_ASD = 5;
					cluster_ASA = 11;
					cluster_ZSA = 9;

					XPR = normal(10, 4);
					cluster_shadowing = 3;
				}
			}
			else
			{
				cout << "carrier_freq ERROR! - Set channel parameter" << endl;
			}
		}		
	}
	//////////////////////////////////////////////////////////////////////// Dense_Urban
	else if (TYPE == 12)
	{
		int lc_Urban_env_Config_Type = Configuration_Type;
		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		{
			if ( self_bs_idx < simple_num_BS )
				lc_Urban_env_Config_Type = 1;
			else
				lc_Urban_env_Config_Type = 3;
		}

		////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////// UMi
		if (lc_Urban_env_Config_Type == 3)
		{
			/////////////////////////////////////////////////////////////////////// Channel Model A
			if (Channel_Model_Type == 0)
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 3.2;

						// K_factor = normal(9, 5);
						mu_K_factor = 9.;
						sigma_K_factor = 5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -7.19;
						sigma_DS = 0.40;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.2;
						sigma_ASD = 0.43;
						mu_ASA = 1.75;
						sigma_ASA = 0.19;

						mu_ZSA = 0.6;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						cluster_DS = -1; /// N/A
						cluster_ASD = 3;
						cluster_ASA = 17;
						cluster_ZSA = 7;

						XPR = normal(9, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 19; // Number of clusters
						r_tau = 3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.89;
						sigma_DS = 0.54;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.41;
						sigma_ASD = 0.17;
						mu_ASA = 1.84;
						sigma_ASA = 0.15;

						mu_ZSA = 0.88;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * MAX((ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height), 0) + 0.9);
						sigma_ZSD = 0.6;
						mu_offset_ZOD = -1 * pow(10, -0.55 * log10(MAX(10, distance)) + 1.6);

						cluster_DS = -1; /// N/A
						cluster_ASD = 10;
						cluster_ASA = 22;
						cluster_ZSA = 7;

						XPR = normal(8, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;

						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.75);
							sigma_ZSD = 0.4;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * MAX((ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height), 0) + 0.9);
							sigma_ZSD = 0.6;
							mu_offset_ZOD = -1 * pow(10, -0.55 * log10(MAX(10, distance)) + 1.6);
						}

						cluster_DS = -1; /// N/A
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 11);
						cluster_shadowing = 4;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 3;

						// K_factor = normal(9, 3.5);
						mu_K_factor = 9.;
						sigma_K_factor = 5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -7.14 - 0.24 * log10(1 + (carrier_freq / 1000000000));
						sigma_DS = 0.38;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.21 - 0.05 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD = 0.41;
						mu_ASA = 1.73 - 0.08 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASA = 0.28 + 0.014 * log10(1 + (carrier_freq / 1000000000));

						mu_ZSA = 0.73 - 0.1 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA = 0.34 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
						sigma_ZSD = 0.35;
						mu_offset_ZOD = 0;

						cluster_DS = 5;
						cluster_ASD = 3;
						cluster_ASA = 17;
						cluster_ZSA = 7;

						XPR = normal(9, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20; // Number of clusters
						r_tau = 2.3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.83 - 0.24 * log10(1 + (carrier_freq / 1000000000));
						sigma_DS = 0.28 + 0.16 * log10(1 + (carrier_freq / 1000000000));
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.53 - 0.23 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD = 0.33 + 0.11 * log10(1 + (carrier_freq / 1000000000));
						mu_ASA = 1.81 - 0.08 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASA = 0.3 + 0.05 * log10(1 + (carrier_freq / 1000000000));

						mu_ZSA = 0.92 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA = 0.41 - 0.07 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
						sigma_ZSD = 0.35;
						mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);

						cluster_DS = 11;
						cluster_ASD = 10;
						cluster_ASA = 22;
						cluster_ZSA = 7;

						XPR = normal(8, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;
						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);
						}

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 20;
						cluster_ZSA = 6;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else
				{
					cout << "carrier_freq ERROR! - Set channel parameter" << endl;
				}
			}
			////////////////////////////////////////////////////////// Channel Model B
			else if (Channel_Model_Type == 1)
			{
				if (carrier_freq > 500000000) /// 0.5GHz < fc <= 100GHz
				// if (carrier_freq > 500000000 && carrier_freq <= 100000000000) /// 0.5GHz < fc <= 100GHz
				{
					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 3;

						// K_factor = normal(9, 5);
						mu_K_factor = 9.;
						sigma_K_factor = 5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -7.14 - 0.24 * log10(1 + (carrier_freq / 1000000000));
						sigma_DS = 0.38;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.21 - 0.05 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD = 0.41;
						mu_ASA = 1.73 - 0.08 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASA = 0.28 + 0.014 * log10(1 + (carrier_freq / 1000000000));

						mu_ZSA = 0.73 - 0.1 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA = 0.34 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
						sigma_ZSD = 0.35;
						mu_offset_ZOD = 0;

						cluster_DS = 5;
						cluster_ASD = 3;
						cluster_ASA = 17;
						cluster_ZSA = 7;

						XPR = normal(9, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 19; // Number of clusters
						r_tau = 2.1;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.83 - 0.24 * log10(1 + (carrier_freq / 1000000000));
						sigma_DS = 0.28 + 0.16 * log10(1 + (carrier_freq / 1000000000));
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.53 - 0.23 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD = 0.33 + 0.11 * log10(1 + (carrier_freq / 1000000000));
						mu_ASA = 1.81 - 0.08 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASA = 0.3 + 0.05 * log10(1 + (carrier_freq / 1000000000));

						mu_ZSA = 0.92 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA = 0.41 - 0.07 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
						sigma_ZSD = 0.35;
						mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.25;
						sigma_ASA = 0.42;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;

						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);
						}

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else
				{
					cout << "carrier_freq ERROR! - Set channel parameter" << endl;
				}
			}
			/*-------------------------- For ETRI Channel model ----------------------------*/
		    else if (Channel_Model_Type == 2) 
			{
				if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 3;

						// K_factor = normal(9, 3.5);
						mu_K_factor    = 9.;
						sigma_K_factor = 5;
						K_factor       = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						//---------------------------------------
						mu_DS          = -8.51;
						sigma_DS       =  0.75;
						DS             = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135
						//---------------------------------------

						mu_ASD         = 1.21 - 0.05 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD      = 0.41;
						
						//---------------------------------------
						mu_ASA         = 1.06;
						sigma_ASA      = 0.33;
						//---------------------------------------

						mu_ZSA         = 0.73 - 0.1 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA      = 0.34 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD         = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
						sigma_ZSD      = 0.35;
						mu_offset_ZOD  = 0;

						cluster_DS     = 5;
						cluster_ASD    = 3;
						cluster_ASA    = 17;
						cluster_ZSA    = 7;

						XPR = normal(9, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20; // Number of clusters
						r_tau = 2.3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

                        //---------------------------------------
						mu_DS     = -7.66;
						sigma_DS  = 0.43;
						DS        = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135
						//---------------------------------------

						mu_ASD    = 1.53 - 0.23 * log10(1 + (carrier_freq / 1000000000));
						sigma_ASD = 0.33 + 0.11 * log10(1 + (carrier_freq / 1000000000));

						//---------------------------------------
						mu_ASA    = 1.36;
						sigma_ASA = 0.32;
						//---------------------------------------

						mu_ZSA = 0.92 - 0.04 * log10(1 + (carrier_freq / 1000000000));
						sigma_ZSA = 0.41 - 0.07 * log10(1 + (carrier_freq / 1000000000));
						mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
						sigma_ZSD = 0.35;
						mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);

						cluster_DS = 11;
						cluster_ASD = 10;
						cluster_ASA = 22;
						cluster_ZSA = 7;

						XPR = normal(8, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;
						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.21, -14.8 * (distance / 1000) + 0.01 * abs(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height) + 0.83);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -3.1 * (distance / 1000) + 0.01 * MAX(ms[self_ms_idx].MS_HEIGHT_FINAL - _bs_height, 0) + 0.2);
							sigma_ZSD = 0.35;
							mu_offset_ZOD = -1 * pow(10, -1.5 * log10(MAX(10, distance)) + 3.3);
						}

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 20;
						cluster_ZSA = 6;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else
				{
					cout << "carrier_freq ERROR! - Set channel parameter" << endl;
					getchar();
				}
			}
		}
		/////////////////////////////////////////////////////////////////////////// UMa
		else
		{
			/////////////////////////////////////////////////////////////////////// Channel Model A
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2)
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 2.5;

						// K_factor = normal(9, 3.5);
						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -7.03;
						sigma_DS = 0.66;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.15;
						sigma_ASD = 0.28;
						mu_ASA = 1.81;
						sigma_ASA = 0.2;

						mu_ZSA = 0.95;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						cluster_DS = -1; /// N/A
						cluster_ASD = 5;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20; // Number of clusters
						r_tau = 2.3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.44;
						sigma_DS = 0.39;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.41;
						sigma_ASD = 0.28;
						mu_ASA = 1.87;
						sigma_ASA = 0.11;

						mu_ZSA = 1.26;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = -1 * pow(10, -0.62 * log10(MAX(10, distance)) + 1.93 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

						cluster_DS = -1; /// N/A
						cluster_ASD = 2;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;

						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
							sigma_ZSD = 0.4;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
							sigma_ZSD = 0.49;
							mu_offset_ZOD = -1 * pow(10, -0.62 * log10(MAX(10, distance)) + 1.93 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

							///// -10���� -�� ������ ����. -10���� �Ǿ������� mu_offset_ZOD���� -1#IND�� ������ ��찡 �־...
						}

						/*
						if (self_ms_idx == 4)
						{
							cout << "Propagation = " << Propagation << endl;
							cout << "LOS = " << LOS << endl;
							cout << "distance = " << distance << endl;
							cout << "ms[self_ms_idx].MS_HEIGHT_FINAL = " << ms[self_ms_idx].MS_HEIGHT_FINAL << endl;
							cout << "mu_offset = " << mu_offset_ZOD << endl;
							getchar();
						}
						*/

						cluster_DS = -1; /// N/A
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 2.5;

						// K_factor = normal(9, 3.5);
						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -6.955 - 0.0963 * log10(carrier_freq / 1000000000);
						sigma_DS = 0.66;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.06 + 0.1114 * log10(carrier_freq / 1000000000);
						sigma_ASD = 0.28;
						mu_ASA = 1.81;
						sigma_ASA = 0.2;

						mu_ZSA = 0.95;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						cluster_DS = MAX(0.25, -3.4084 * log10(carrier_freq / 1000000000) + 6.5622);
						cluster_ASD = 5;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20; // Number of clusters
						r_tau = 2.3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.28 - 0.204 * log10(carrier_freq / 1000000000);
						sigma_DS = 0.39;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.5 - 0.1144 * log10(carrier_freq / 1000000000);
						sigma_ASD = 0.28;
						mu_ASA = 2.08 - 0.27 * log10(carrier_freq / 1000000000);
						sigma_ASA = 0.11;

						mu_ZSA = -0.3236 * log10(carrier_freq / 1000000000) + 1.512;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = 7.66 * log10(carrier_freq / 1000000000) - 5.96 - pow(10, (0.208 * log10(carrier_freq / 1000000000) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(carrier_freq / 1000000000) + 2.03);

						cluster_DS = MAX(0.25, -3.4084 * log10(carrier_freq / 1000000000) + 6.5622);
						cluster_ASD = 2;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.63;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;
						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
							sigma_ZSD = 0.4;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
							sigma_ZSD = 0.49;
							mu_offset_ZOD = 7.66 * log10(carrier_freq / 1000000000) - 5.96 - pow(10, (0.208 * log10(carrier_freq / 1000000000) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(carrier_freq / 1000000000) + 2.03);
						}

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else
				{
					cout << "carrier_freq ERROR! - Set channel parameter" << endl;
				}
			}
			////////////////////////////////////////////////////////// Channel Model B
			else if (Channel_Model_Type == 1)
			{

				if (carrier_freq > 500000000) /// 0.5GHz < fc <= 100GHz
				// if (carrier_freq > 500000000 && carrier_freq <= 100000000000) /// 0.5GHz < fc <= 100GHz
				{
					Real fc;

					if (carrier_freq < 6000000000) // below 6GHz
					{
						fc = 6;
					}
					else
					{
						fc = (carrier_freq / 1000000000);
					}

					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12; // Number of clusters
						r_tau = 2.5;

						// K_factor = normal(9, 3.5);
						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

						mu_DS = -6.955 - 0.0963 * log10(fc);
						sigma_DS = 0.66;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.06 + 0.1114 * log10(fc);
						sigma_ASD = 0.28;
						mu_ASA = 1.81;
						sigma_ASA = 0.2;

						mu_ZSA = 0.95;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = 5;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20; // Number of clusters
						r_tau = 2.3;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.28 - 0.204 * log10(fc);
						sigma_DS = 0.39;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.5 - 0.1144 * log10(fc);
						sigma_ASD = 0.28;
						mu_ASA = 2.08 - 0.27 * log10(fc);
						sigma_ASA = 0.11;

						mu_ZSA = -0.3236 * log10(fc) + 1.512;
						sigma_ZSA = 0.16;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = 7.66 * log10(fc) - 5.96 - pow(10, (0.208 * log10(fc) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = 2;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12; // Number of clusters
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

						mu_ASD = 1.25;
						sigma_ASD = 0.42;
						mu_ASA = 1.76;
						sigma_ASA = 0.16;

						mu_ZSA = 1.01;
						sigma_ZSA = 0.43;
						if (LOS == 1)
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
							sigma_ZSD = 0.4;
							mu_offset_ZOD = 0;
						}
						else
						{
							mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
							sigma_ZSD = 0.49;
							mu_offset_ZOD = 7.66 * log10(fc) - 5.96 - pow(10, (0.208 * log10(fc) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));
						}

						cluster_DS = 11;
						cluster_ASD = 5;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else
				{
					cout << "carrier_freq ERROR! - Set channel parameter" << endl;
				}
			}
		}
	}
	/////////////////////////////////////////////////// RMa
	else if (TYPE == 13)
	{
		if (carrier_freq > 500000000) /// RMm_A:0.5GHz < fc <= 6GHz, RMa_B:0.5GHz < fc <= 7GHz
		// if (carrier_freq > 500000000 && carrier_freq <= 7000000000) /// RMm_A:0.5GHz < fc <= 6GHz, RMa_B:0.5GHz < fc <= 7GHz
		{
			XPR_mean[LOS_propagation] = 12;
			XPR_mean[NLOS_propagation] = 7;
			XPR_mean[OUT2IN_propagation] = 7;

			XPR_std[LOS_propagation] = 4;
			XPR_std[NLOS_propagation] = 3;
			XPR_std[OUT2IN_propagation] = 3;

			if (Propagation == 1) // LOS
			{
				num_path = 11; // Number of clusters
				r_tau = 3.8;

				// K_factor = normal(7, 4);
				mu_K_factor = 7.;
				sigma_K_factor = 4.;
				K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); //

				mu_DS = -7.49;
				sigma_DS = 0.55;
				DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

				mu_ASD = 0.9;
				sigma_ASD = 0.38;
				mu_ASA = 1.52;
				sigma_ASA = 0.24;

				mu_ZSA = 0.47;
				sigma_ZSA = 0.4;
				mu_ZSD = MAX(-1, -0.17 * (distance / 1000.) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.22);
				sigma_ZSD = 0.34;
				mu_offset_ZOD = 0;

				cluster_DS = -1; /// N/A
				cluster_ASD = 2;
				cluster_ASA = 3;
				cluster_ZSA = 3;

				XPR = normal(12, 4);
				cluster_shadowing = 3;
			}
			else if (Propagation == 0) //// NLOS
			{
				num_path = 10; // Number of clusters
				r_tau = 1.7;

				K_factor = -1; // N/A
				mu_K_factor = -1;
				sigma_K_factor = -1;

				mu_DS = -7.43;
				sigma_DS = 0.48;
				DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

				mu_ASD = 0.95;
				sigma_ASD = 0.45;
				mu_ASA = 1.52;
				sigma_ASA = 0.13;

				mu_ZSA = 0.58;
				sigma_ZSA = 0.37;
				mu_ZSD = MAX(-1, -0.19 * (distance / 1000.) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.28);
				sigma_ZSD = 0.3;
				mu_offset_ZOD = atan((35. - 3.5) / distance) - atan((35. - 1.5) / distance);

				cluster_DS = -1; /// N/A
				cluster_ASD = 2;
				cluster_ASA = 3;
				cluster_ZSA = 3;

				XPR = normal(7, 3);
				cluster_shadowing = 3;
			}
			else if (Propagation == 2) ///// OUT2IN
			{
				num_path = 10; // Number of clusters
				r_tau = 1.7;

				K_factor = -1; // N/A
				mu_K_factor = -1;
				sigma_K_factor = -1;

				mu_DS = -7.47;
				sigma_DS = 0.24;
				DS = pow(10, normal(mu_DS, sigma_DS)); // = sigma_tau in M.2135

				mu_ASD = 0.67;
				sigma_ASD = 0.18;
				mu_ASA = 1.66;
				sigma_ASA = 0.21;

				mu_ZSA = 0.93;
				sigma_ZSA = 0.22;
				mu_ZSD = MAX(-1, -0.19 * (distance / 1000.) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.28);
				sigma_ZSD = 0.3;
				mu_offset_ZOD = atan((35. - 3.5) / distance) - atan((35. - 1.5) / distance);

				cluster_DS = -1; /// N/A
				cluster_ASD = 2;
				cluster_ASA = 3;
				cluster_ZSA = 3;

				XPR = normal(7, 3);
				cluster_shadowing = 3;
			}
		}
		else
		{
			cout << "carrier_freq ERROR! - Set channel parameter" << endl;
		}
	}
}

void CHANNEL::Set_DELAY()
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter -- delay
	///
	/////////////////////////////////////////////////

	Real K_factor_dB = K_dB;

	D = 0.7705 - 0.0433 * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB + 0.000017 * K_factor_dB * K_factor_dB * K_factor_dB;

	// delay = new Real[MAX_NUM_CLUSTERS];
	// delay_LOS = new Real[MAX_NUM_CLUSTERS];

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		delay[i] = 0.;
		delay_LOS[i] = 0.;

	} /// �̸� ����ϴ� Hm�κп��� MAX_NUM_CLUSTER�� numpath�� ���� �ʾƼ� ������ ���� // idx�� num_path���� Ŭ �� delay���� �����Ⱚ�� ���� �̷��� ���� 0���� �ʱ�ȭ�ϵ����ߴ�

	for (int i = 0; i < num_path; i++)
	{
		delay[i] = -r_tau * _DS * log(randnum.u());
	}

	// cout << "delay = " << delay[0] << endl;

	for (int i = 0; i < num_path; i++)
	{
		Real temp;
		if (i == (num_path - 1))
		{
			break;
		}

		for (int j = i + 1; j < num_path; j++)
		{
			temp = delay[i];

			if (delay[j] < temp)
			{
				temp = delay[j];
				delay[j] = delay[i];
				delay[i] = temp;
			}
		}
	} // sort delay, small -> large

	Real temp_delay = delay[0]; // minimum

	for (int i = 0; i < num_path; i++)
	{

		delay[i] = delay[i] - temp_delay;
	}

	// LOS CASE -> scaled delay for delay spread, not to be used in cluster power generation (M.2412 Small scale generation Step 5)
	for (int i = 0; i < num_path; i++)
	{

		delay_LOS[i] = delay[i] / D;
		delay[i] = delay[i];
	}
}

void CHANNEL::Set_POWER()
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter -- power
	///
	/////////////////////////////////////////////////

	// power = new Real[MAX_NUM_CLUSTERS];
	// power_LOS = new Real[MAX_NUM_CLUSTERS];
	// power_NLOS = new Real[MAX_NUM_CLUSTERS];   /// LOS, NLOS separate for pre-calc

	Real sum_power = 0.;
	Real sum_power_LOS = 0.;
	Real sum_power_NLOS = 0.;

	// Real k_R = pow(10, K_factor / 10);
	Real k_R = K_linear; /////

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		power[i] = 0.;
		power_LOS[i] = 0.;
		power_NLOS[i] = 0.;

	} //// initialize

	for (int i = 0; i < num_path; i++)
	{
		power[i]        = exp(-delay[i] * (r_tau - 1.) / (r_tau * _DS)) * pow(10., -1 * normal(0, cluster_shadowing) / 10.);

		power_LOS[i]    = power[i]; // delay LOS -> cluster power generation
		power_NLOS[i]   = power[i];

		sum_power      += power[i];
		sum_power_LOS  += power_LOS[i];
		sum_power_NLOS += power_NLOS[i];
	}

	// cout << "power = " << power[0] << endl;

	if (LOS == 1)
	{
		for (int i = 0; i < num_path; i++)
		{
			power[i] = power[i] / (sum_power * (k_R + 1));
		}
		power[0] = power[0] + k_R / (k_R + 1);
	}
	else if (LOS == 0)
	{
		for (int i = 0; i < num_path; i++)
		{
			power[i] = power[i] / sum_power;
		}
	}

	///// ================================================================== ////
	////  LOS, NLOS case for precalculation                                  ////
	////  ================================================================== ////
	// for (int i = 0; i < num_path; i++)
	//{
	//	power_LOS[i] = power_LOS[i] / (sum_power_LOS * (k_R + 1));
	// }
	// power_LOS[0] = power_LOS[0] + k_R / (k_R + 1);

	// for (int i = 0; i < num_path; i++)
	//{
	//	power_NLOS[i] = power_NLOS[i] / sum_power_NLOS;
	// }
	////  ================================================================== ////
}

void CHANNEL::Find_Strong2Clusters()
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter --  find the 2 strongest clusters
	///
	/////////////////////////////////////////////////

	strongest_power = power[0];
	strongest_power2 = power[1];
	strongest_power_idx = 0;
	strongest_power_idx2 = 1;

	if (power[0] < power[1])
	{
		strongest_power = power[1];
		strongest_power2 = power[0];
		strongest_power_idx = 1;
		strongest_power_idx2 = 0;
	}

	for (int i = 2; i < num_path; i++)
	{
		if (power[i] > strongest_power2)
		{
			if (power[i] > strongest_power)
			{
				strongest_power2 = strongest_power;
				strongest_power_idx2 = strongest_power_idx;
				strongest_power = power[i];
				strongest_power_idx = i;
			}
			else
			{
				strongest_power2 = power[i];
				strongest_power_idx2 = i;
			}
		}
	}
}

void CHANNEL::Set_AOAAOD(int _bs_idx, int _ue_idx)
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter -- AOA, AOD
	///
	/////////////////////////////////////////////////

	Real K_factor_dB = K_dB;
	// Real K_factor_dB = 10.;

	// Real ASD = pow(10, normal(mu_ASD, sigma_ASD));
	// Real ASA = pow(10, normal(mu_ASA, sigma_ASA));
	Real C = 0;

	//////////////////////////////////////////////////////// IMT-2020 /////////////////////////////////////////////////////////////////////////////////////////////////
	////// IMT 2020 EVAL Page 51 ~
	////////////////////////////////////////
	/////////////////////////////////////////////////////////// InH
	if (TYPE == 11)
	{
		if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			{
				if (Propagation == 1) ///// LOS  Laplacian    num_clusters = 15
				{
					C = 1.434 * (0.9275 + 0.0439 * K_factor_dB - 0.0071 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS   Laplacian    num_clusters = 19
				{
					C = 1.501;
				}
			}
			// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
			else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
			{
				if (Propagation == 1) ///// LOS  Gaussian    num_clusters = 15
				{
					C = 1.211 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS   Gaussian   num_clusters = 19
				{
					C = 1.273;
				}
			}
		}
		else if (Channel_Model_Type == 1) // Channel Model B
		{
			if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
			// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
			{
				if (Propagation == 1) ///// LOS  Gaussian   num_clusters = 15
				{
					C = 1.211 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS  Gaussian    num_clusters = 19
				{
					C = 1.273;
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////// Dense_Urban
	else if (TYPE == 12)
	{
		int lc_Urban_env_Config_Type = Configuration_Type;
		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		{
			if ( self_bs_idx < simple_num_BS )
				lc_Urban_env_Config_Type = 1;
			else
				lc_Urban_env_Config_Type = 3;
		}		
		////////////////////////////////////////////// UMi
		if (lc_Urban_env_Config_Type == 3)
		{
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian  num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 19
					{
						C = 1.273;
					}
					else if (Propagation == 2) ////OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian    num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 19
					{
						C = 1.273;
					}
					else if (Propagation == 2) ////  OUT2IN  Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
			}
			else if (Channel_Model_Type == 1) // Channel Model B
			{
				if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
				// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian   num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS  Gaussian   num_clusters = 19
					{
						C = 1.273;
					}
					else if (Propagation == 2) //// OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
			}
		}
		////////////////////////////////////////////// UMa
		else
		{
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian  num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 20
					{
						C = 1.289;
					}
					else if (Propagation == 2) ////OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian    num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 20
					{
						C = 1.289;
					}
					else if (Propagation == 2) ////  OUT2IN  Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
			}
			else if (Channel_Model_Type == 1) // Channel Model B
			{
				// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
				if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian   num_clusters = 12
					{
						C = 1.146 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS  Gaussian   num_clusters = 20
					{
						C = 1.289;
					}
					else if (Propagation == 2) //// OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.146;
					}
				}
			}
		}
	}
	///////////////////////////////////////////////////// RMa
	else if (TYPE == 13)
	{
		if (Channel_Model_Type == 0) /// Channel Model A
		{
			//if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			if (carrier_freq >= 500000000)
			{
				if (Propagation == 1) ///// LOS Laplacian?? Gaussian??   num_clusters = 11
				{
					C = 1.123 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS      num_clusters = 10
				{
					C = 1.090;
				}
				else if (Propagation == 2) ////OUT2IN    num_clusters = 10
				{
					C = 1.090;
				}
			}
		}
		else if (Channel_Model_Type == 1) // Channel Model B
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 7000000000) //// 0.5GHz <= fc <= 7GHz
			{
				if (Propagation == 1) ///// LOS   num_clusters = 11
				{
					C = 1.123 * (1.1035 - 0.028 * K_factor_dB - 0.002 * K_factor_dB * K_factor_dB + 0.0001 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS    num_clusters = 10
				{
					C = 1.090;
				}
				else if (Propagation == 2) //// OUT2IN   num_clusters = 10
				{
					C = 1.090;
				}
			}
		}
	}

	// cout << "C = " << C << endl;

	// AOA = new Real[MAX_NUM_CLUSTERS];
	// AOD = new Real[MAX_NUM_CLUSTERS];

	Real first_AOD;
	Real first_AOA;

	Real X1_AOD;
	Real Y1_AOD;
	Real X1_AOA;
	Real Y1_AOA;

	for (int i = 0; i < num_path; i++)
	{

		Real R = randnum.u();
		Real X;
		if (R > 0.5)
		{
			X = 1;
		}
		else
		{
			X = -1;
		}

		// Real Y_AOA = sigma_AOA / 7. * randnum.n();
		// Real Y_AOD = sigma_AOD / 7. * randnum.n();

		Real Y_AOA = normal(0, ASA / 7.);
		Real Y_AOD = normal(0, ASD / 7.);

		// Real Y_AOA = normal(0, (ASA / 7. , 2.));
		// Real Y_AOD = normal(0, (ASD / 7. , 2.));

		if (i == 0)
		{
			X1_AOD = X;
			Y1_AOD = Y_AOD;

			X1_AOA = X;
			Y1_AOA = Y_AOA;
		}

		if (power[i] == 0.)
		{
			AOD[i] = 0.;
			AOA[i] = 0.;
		}
		else
		{
			if ((TYPE == 11) &&  ((Channel_Model_Type == 0) || (Channel_Model_Type == 2)) && (carrier_freq >= 500000000 && carrier_freq <= 6000000000)) //// 0.5GHz <= fc <= 6GHz)   //// Laplacian (InH_A && F<=6GHz)
			{
				AOD[i] = -1 * ASD * log(power[i] / strongest_power) / C;
				AOA[i] = -1 * ASA * log(power[i] / strongest_power) / C;
			}
			else ///// Gaussian
			{
				AOD[i] = 2. * (ASD / 1.4) * sqrt(-1 * log(power[i] / strongest_power)) / C;
				AOA[i] = 2. * (ASA / 1.4) * sqrt(-1 * log(power[i] / strongest_power)) / C;
			}
		}

		if (i == 0)
		{
			first_AOD = AOD[0];
			first_AOA = AOA[0];
		}

		if (channel[_bs_idx][_ue_idx].self_ms.y > channel[_bs_idx][_ue_idx].self_bs.y)
		{
			if (channel[_bs_idx][_ue_idx].self_ms.x > channel[_bs_idx][_ue_idx].self_bs.x)
			{
				// cout << "case1" << endl;
				LOS_AOD_GCS = atan((channel[_bs_idx][_ue_idx].self_ms.y - channel[_bs_idx][_ue_idx].self_bs.y) / (channel[_bs_idx][_ue_idx].self_ms.x - channel[_bs_idx][_ue_idx].self_bs.x)) * (180. / pi);
				LOS_AOA_GCS = (180. + LOS_AOD_GCS) - 360.;
			}
			else
			{
				// cout << "case2" << endl;
				LOS_AOD_GCS = 180. - atan((channel[_bs_idx][_ue_idx].self_ms.y - channel[_bs_idx][_ue_idx].self_bs.y) / (channel[_bs_idx][_ue_idx].self_bs.x - channel[_bs_idx][_ue_idx].self_ms.x)) * (180. / pi);
				LOS_AOA_GCS = -1. * atan((channel[_bs_idx][_ue_idx].self_ms.y - channel[_bs_idx][_ue_idx].self_bs.y) / (channel[_bs_idx][_ue_idx].self_bs.x - channel[_bs_idx][_ue_idx].self_ms.x)) * (180. / pi);
			}
		}
		else ///// (ms[_ue_idx].loc.y < bs[_bs_idx].loc.y)
		{
			if (channel[_bs_idx][_ue_idx].self_ms.x > channel[_bs_idx][_ue_idx].self_bs.x)
			{
				// cout << "case3" << endl;
				LOS_AOD_GCS = -1. * atan((channel[_bs_idx][_ue_idx].self_bs.y - channel[_bs_idx][_ue_idx].self_ms.y) / (channel[_bs_idx][_ue_idx].self_ms.x - channel[_bs_idx][_ue_idx].self_bs.x)) * (180. / pi);
				LOS_AOA_GCS = 180. - atan((channel[_bs_idx][_ue_idx].self_bs.y - channel[_bs_idx][_ue_idx].self_ms.y) / (channel[_bs_idx][_ue_idx].self_ms.x - channel[_bs_idx][_ue_idx].self_bs.x)) * (180. / pi);
			}
			else
			{
				// cout << "case4" << endl;

				LOS_AOA_GCS = atan((channel[_bs_idx][_ue_idx].self_bs.y - channel[_bs_idx][_ue_idx].self_ms.y) / (channel[_bs_idx][_ue_idx].self_bs.x - channel[_bs_idx][_ue_idx].self_ms.x)) * (180. / pi);

				LOS_AOD_GCS = LOS_AOA_GCS - 180.;
			}
		}
		/*
		cout << "channel[_bs_idx][_ue_idx].self_ms.x = " << channel[_bs_idx][_ue_idx].self_ms.x << endl;
		cout << "channel[_bs_idx][_ue_idx].self_ms.y = " << channel[_bs_idx][_ue_idx].self_ms.y << endl;
		cout << "channel[_bs_idx][_ue_idx].self_bs.x = " << channel[_bs_idx][_ue_idx].self_bs.x << endl;
		cout << "channel[_bs_idx][_ue_idx].self_bs.y = " << channel[_bs_idx][_ue_idx].self_bs.y << endl;

		cout << "LOS_AOD_GCS = " << LOS_AOD_GCS << endl;
		cout << "LOS_AOA_GCS = " << LOS_AOA_GCS << endl;
		getchar();
		*/

		if (LOS_AOD_GCS > 180. || LOS_AOD_GCS < -180.)
		{
			cout << "LOS_AOD_GCS Error!  channel.cpp" << endl;
			getchar();
		}

		if (LOS_AOA_GCS > 180. || LOS_AOA_GCS < -180.)
		{
			cout << "LOS_AOA_GCS Error!  channel.cpp" << endl;
			getchar();
		}

		if (Propagation == 1) // LOS
		{
			/*
			LOCATION vector;

			// bs_idx, ms_idx -> sector_idx
			Real sector_idx;

			if (self_ms_idx % 30 < 10)
			{
				sector_idx = 1;
			}
			else if (self_ms_idx % 30 < 20)
			{
				sector_idx = 2;
			}
			else if (self_ms_idx % 30 < 30)
			{
				sector_idx = 3;
			}

			vector.x = bs_location.x + cos((30. + 120. * (sector_idx - 1.))*pi / 180.);
			vector.y = bs_location.y + sin((30. + 120. * (sector_idx - 1.))*pi / 180.);

			Real angle = links[self_ms_idx].Get_angle(bs_location, ms_location, vector);
			*/

			AOD[i] = X * AOD[i] + Y_AOD - (X1_AOD * first_AOD + Y1_AOD - LOS_AOD_GCS);
			AOA[i] = X * AOA[i] + Y_AOA - (X1_AOA * first_AOA + Y1_AOA - LOS_AOA_GCS);
		}
		else // NLOS, O2I
		{
			AOD[i] = X * AOD[i] + Y_AOD + LOS_AOD_GCS;
			AOA[i] = X * AOA[i] + Y_AOA + LOS_AOA_GCS;
		}
	}
}

void CHANNEL::Set_ZOAZOD(int _bs_idx, int _ue_idx)
{

	/////////////////////////////////////////////////
	///
	/// Small Scale parameter -- ZOA, ZOD
	///
	/////////////////////////////////////////////////

	// Real K_factor_dB = 10;
	Real K_factor_dB = K_dB;
	// Real ZSD = pow(10, normal(mu_ZSD, sigma_ZSD));
	// Real ZSA = pow(10, normal(mu_ZSA, sigma_ZSA));
	Real C = 1;

	//////////////////////////////////////////////////////// IMT-2020 /////////////////////////////////////////////////////////////////////////////////////////////////
	////// IMT 2020 EVAL Page 51 ~
	if (TYPE == 11) // InH
	{
		if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
		{

			if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			{
				if (Propagation == 1) ///// LOS  Laplacian    num_clusters = 15
				{
					C = 1.1088 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS   Laplacian    num_clusters = 19
				{
					C = 1.1764;
				}
			}
			else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
			// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
			{
				if (Propagation == 1) ///// LOS  Gaussian    num_clusters = 15
				{
					C = 1.1088 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS   Gaussian   num_clusters = 19
				{
					C = 1.1764;
				}
			}
		}
		else if (Channel_Model_Type == 1) // Channel Model B
		{
			if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
			// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
			{
				if (Propagation == 1) ///// LOS  Gaussian   num_clusters = 15
				{
					C = 1.1088 * (1.3086 + 0.0339 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS  Gaussian    num_clusters = 19
				{
					C = 1.184;
				}
			}
		}
	}
	else if (TYPE == 12) // Dense_Urban
	{
		int lc_Urban_env_Config_Type = Configuration_Type;
		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		{
			if ( self_bs_idx < simple_num_BS )
				lc_Urban_env_Config_Type = 1;
			else
				lc_Urban_env_Config_Type = 3;
		}
		///////////////////////////////////////////// UMi
		if (lc_Urban_env_Config_Type == 3)
		{
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					if (Propagation == 1) ///// LOS    num_clusters = 12
					{
						C = 1.04 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS      num_clusters = 19
					{
						C = 1.1764;
					}
					else if (Propagation == 2) ////OUT2IN   num_clusters = 12
					{
						C = 1.04;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS      num_clusters = 12
					{
						C = 1.04 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS     num_clusters = 19
					{
						C = 1.1764;
					}
					else if (Propagation == 2) ////  OUT2IN    num_clusters = 12
					{
						C = 1.04;
					}
				}
			}
			else if (Channel_Model_Type == 1) // Channel Model B
			{
				if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
				// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS     num_clusters = 12
					{
						C = 1.104 * (1.3086 + 0.0339 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS    num_clusters = 19
					{
						C = 1.184;
					}
					else if (Propagation == 2) //// OUT2IN    num_clusters = 12
					{
						C = 1.104;
					}
				}
			}
		}
		///////////////////////////////////////////// UMa
		else
		{
			if (Channel_Model_Type == 0 || Channel_Model_Type == 2) /// Channel Model A
			{
				if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian  num_clusters = 12
					{
						C = 1.04 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 20
					{
						C = 1.1918;
					}
					else if (Propagation == 2) ////OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.04;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				// else if (carrier_freq > 6000000000 && carrier_freq <= 100000000000) /// 6GHz < fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian    num_clusters = 12
					{
						C = 1.04 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS   Gaussian    num_clusters = 20
					{
						C = 1.1918;
					}
					else if (Propagation == 2) ////  OUT2IN  Gaussian   num_clusters = 12
					{
						C = 1.04;
					}
				}
			}
			else if (Channel_Model_Type == 1) // Channel Model B
			{
				if (carrier_freq >= 500000000) //// 0.5GHz <= fc <= 100GHz
				// if (carrier_freq >= 500000000 && carrier_freq <= 100000000000) //// 0.5GHz <= fc <= 100GHz
				{
					if (Propagation == 1) ///// LOS  Gaussian   num_clusters = 12
					{
						C = 1.104 * (1.3086 + 0.0339 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
					}
					else if (Propagation == 0) /////////// NLOS  Gaussian   num_clusters = 20
					{
						C = 1.178;
					}
					else if (Propagation == 2) //// OUT2IN Gaussian   num_clusters = 12
					{
						C = 1.104;
					}
				}
			}
		}
	}
	else if (TYPE == 13) // RMa
	{
		if (Channel_Model_Type == 0) /// Channel Model A
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 6000000000) //// 0.5GHz <= fc <= 6GHz
			{
				if (Propagation == 1) ///// LOS Laplacian?? Gaussian??   num_clusters = 11
				{
					C = 1.013 * (1.35 + 0.0202 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS      num_clusters = 10
				{
					C = 0.9854;
				}
				else if (Propagation == 2) ////OUT2IN    num_clusters = 10
				{
					C = 0.9854;
				}
			}
		}
		else if (Channel_Model_Type == 1) // Channel Model B
		{
			if (carrier_freq >= 500000000 && carrier_freq <= 7000000000) //// 0.5GHz <= fc <= 7GHz
			{
				if (Propagation == 1) ///// LOS   num_clusters = 11
				{
					C = 1.031 * (1.3086 + 0.0339 * K_factor_dB - 0.0077 * K_factor_dB * K_factor_dB + 0.0002 * K_factor_dB * K_factor_dB * K_factor_dB);
				}
				else if (Propagation == 0) /////////// NLOS    num_clusters = 10
				{
					C = 0.957;
				}
				else if (Propagation == 2) //// OUT2IN   num_clusters = 10
				{
					C = 0.957;
				}
			}
		}
	}

	// cout << "C_zoazod = " << C << endl;

	// ZOD = new Real[MAX_NUM_CLUSTERS];
	// ZOA = new Real[MAX_NUM_CLUSTERS];

	Real first_ZOD;
	Real first_ZOA;

	Real X1_ZOD;
	Real Y1_ZOD;
	Real X1_ZOA;
	Real Y1_ZOA;

	for (int i = 0; i < num_path; i++)
	{

		Real R = randnum.u();
		Real X;
		if (R > 0.5)
		{
			X = 1;
		}
		else
		{
			X = -1;
		}

		// Real Y_ZOA = normal(0, (ZSA / 7., 2.)); //
		// Real Y_ZOD = normal(0, (ZSD / 7., 2.));

		Real Y_ZOA = normal(0, ZSA / 7.); //
		Real Y_ZOD = normal(0, ZSD / 7.);

		if (i == 0)
		{
			X1_ZOD = X;
			Y1_ZOD = Y_ZOD;

			X1_ZOA = X;
			Y1_ZOA = Y_ZOA;
		}

		if (power[i] == 0.)
		{
			ZOD[i] = 0.;
			ZOA[i] = 0.;
		}
		else
		{
			ZOD[i] = -1 * ZSD * log(power[i] / strongest_power) / C;
			ZOA[i] = -1 * ZSA * log(power[i] / strongest_power) / C;
		}
		if (i == 0)
		{
			first_ZOD = ZOD[0];
			first_ZOA = ZOA[0];
		}

		LOS_ZOA_GCS = 90. - (atan((_bs_height - ms[_ue_idx].MS_HEIGHT_FINAL) / distance) * (180. / pi));
		LOS_ZOD_GCS = 180. - LOS_ZOA_GCS;

		/*
		cout << "distance = " << distance << endl;
		cout << "(_bs_height - ms[_ue_idx].MS_HEIGHT_FINAL) = " << (_bs_height - ms[_ue_idx].MS_HEIGHT_FINAL) << endl;
		cout << "LOS_ZOD_GCS = " << LOS_ZOD_GCS << endl;
		cout << "LOS_ZOA_GCS = " << LOS_ZOA_GCS << endl;
		getchar();
		*/

		if (LOS_ZOD_GCS > 180. || LOS_ZOD_GCS < 0.)
		{
			cout << "LOS_ZOD_GCS Error!  channel.cpp" << endl;
			getchar();
		}

		if (LOS_ZOA_GCS > 180. || LOS_ZOA_GCS < 0.)
		{
			cout << "LOS_ZOD_GCS Error!  channel.cpp" << endl;
			getchar();
		}

		if (Propagation == 1) // LOS
		{
			/*
			LOCATION vector;

			// bs ms -> angle
			Real sector_idx;

			if (self_ms_idx % 30 < 10)
			{
				sector_idx = 1;
			}
			else if (self_ms_idx % 30 < 20)
			{
				sector_idx = 2;
			}
			else if (self_ms_idx % 30 < 30)
			{
				sector_idx = 3;
			}
			vector.x = bs_location.x + cos((30. + 120. * (sector_idx - 1.))*pi / 180.);
			vector.y = bs_location.y + sin((30. + 120. * (sector_idx - 1.))*pi / 180.);

			Real angle = links[self_ms_idx].Get_angle(bs_location, ms_location, vector);
			*/

			ZOA[i] = (X * ZOA[i] + Y_ZOA) - (X1_ZOA * first_ZOA + Y1_ZOA - LOS_ZOA_GCS);
			ZOD[i] = (X * ZOD[i] + Y_ZOD) - (X1_ZOD * first_ZOD + Y1_ZOD - LOS_ZOD_GCS);
		}
		else // NLOS, O2I
		{
			if (Propagation == 2)
			{
				ZOA[i] = X * ZOA[i] + Y_ZOA + 90.;
				ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
			}

			else if (Propagation == 0)
			{
				ZOA[i] = X * ZOA[i] + Y_ZOA + LOS_ZOA_GCS;
				ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
			}
		}
		/*
		if (_ue_idx == 4)
		{
			cout << "LOS = " << LOS << endl;
			cout << "mu_offset_ZOD = " << mu_offset_ZOD << endl;
			cout << "X = " << X << endl;
			cout << "Y_ZOD = " << Y_ZOD << endl;
			cout << "X1_ZOD = " << X1_ZOD << endl;
			cout << "first_ZOD = " << first_ZOD << endl;
			cout << "Y1_ZOD = " << Y1_ZOD << endl;
			cout << "LOS_ZOD = " << LOS_ZOD << endl;
			cout << "ZOD[0] = " << ZOD[0] << endl;
			getchar();
		}
		*/
	}
}

void CHANNEL::Set_InitialPhase()
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter -- initial phase
	///
	/////////////////////////////////////////////////

	/*
	random_phase_vv = new Real *[MAX_NUM_CLUSTERS];
	random_phase_vh = new Real *[MAX_NUM_CLUSTERS];
	random_phase_hv = new Real *[MAX_NUM_CLUSTERS];
	random_phase_hh = new Real *[MAX_NUM_CLUSTERS];

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		random_phase_vv[i] = new Real[MAX_NUM_RAYS];
		random_phase_vh[i] = new Real[MAX_NUM_RAYS];
		random_phase_hv[i] = new Real[MAX_NUM_RAYS];
		random_phase_hh[i] = new Real[MAX_NUM_RAYS];
	}
	*/

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		for (int j = 0; j < MAX_NUM_RAYS; j++)
		{

			random_phase_vv[i][j] = 360. * randnum.u() - 180.;
			random_phase_vh[i][j] = 360. * randnum.u() - 180.;
			random_phase_hv[i][j] = 360. * randnum.u() - 180.;
			random_phase_hh[i][j] = 360. * randnum.u() - 180.;
		}
	}

	random_phase_vv_LOS = 360. * randnum.u() - 180.;
	random_phase_hh_LOS = 360. * randnum.u() - 180.;
}

void CHANNEL::Sampling_DelaySpread()
{
	Real sampling_rate   = subcarrier_spacing*fft_size;
	Real sampling_period = 1/sampling_rate;

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++ )
	{
		sampled_delay[i] = (int)round(delay[i]/sampling_period);
	}
}

void CHANNEL::Set_SUBCLUSTER()
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////
	////// offset_angle
	///////////////////////////////////

	/*
	offset_angle = new Real *[MAX_NUM_CLUSTERS];

	for (int path_idx = 0; path_idx < MAX_NUM_CLUSTERS; path_idx++)
	{
		offset_angle[path_idx] = new Real[MAX_NUM_RAYS];

	}
	*/

	for (int path_idx = 0; path_idx < MAX_NUM_CLUSTERS; path_idx++)
	{
		NUM_RAY_per_ClusterNUM[path_idx] = 0;

		for (int ray_idx = 0; ray_idx < MAX_NUM_RAYS; ray_idx++)
		{
			offset_angle[path_idx][ray_idx] = 0.;
			offset_angle_rand_coupling[path_idx][ray_idx] = 0.;
		}
	}

	NUM_PATH_for_channelcoeff = 0; //// find SSP per time, so initialize

	//////////////////////////////////////////////
	/////
	///// Set_subcluster
	////  strong 2 cluster -> divide in 3 cluster and add offset, -> In my code, not divide -> add 4 cluster... the number of ray is different
	////
	/////////////////////////////////////////////

	// modified by jhnoh (multiplying 1e-9 for making nanosec)
	if (cluster_DS < 0) ///// cluster_DS == N/A
	{
		delay[num_path] = delay[strongest_power_idx]      + 1.28 * 3.91 * 1e-9;
		delay[num_path + 1] = delay[strongest_power_idx]  + 2.56 * 3.91 * 1e-9;
		delay[num_path + 2] = delay[strongest_power_idx2] + 1.28 * 3.91 * 1e-9;
		delay[num_path + 3] = delay[strongest_power_idx2] + 2.56 * 3.91 * 1e-9;
	}
	else
	{
		delay[num_path] = delay[strongest_power_idx]      + 1.28 * cluster_DS * 1e-9;
		delay[num_path + 1] = delay[strongest_power_idx]  + 2.56 * cluster_DS * 1e-9;
		delay[num_path + 2] = delay[strongest_power_idx2] + 1.28 * cluster_DS * 1e-9;
		delay[num_path + 3] = delay[strongest_power_idx2] + 2.56 * cluster_DS * 1e-9;
	}

	power[num_path]             = power[strongest_power_idx ] * 6. / 20.;
	power[num_path + 1]         = power[strongest_power_idx ] * 4. / 20.;
	power[num_path + 2]         = power[strongest_power_idx2] * 6. / 20.;
	power[num_path + 3]         = power[strongest_power_idx2] * 4. / 20.;
	power[strongest_power_idx]  = power[strongest_power_idx ] * 10. / 20.;
	power[strongest_power_idx2] = power[strongest_power_idx2] * 10. / 20.;

	AOA[num_path] = AOA[strongest_power_idx];
	AOA[num_path + 1] = AOA[strongest_power_idx];
	AOA[num_path + 2] = AOA[strongest_power_idx2];
	AOA[num_path + 3] = AOA[strongest_power_idx2];

	AOD[num_path] = AOD[strongest_power_idx];
	AOD[num_path + 1] = AOD[strongest_power_idx];
	AOD[num_path + 2] = AOD[strongest_power_idx2];
	AOD[num_path + 3] = AOD[strongest_power_idx2];

	ZOA[num_path] = ZOA[strongest_power_idx];
	ZOA[num_path + 1] = ZOA[strongest_power_idx];
	ZOA[num_path + 2] = ZOA[strongest_power_idx2];
	ZOA[num_path + 3] = ZOA[strongest_power_idx2];

	ZOD[num_path] = ZOD[strongest_power_idx];
	ZOD[num_path + 1] = ZOD[strongest_power_idx];
	ZOD[num_path + 2] = ZOD[strongest_power_idx2];
	ZOD[num_path + 3] = ZOD[strongest_power_idx2];

	for (int path_idx = 0; path_idx < num_path + 4; path_idx++)
	{
		if (path_idx == strongest_power_idx || path_idx == strongest_power_idx2) /// Subcluster No.1 (1,2,3,4,5,6,7,8,19,20 rays)
		{
			// num_ray = 10;
			NUM_RAY_per_ClusterNUM[path_idx] = 10;

			offset_angle[path_idx][0] = 0.0447;
			offset_angle[path_idx][1] = -0.0447;
			offset_angle[path_idx][2] = 0.1413;
			offset_angle[path_idx][3] = -0.1413;
			offset_angle[path_idx][4] = 0.2492;
			offset_angle[path_idx][5] = -0.2492;
			offset_angle[path_idx][6] = 0.3715;
			offset_angle[path_idx][7] = -0.3715;
			offset_angle[path_idx][8] = 2.1551;
			offset_angle[path_idx][9] = -2.1551;

			vector<int> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
			random_shuffle(v.begin(), v.end());

			offset_angle_rand_coupling[path_idx][0] = offset_angle[path_idx][v.at(0)]; //// random coupling AOD-AOA
			offset_angle_rand_coupling[path_idx][1] = offset_angle[path_idx][v.at(1)];
			offset_angle_rand_coupling[path_idx][2] = offset_angle[path_idx][v.at(2)];
			offset_angle_rand_coupling[path_idx][3] = offset_angle[path_idx][v.at(3)];
			offset_angle_rand_coupling[path_idx][4] = offset_angle[path_idx][v.at(4)];
			offset_angle_rand_coupling[path_idx][5] = offset_angle[path_idx][v.at(5)];
			offset_angle_rand_coupling[path_idx][6] = offset_angle[path_idx][v.at(6)];
			offset_angle_rand_coupling[path_idx][7] = offset_angle[path_idx][v.at(7)];
			offset_angle_rand_coupling[path_idx][8] = offset_angle[path_idx][v.at(8)];
			offset_angle_rand_coupling[path_idx][9] = offset_angle[path_idx][v.at(9)];
		}

		else if (path_idx == num_path || path_idx == num_path + 2) /// Subcluster No.2 (9,10,11,12,17,18 rays)
		{
			// num_ray = 6;
			NUM_RAY_per_ClusterNUM[path_idx] = 6;
			offset_angle[path_idx][0] = 0.5129;
			offset_angle[path_idx][1] = -0.5129;
			offset_angle[path_idx][2] = 0.6797;
			offset_angle[path_idx][3] = -0.6797;
			offset_angle[path_idx][4] = 1.5195;
			offset_angle[path_idx][5] = -1.5195;

			vector<int> v = {0, 1, 2, 3, 4, 5};
			random_shuffle(v.begin(), v.end());

			offset_angle_rand_coupling[path_idx][0] = offset_angle[path_idx][v.at(0)];
			offset_angle_rand_coupling[path_idx][1] = offset_angle[path_idx][v.at(1)];
			offset_angle_rand_coupling[path_idx][2] = offset_angle[path_idx][v.at(2)];
			offset_angle_rand_coupling[path_idx][3] = offset_angle[path_idx][v.at(3)];
			offset_angle_rand_coupling[path_idx][4] = offset_angle[path_idx][v.at(4)];
			offset_angle_rand_coupling[path_idx][5] = offset_angle[path_idx][v.at(5)];
		}

		else if (path_idx == num_path + 1 || path_idx == num_path + 3) /// Subcluster No.3 (13,14,15,16 rays)
		{
			// num_ray = 4;
			NUM_RAY_per_ClusterNUM[path_idx] = 4;
			offset_angle[path_idx][0] = 0.8844;
			offset_angle[path_idx][1] = -0.8844;
			offset_angle[path_idx][2] = 1.1481;
			offset_angle[path_idx][3] = -1.1481;

			vector<int> v = {0, 1, 2, 3};
			random_shuffle(v.begin(), v.end());

			offset_angle_rand_coupling[path_idx][0] = offset_angle[path_idx][v.at(0)];
			offset_angle_rand_coupling[path_idx][1] = offset_angle[path_idx][v.at(1)];
			offset_angle_rand_coupling[path_idx][2] = offset_angle[path_idx][v.at(2)];
			offset_angle_rand_coupling[path_idx][3] = offset_angle[path_idx][v.at(3)];
		}

		else /// Not 2 strong cluster
		{
			// num_ray = 20;
			NUM_RAY_per_ClusterNUM[path_idx] = 20;
			offset_angle[path_idx][0] = 0.0447;
			offset_angle[path_idx][1] = -0.0447;
			offset_angle[path_idx][2] = 0.1413;
			offset_angle[path_idx][3] = -0.1413;
			offset_angle[path_idx][4] = 0.2492;
			offset_angle[path_idx][5] = -0.2492;
			offset_angle[path_idx][6] = 0.3715;
			offset_angle[path_idx][7] = -0.3715;
			offset_angle[path_idx][8] = 0.5129;
			offset_angle[path_idx][9] = -0.5129;
			offset_angle[path_idx][10] = 0.6797;
			offset_angle[path_idx][11] = -0.6797;
			offset_angle[path_idx][12] = 0.8844;
			offset_angle[path_idx][13] = -0.8844;
			offset_angle[path_idx][14] = 1.1481;
			offset_angle[path_idx][15] = -1.1481;
			offset_angle[path_idx][16] = 1.5195;
			offset_angle[path_idx][17] = -1.5195;
			offset_angle[path_idx][18] = 2.1551;
			offset_angle[path_idx][19] = -2.1551;

			vector<int> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
			random_shuffle(v.begin(), v.end());

			offset_angle_rand_coupling[path_idx][0] = offset_angle[path_idx][v.at(0)];
			offset_angle_rand_coupling[path_idx][1] = offset_angle[path_idx][v.at(1)];
			offset_angle_rand_coupling[path_idx][2] = offset_angle[path_idx][v.at(2)];
			offset_angle_rand_coupling[path_idx][3] = offset_angle[path_idx][v.at(3)];
			offset_angle_rand_coupling[path_idx][4] = offset_angle[path_idx][v.at(4)];
			offset_angle_rand_coupling[path_idx][5] = offset_angle[path_idx][v.at(5)];
			offset_angle_rand_coupling[path_idx][6] = offset_angle[path_idx][v.at(6)];
			offset_angle_rand_coupling[path_idx][7] = offset_angle[path_idx][v.at(7)];
			offset_angle_rand_coupling[path_idx][8] = offset_angle[path_idx][v.at(8)];
			offset_angle_rand_coupling[path_idx][9] = offset_angle[path_idx][v.at(9)];
			offset_angle_rand_coupling[path_idx][10] = offset_angle[path_idx][v.at(10)];
			offset_angle_rand_coupling[path_idx][11] = offset_angle[path_idx][v.at(11)];
			offset_angle_rand_coupling[path_idx][12] = offset_angle[path_idx][v.at(12)];
			offset_angle_rand_coupling[path_idx][13] = offset_angle[path_idx][v.at(13)];
			offset_angle_rand_coupling[path_idx][14] = offset_angle[path_idx][v.at(14)];
			offset_angle_rand_coupling[path_idx][15] = offset_angle[path_idx][v.at(15)];
			offset_angle_rand_coupling[path_idx][16] = offset_angle[path_idx][v.at(16)];
			offset_angle_rand_coupling[path_idx][17] = offset_angle[path_idx][v.at(17)];
			offset_angle_rand_coupling[path_idx][18] = offset_angle[path_idx][v.at(18)];
			offset_angle_rand_coupling[path_idx][19] = offset_angle[path_idx][v.at(19)];
		}
	}

	NUM_PATH_for_channelcoeff = num_path + 4;
}

void CHANNEL::Set_RMS_delay_spread()
{
	if (LOS == 0) // NLOS
	{
		// delay_power   = new Real[MAX_NUM_CLUSTERS];
		// power_ray_sum = new Real[MAX_NUM_CLUSTERS];

		for (int i = 0; i < num_path; i++)
		{
			// delay_power[i] = num_ray * delay[i] * power[i];
			// power_ray_sum[i] = num_ray * power[i];

			delay_power[i] = delay[i] * power[i];
			power_ray_sum[i] = power[i];
		}

		Real avg_delay;
		Real delay_power_sum = 0;
		Real power_sum = 0;

		for (int i = 0; i < num_path; i++)
		{
			delay_power_sum += delay_power[i];
			power_sum += power_ray_sum[i];
		}

		avg_delay = delay_power_sum / power_sum;
		for (int i = 0; i < num_path; i++)
		{
			RMS_delay_power[i] = (delay[i] - avg_delay) * (delay[i] - avg_delay) * power[i];
		}

		// Real RMS_delay_spread = 0;
		Real RMS_delay_power_sum = 0;

		for (int i = 0; i < num_path; i++)
		{
			RMS_delay_power_sum += RMS_delay_power[i];
		}

		RMS_delay_spread = sqrt(RMS_delay_power_sum / power_sum);
	}

	else if (LOS == 1) // LOS
	{
		// delay_power   = new Real[MAX_NUM_CLUSTERS];
		// power_ray_sum = new Real[MAX_NUM_CLUSTERS];

		for (int i = 0; i < num_path; i++)
		{
			// delay_power[i] = num_ray * delay[i] * power[i];
			// power_ray_sum[i] = num_ray * power[i];

			delay_power[i] = delay_LOS[i] * power[i];
			power_ray_sum[i] = power[i];
		}

		Real avg_delay;
		Real delay_power_sum = 0;
		Real power_sum = 0;

		for (int i = 0; i < num_path; i++)
		{
			delay_power_sum += delay_power[i];
			power_sum += power_ray_sum[i];
		}

		avg_delay = delay_power_sum / power_sum;

		// RMS_delay_power = new Real[MAX_NUM_CLUSTERS];

		for (int i = 0; i < num_path; i++)
		{
			RMS_delay_power[i] = (delay_LOS[i] - avg_delay) * (delay_LOS[i] - avg_delay) * power[i];
		}

		// Real RMS_delay_spread = 0;
		Real RMS_delay_power_sum = 0;

		for (int i = 0; i < num_path; i++)
		{
			RMS_delay_power_sum += RMS_delay_power[i];
		}

		RMS_delay_spread = sqrt(RMS_delay_power_sum / power_sum);
	}
}

void CHANNEL::Set_circular_angle_spread()
{
	// Real mu_theta = 0;

	// 3GPP 25.996 Annex A
	for (int i = 0; i < num_path; i++)
	{
		// delay_power[i] = num_ray * delay[i] * power[i];
		// power_ray_sum[i] = num_ray * power[i];

		AOA_power[i] = AOA[i] * power[i];
		AOD_power[i] = AOD[i] * power[i];
		power_ray_sum_circular[i] = power[i];
	}

	Real mu_theta_AOA;
	Real mu_theta_AOD;
	Real AOA_power_sum = 0;
	Real AOD_power_sum = 0;
	Real power_sum = 0;

	// AOA AOD�� degree

	for (int i = 0; i < num_path; i++)
	{
		AOA_power_sum += AOA_power[i];
		AOD_power_sum += AOD_power[i];
		power_sum += power_ray_sum_circular[i];
	}

	mu_theta_AOA = AOA_power_sum / power_sum; // �̰͵� degree�ε�
	mu_theta_AOD = AOD_power_sum / power_sum;

	/*
	cout << "A = " << mu_theta_AOA << endl;
	cout << "D = " << mu_theta_AOD << endl;
	cout << pi << endl;
	getchar();
	*/

	// Calibration ������ �߰��� ����
	if (mu_theta_AOA < -180.) // pi��� �Ǿ��ִµ� degree�� ������� => -180 ~ 180���� �����ִ°���
	{
		mu_theta_AOA = 2. * 180. + mu_theta_AOA;
	}
	else if (mu_theta_AOA > 180.)
	{
		mu_theta_AOA = mu_theta_AOA - 2. * 180.;
	}
	else
	{
		mu_theta_AOA = mu_theta_AOA;
	}

	if (mu_theta_AOD < -180.)
	{
		mu_theta_AOD = 2. * 180. + mu_theta_AOD;
	}
	else if (mu_theta_AOD > 180.)
	{
		mu_theta_AOD = mu_theta_AOD - 2. * 180.;
	}
	else
	{
		mu_theta_AOD = mu_theta_AOD;
	}

	for (int i = 0; i < num_path; i++)
	{
		if ((AOA[i] - mu_theta_AOA) < -180.) // pi -> 180
		{
			theta_n_m_mu_AOA[i] = 2. * 180. + (AOA[i] - mu_theta_AOA);
		}
		else if (((AOA[i] - mu_theta_AOA) >= -180.) && ((AOA[i] - mu_theta_AOA) <= 180.))
		{
			theta_n_m_mu_AOA[i] = (AOA[i] - mu_theta_AOA);
		}
		else if ((AOA[i] - mu_theta_AOA) > 180.)
		{
			theta_n_m_mu_AOA[i] = 2. * 180. - (AOA[i] - mu_theta_AOA);
		}
	}

	for (int i = 0; i < num_path; i++)
	{
		if ((AOD[i] - mu_theta_AOD) < -180.)
		{
			theta_n_m_mu_AOD[i] = 2. * 180. + (AOD[i] - mu_theta_AOD);
		}
		else if (((AOD[i] - mu_theta_AOD) >= -180.) && ((AOD[i] - mu_theta_AOD) <= 180.))
		{
			theta_n_m_mu_AOD[i] = (AOD[i] - mu_theta_AOD);
		}
		else // if ((AOD[i] - mu_theta_AOD) > 180.)
		{
			theta_n_m_mu_AOD[i] = 2. * 180. - (AOD[i] - mu_theta_AOD);
		}
	}

	for (int i = 0; i < num_path; i++)
	{
		theta_power_AOA[i] = theta_n_m_mu_AOA[i] * theta_n_m_mu_AOA[i] * power[i];
		theta_power_AOD[i] = theta_n_m_mu_AOD[i] * theta_n_m_mu_AOD[i] * power[i];
	}

	Real theta_power_sum_AOA = 0;
	Real theta_power_sum_AOD = 0;

	for (int i = 0; i < num_path; i++)
	{
		theta_power_sum_AOA += theta_power_AOA[i];
		theta_power_sum_AOD += theta_power_AOD[i];
	}

	circular_angle_spread_AOA = sqrt(theta_power_sum_AOA / power_sum);
	circular_angle_spread_AOD = sqrt(theta_power_sum_AOD / power_sum);
}

void CHANNEL::Precalculate(int _time_idx) // �̸� ����ϴ� CH ��Ʈ���� �κ�
{
	/*
	//// =============������¿�=================================================//
	//ofstream CHIR_Precalc("Precalc_CHIR.dat");
	//CHIR_Precalc << "time_idx" << "\t" << "LOS" << "\t" << "INDOOR" << "\t" << "V_ANGLE_MAPPING" << "\t" << "NUM_TX" << "\t" << "NUM_RX" << "\t" << "NUM_CLUSTER" << "\t" << "CHIR" << endl;
	//// =============������¿�=================================================//

	if (TYPE == 1 || TYPE == 2 || TYPE == 3 || TYPE == 4 || TYPE == 5)
	{
		//ms_idx = _ms_idx;

		Real ray_AOD = 0.; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
		Real ray_AOA = 0.;
		Real doppler_freq = 0.;


		//for (int time = 0; time < precalculate_time; time++)
		//{
			for (int los = 0; los < 2; los++)
			{
				for (int indoor = 0; indoor < 2; indoor++)
				{
					for (int v = 0; v < anglev_mapping_value; v++)
					{
						for (int u = 0; u < NUM_TX; u++)
						{
							for (int s = 0; s < NUM_RX; s++)
							{
								for (int i = 0; i < num_path; i++)
								{
									for (int j = 0; j < num_ray; j++)
									{
										if (los == 0)
										{
											ray_AOD = AOD_NLOS[i] + cluster_ASD * offset_angle[i][j]; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
											ray_AOA = AOA_NLOS[i] + cluster_ASA * offset_angle[i][j];
										}
										else if (los == 1)
										{
											ray_AOD = AOD_LOS[i] + cluster_ASD * offset_angle[i][j]; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
											ray_AOA = AOA_LOS[i] + cluster_ASA * offset_angle[i][j];
										}

										// indoor�� ���� user_speed �ٸ����� ��ġ��
										doppler_freq = (user_speed * 1000. / 3600.) * cos((ray_AOA - (2. * 180. * randnum.u() - 180.))*pi / 180) / Wavelength;



										Real xpr = pow(10, XPR / 10);


										// angle_v ��� �̾Ƽ� 10�� �� �ϳ� �����ϵ��� ��ġ��
										// ���� ��� angle_v = (10 20 30 40 ...) ���� �ϳ� ������ ���..


										//Real angle_v = atan2((_bs_height - ms_height), distance) * 180 / pi;

										//Real angle_v = atan2((_bs_height - ms_height), min_distance) * 180 / pi;

										//Real angle_v = 10. * (anglev_mapping_value + 1.);

										Real angle_v = ((MAX_v_angle / anglev_mapping_value) * (v + 1)) * randnum.u();

										if (angle_v < ((MAX_v_angle / anglev_mapping_value) * v))
										{
											angle_v = angle_v + ((MAX_v_angle / anglev_mapping_value) * v);
										}


										///////////F_tx_v
										F_tx_v = 12. * (angle_v - ANGLE_tilt)*(angle_v - ANGLE_tilt) / (15. * 15.); // 15 = pi_3dB // vertical(elevation) antenna pattern
										if (F_tx_v > 20.)
										{
											F_tx_v = 20.;
										}

										////////////F_tx_h
										F_tx_h = 12. * ray_AOD*ray_AOD / (70. * 70.); // 70 = theta_3dB // horizental antenna pattern
										if (F_tx_h > 20.)
										{
											F_tx_h = 20.;
										}

										//////////F_rx_v
										F_rx_v = 0.;

										//////////F_rx_h
										F_rx_h = 0.;

										F_vv = pow(10., -(F_rx_v + F_tx_v) / 10.);
										F_vh = pow(10., -(F_rx_v + F_tx_h) / 10.);
										F_hv = pow(10., -(F_rx_h + F_tx_v) / 10.);
										F_hh = pow(10., -(F_rx_h + F_tx_h) / 10.);



										if (TYPE == 1)
										{
											F_vv = 1;
											F_vh = 1;
											F_hv = 1;
											F_hh = 1;
										}



										ComplexReal jay(0, 1);
										ComplexReal  antenna_gain = F_vv*exp(jay*random_phase_vv[i][j]) + sqrt(1 / xpr)*F_hv*exp(jay*random_phase_hv[i][j]) + sqrt(1 / xpr)*F_vh*exp(jay*random_phase_vh[i][j]) + F_hh*exp(jay*random_phase_hh[i][j]);


										if (los == 1)  // LOS
										{
											CHIR_precalc[los][indoor][v][u][s][i] += sqrt(power_LOS[i] / num_ray) * antenna_gain * exp(jay * (d_tx * u) * 2. * pi * sin(ray_AOD*pi / 180.) / Wavelength) * exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA*pi / 180.) / Wavelength) * exp(jay * 2. * pi * doppler_freq * (Real)(_time_idx / 1000.));
											// d_tx �� uniform distance between tx element (m.2135-1)
										}
										else           // NLOS
										{
											CHIR_precalc[los][indoor][v][u][s][i] += sqrt(power_NLOS[i] / num_ray) * antenna_gain * exp(jay * (d_tx * u) * 2. * pi * sin(ray_AOD*pi / 180.) / Wavelength) * exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA*pi / 180.) / Wavelength) * exp(jay * 2. * pi * doppler_freq * (Real)(_time_idx / 1000.));
											// d_tx �� uniform distance between tx element (m.2135-1)
										}


										if (los == 1)
										{
											Real k_factor = pow(10, K_factor / 10);

											//CHIR_precalculate[time][los][indoor][v][u][s][i] *= sqrt(1 / (k_factor + 1));
											CHIR_precalc[los][indoor][v][u][s][i] *= sqrt(1 / (k_factor + 1));

											Real ray_delay = delay[0];
											Real ray_AOD_los = AOD_LOS[0];
											Real ray_AOA_los = AOA_LOS[0];

											// �����ʿ�
											doppler_freq = (user_speed * 1000. / 3600.) * cos(ray_AOA_los - (2. * 180. * randnum.u() - 180.)) / Wavelength;

											///////////F_tx_v
											F_tx_v = 12. * (angle_v - ANGLE_tilt)*(angle_v - ANGLE_tilt) / (15. * 15.);
											if (F_tx_v > 20.)
											{
												F_tx_v = 20.;
											}

											////////////F_tx_h
											F_tx_h = 12. * ray_AOD_los * ray_AOD_los / (70. * 70.);
											if (F_tx_h > 20.)
											{
												F_tx_h = 20.;
											}


											//////////F_rx_v
											F_rx_v = 0.;

											//////////F_rx_h
											F_rx_h = 0.;

											F_vv = pow(10., -(F_rx_v + F_tx_v) / 10.);
											F_vh = pow(10., -(F_rx_v + F_tx_h) / 10.);
											F_hv = pow(10., -(F_rx_h + F_tx_v) / 10.);
											F_hh = pow(10., -(F_rx_h + F_tx_h) / 10.);



											if (TYPE == 1)
											{
												F_vv = 1;
												F_vh = 1;
												F_hv = 1;
												F_hh = 1;
											}


											// LOS�϶� delta(n-1) * .... �����ִ°� �ִ�. i = 0 �϶��� �����ִ°�.. i �� cluster index
											if (i == 0)
											{
												//CHIR_precalculate[time][los][indoor][v][u][s][0] += sqrt(k_factor / (k_factor + 1)) * antenna_gain *exp(jay * (d_tx * u) * 2. * pi*sin(ray_AOD_los*pi / 180) / Wavelength)*exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA_los*pi / 180) / Wavelength)* exp(jay*2.*pi*doppler_freq * (Real)(time/1000.));
												CHIR_precalc[los][indoor][v][u][s][0] += sqrt(k_factor / (k_factor + 1)) * antenna_gain *exp(jay * (d_tx * u) * 2. * pi*sin(ray_AOD_los*pi / 180.) / Wavelength)*exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA_los*pi / 180.) / Wavelength)* exp(jay*2.*pi*doppler_freq * (Real)(_time_idx / 1000.));
											}
										}
									}
									//// =============������¿�=================================================//
									//CHIR_Precalc << time << "\t" << los << "\t" << indoor << "\t" << v << "\t" << u << "\t" << s << "\t" << i << "\t" << CHIR_precalc[time][los][indoor][v][u][s][i] << endl;
									//// =============������¿�=================================================//
								}
							}
						}
					}
				}
			}

			//cout << "CHIR Time = " << time + 1 << "ms" << endl;

		//}
	}



	///////////////////////////
	///
	/// 3D channel model
	///
	///////////////////////////
	else
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//ms_idx = _ms_idx;

		Real ray_AOD = 0.; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
		Real ray_AOA = 0.;
		Real doppler_freq = 0.;
		Real ray_ZOA = 0.;
		Real ray_ZOD = 0.;


		//for (int time = 0; time < precalculate_time; time++)
		//{
			for (int los = 0; los < 2; los++)
			{
				for (int indoor = 0; indoor < 2; indoor++)
				{
					for (int v = 0; v < anglev_mapping_value; v++)
					{
						for (int u = 0; u < NUM_TX; u++)
						{
							for (int s = 0; s < NUM_RX; s++)
							{
								for (int i = 0; i < num_path; i++)
								{
									for (int j = 0; j < num_ray; j++)
									{
										if (los == 0)
										{
											ray_AOD = AOD_NLOS[i] + cluster_ASD * offset_angle[i][j]; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
											ray_AOA = AOA_NLOS[i] + cluster_ASA * offset_angle[i][j];

											ray_ZOA = ZOA_NLOS[i] + cluster_ZSA * offset_angle[i][j];
											ray_ZOD = ZOD_NLOS[i] + (3 / 8) * pow(10, mu_ZSD) * offset_angle[i][j];
										}
										else if (los == 1)
										{
											ray_AOD = AOD_LOS[i] + cluster_ASD * offset_angle[i][j]; //ZOD�� ���� ZOD[i] + (3/8) * pow(10,mu_ZSD) * offset_angle[i][j]
											ray_AOA = AOA_LOS[i] + cluster_ASA * offset_angle[i][j];

											ray_ZOA = ZOA_LOS[i] + cluster_ZSA * offset_angle[i][j];
											ray_ZOD = ZOD_LOS[i] + (3 / 8) * pow(10, mu_ZSD) * offset_angle[i][j];
										}

										// indoor�� ���� user_speed �ٸ����� ��ġ��
										doppler_freq = (user_speed * 1000. / 3600.) * (sin(ray_ZOA * pi / 180) * cos(ray_AOA*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180) * cos((2. * 180. * randnum.u() - 180.)*pi / 180))
											+ (user_speed * 1000. / 3600.) * (sin(ray_ZOA*pi / 180) * sin(ray_AOA*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180))
											+ (user_speed * 1000. / 3600.) * (cos(ray_ZOA*pi / 180) * cos((2. * 180. * randnum.u() - 180.)*pi / 180));

										doppler_freq = doppler_freq / Wavelength;

										Real xpr = pow(10, XPR / 10);

										///////////F_tx_v

										Real ant_ZOD_v = 12. * (ray_ZOD - 90)*(ray_ZOD - 90) / (65. * 65.);
										Real ant_AOD_v = 12. * (ray_AOD - 90)*(ray_AOD - 90) / (65. * 65.);

										F_tx_v = -MIN(-(ant_ZOD_v + ant_AOD_v), 30.);

										////////////F_tx_h

										Real ant_ZOD_h = 12. * (ray_ZOD)*(ray_ZOD) / (65. * 65.);
										Real ant_AOD_h = 12. * (ray_AOD)*(ray_AOD) / (65. * 65.);

										F_tx_h = -MIN(-(ant_ZOD_h + ant_AOD_h), 30.);

										//////////F_rx_v
										F_rx_v = 0.;

										//////////F_rx_h
										F_rx_h = 0.;

										F_vv = pow(10., -(F_rx_v + F_tx_v) / 10.);
										F_vh = pow(10., -(F_rx_v + F_tx_h) / 10.);
										F_hv = pow(10., -(F_rx_h + F_tx_v) / 10.);
										F_hh = pow(10., -(F_rx_h + F_tx_h) / 10.);


										ComplexReal jay(0, 1);
										ComplexReal  antenna_gain = F_vv*exp(jay*random_phase_vv[i][j]) + sqrt(1 / xpr)*F_hv*exp(jay*random_phase_hv[i][j]) + sqrt(1 / xpr)*F_vh*exp(jay*random_phase_vh[i][j]) + F_hh*exp(jay*random_phase_hh[i][j]);

										if (los == 1)  // LOS
										{
											CHIR_precalc[los][indoor][v][u][s][i] += sqrt(power_LOS[i] / num_ray) * antenna_gain * exp(jay * (d_tx * u) * 2. * pi * sin(ray_AOD*pi / 180.) / Wavelength) * exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA*pi / 180.) / Wavelength) * exp(jay * 2. * pi * doppler_freq * (Real)(_time_idx / 1000.));

											// d_tx �� uniform distance between tx element (m.2135-1)
											//d_tx, d_rx �� �̻��ϴ�. �����ڵ带 ���� d_tx = 10 * wavelength * num_tx_antenna�ε� // num_tx_antenna�� ���������� �ʴ�(�����ڵ尡), �����ݿ� �ʿ��ϴ� // �����ߴ�

											//���⼭ d_tx�� location vector of tx antenna element (36.873) // m.2135-1�� �Ȱ��� ����ϵ��� �س��µ� �ٸ����� �ִ�(Ȯ��ġ ���� �κ�..)
											//  exp �����ִ� �κ� ���� �κ� ray_AOD ��� Ȯ��
										}
										else           // NLOS
										{
											CHIR_precalc[los][indoor][v][u][s][i] += sqrt(power_NLOS[i] / num_ray) * antenna_gain * exp(jay * (d_tx * u) * 2. * pi * sin(ray_AOD*pi / 180.) / Wavelength) * exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA*pi / 180.) / Wavelength) * exp(jay * 2. * pi * doppler_freq * (Real)(_time_idx / 1000.));

											// d_tx �� uniform distance between tx element (m.2135-1)
											//d_tx, d_rx �� �̻��ϴ�. �����ڵ带 ���� d_tx = 10 * wavelength * num_tx_antenna�ε� // num_tx_antenna�� ���������� �ʴ�(�����ڵ尡), �����ݿ� �ʿ��ϴ� // �����ߴ�

											//���⼭ d_tx�� location vector of tx antenna element (36.873) // m.2135-1�� �Ȱ��� ����ϵ��� �س��µ� �ٸ����� �ִ�(Ȯ��ġ ���� �κ�..)
											//  exp �����ִ� �κ� ���� �κ� ray_AOD ��� Ȯ��
										}



										if (los == 1)
										{
											Real k_factor = pow(10, K_factor / 10);

											//CHIR_precalculate[time][los][indoor][v][u][s][i] *= sqrt(1 / (k_factor + 1));

											CHIR_precalc[los][indoor][v][u][s][i] *= sqrt(1 / (k_factor + 1));

											Real ray_delay = delay[0];
											Real ray_AOD_los = AOD_LOS[0];
											Real ray_AOA_los = AOA_LOS[0];

											Real ray_ZOA_los = ZOA_LOS[0];
											Real ray_ZOD_los = ZOD_LOS[0];

											// �����ʿ�
											// indoor�� ���� user_speed �ٸ����� ��ġ��
											doppler_freq = user_speed * (sin(ray_ZOA * pi / 180) * cos(ray_AOA*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180) * cos((2. * 180. * randnum.u() - 180.)*pi / 180))
												+ user_speed * (sin(ray_ZOA*pi / 180) * sin(ray_AOA*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180) * sin((2. * 180. * randnum.u() - 180.)*pi / 180))
												+ user_speed * (cos(ray_ZOA*pi / 180) * cos((2. * 180. * randnum.u() - 180.)*pi / 180));

											doppler_freq = doppler_freq / Wavelength;

											///////////F_tx_v

											Real ant_ZOD_v = 12. * (ray_ZOD - 90)*(ray_ZOD - 90) / (65. * 65.);
											Real ant_AOD_v = 12. * (ray_AOD - 90)*(ray_AOD - 90) / (65. * 65.);

											F_tx_v = -MIN(-(ant_ZOD_v + ant_AOD_v), 30.);

											////////////F_tx_h

											Real ant_ZOD_h = 12. * (ray_ZOD)*(ray_ZOD) / (65. * 65.);
											Real ant_AOD_h = 12. * (ray_AOD)*(ray_AOD) / (65. * 65.);

											F_tx_h = -MIN(-(ant_ZOD_h + ant_AOD_h), 30.);

											//////////F_rx_v
											F_rx_v = 0.;

											//////////F_rx_h
											F_rx_h = 0.;

											F_vv = pow(10., -(F_rx_v + F_tx_v) / 10.);
											F_vh = pow(10., -(F_rx_v + F_tx_h) / 10.);
											F_hv = pow(10., -(F_rx_h + F_tx_v) / 10.);
											F_hh = pow(10., -(F_rx_h + F_tx_h) / 10.);



											// LOS�϶� delta(n-1) * .... �����ִ°� �ִ�. i =0 �϶��� �����ִ°�.. i �� cluster index
											if (i == 0)
											{
												//CHIR_precalculate[time][los][indoor][v][u][s][0] += sqrt(k_factor / (k_factor + 1)) * antenna_gain *exp(jay * (d_tx * u) * 2. * pi*sin(ray_AOD_los*pi / 180) / Wavelength)*exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA_los*pi / 180) / Wavelength)* exp(jay*2.*pi*doppler_freq * (Real)(time/1000.));
												CHIR_precalc[los][indoor][v][u][s][0] += sqrt(k_factor / (k_factor + 1)) * antenna_gain *exp(jay * (d_tx * u) * 2. * pi*sin(ray_AOD_los*pi / 180.) / Wavelength)*exp(jay * (d_rx * s) * 2. * pi*sin(ray_AOA_los*pi / 180.) / Wavelength)* exp(jay*2.*pi*doppler_freq * (Real)(_time_idx / 1000.));
											}
										}
									}

									//// =============������¿�=================================================//
									//CHIR_Precalc << time << "\t" << los << "\t" << indoor << "\t" << v << "\t" << u << "\t" << s << "\t" << i << "\t" << CHIR_precalc[time][los][indoor][v][u][s][i] << endl;
									//// =============������¿�=================================================//
								}
							}
						}
					}
				}
			}

			//cout << "CHIR Time = " << time + 1 << "ms" << endl;
		//}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	//// =============������¿�=================================================//
	//CHIR_Precalc.close();
	//// =============������¿�=================================================//
	*/
}

void CHANNEL::Load_precalculate(Real t, int _ms_idx, int _initial_time_offset)
{
	ms_idx = _ms_idx;

	int los_idx;
	int indoor_idx;
	int v_idx;

	if (LOS == 1)
	{
		los_idx = 1;
	}
	else
	{
		los_idx = 0;
	}

	if (Indoor == true) ////////
	{
		indoor_idx = 1;
	}
	else
	{
		indoor_idx = 0;
	}

	if (TYPE == 1 || TYPE == 2 || TYPE == 3 || TYPE == 4 || TYPE == 5)
	{
		Real angle_v = atan2((_bs_height - ms_height), distance) * 180. / pi;

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
		// v_idx = (int)(anglev_mapping_value * randnum.u());
		v_idx = 0; //
	}
}

void CHANNEL::Update(Real t, int _ms_idx, int adj_sector)
{
	/*//////////////////////////////////////////////////////////
	// INPUT:
	// OUTPUT: CHIR
	// Notes: Channel coefficient generation - M.2412 page.59
	*/
	//////////////////////////////////////////////////////////

	// 220809 jhnoh
	int _bs_idx;
	if (TYPE == 11 && num_Indoor_TRxP == 1) //// Get all BS-MS channel coefficient
		_bs_idx = adj_sector;
	else
		_bs_idx = (int)(adj_sector / 3);

	// CHIR 메모리 지연 할당 (필요한 섹터만)
	int sector_num_idx = adj_sector % 3;
	Allocate_CHIR_memory(sector_num_idx);

	Real F_theta_GCS[2]; // 2 polar
	Real F_pi_GCS[2];
	Real NLOS_F_theta_GCS[2]; // 2 polar
	Real NLOS_F_pi_GCS[2];

	Real UE_F_theta_GCS[2]; // 2
	Real UE_F_pi_GCS[2];
	Real UE_NLOS_F_theta_GCS[2]; // 2
	Real UE_NLOS_F_pi_GCS[2];

	complex<Real> F_tx_theta = {0, 0};
	complex<Real> F_tx_pi = {0, 0};

	complex<Real> NLOS_F_tx_theta = {0, 0};
	complex<Real> NLOS_F_tx_pi = {0, 0};

	LOCATION3D r_tx;
	LOCATION3D d_tx;
	LOCATION3D r_rx;
	LOCATION3D d_rx;
	complex<Real> w;
	complex<Real> weight = 0;
	complex<Real> jay(0.0, 1.0);
	complex<Real> F_theta_temp(0, 0);
	complex<Real> F_pi_temp(0, 0);

	int K = BS_M / BS_Mp;	 /////   4 / 4 = 1
	int L = BS_N / BS_Np;	 //////   4/4 = 1
	int ue_K = MS_M / MS_Mp; /////   4 / 4 = 1
	int ue_L = MS_N / MS_Np; //////   4/4 = 1

	//////////////////////////////////// RX LOS
	for (int m_idx = 0; m_idx < MS_Mp; m_idx++)
	{
		for (int n_idx = 0; n_idx < MS_Np; n_idx++)
		{
			for (int p = 0; p < 2; p++)
			{
				Real UE_LOS_h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOA_GCS * (pi / 180.);	   //// pi in GCS, rad
				Real UE_LOS_v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOA_GCS * (pi / 180.); ///// theta in GCS, rad

				Real UE_LOS_combined_antenna_gain = Get_UE_antenna_pattern(p, UE_LOS_v_angle_theta, UE_LOS_h_angle_pi, _ms_idx, 0, UE_F_theta_GCS[0], UE_F_pi_GCS[0], UE_F_theta_GCS[1], UE_F_pi_GCS[1]);

				r_rx.x = sin(UE_LOS_v_angle_theta) * cos(UE_LOS_h_angle_pi);
				r_rx.y = sin(UE_LOS_v_angle_theta) * sin(UE_LOS_h_angle_pi);
				r_rx.z = cos(UE_LOS_v_angle_theta);

				weight = {0, 0};

				for (int k = 0; k < ue_K; k++) // vertical element
				{
					for (int l = 0; l < ue_L; l++) // horizental element
					{

						// w = ue_virtualization_weight_wv[links[_ms_idx].zenith_angle_idx_selected_for_interference][links[_ms_idx].azimuth_angle_idx_selected_for_interference][k][l];

						// Generate channel for comp signal, not
						w = ue_virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].z]
													   [links[_ms_idx].analog_beam_selection[adj_sector].a]
													   [k][l];

						// d_rx = ms[_ms_idx].d_rx[m_idx*ue_K + k][n_idx*ue_L + l][p][0][links[_ms_idx].panel_idx_selected_for_interference];  //// sector(rx) [M][N][P][Mg][Ng]
						d_rx = ms[_ms_idx].d_rx[m_idx * ue_K + k][n_idx * ue_L + l][p][0][links[_ms_idx].analog_beam_selection[adj_sector].p];
						weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_rx, d_rx));
					}
				}
				F_theta_temp = (UE_F_theta_GCS[p] * weight);
				F_pi_temp = (UE_F_pi_GCS[p] * weight);

				if (ue_antenna_element_gain == 0)
				{
					F_theta_temp = (UE_F_theta_GCS[p]);
					F_pi_temp = (UE_F_pi_GCS[p]);
				}

				RX_LOS_gain_theta[m_idx * MS_Np * 2 + n_idx * 2 + p] = F_theta_temp;
				RX_LOS_gain_pi[m_idx * MS_Np * 2 + n_idx * 2 + p] = F_pi_temp;
			}
		}
	}

	//////////////////////////////////// RX NLOS
	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
		{
			for (int m_idx = 0; m_idx < MS_Mp; m_idx++)
			{
				for (int n_idx = 0; n_idx < MS_Np; n_idx++)
				{
					for (int p = 0; p < 2; p++)
					{
						Real UE_NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].AOA[i] + channel[_bs_idx][_ms_idx].cluster_ASA * channel[_bs_idx][_ms_idx].offset_angle[i][j];
						Real UE_NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ZOA[i] + channel[_bs_idx][_ms_idx].cluster_ZSA * channel[_bs_idx][_ms_idx].offset_angle[i][j];

						Transform_angle_minus_180_to_plus_180(UE_NLOS_h_angle_pi);
						Transform_angle_0_to_plus_180(UE_NLOS_v_angle_theta);

						UE_NLOS_h_angle_pi = UE_NLOS_h_angle_pi * (pi / 180.);
						UE_NLOS_v_angle_theta = UE_NLOS_v_angle_theta * (pi / 180.);

						Real UE_NLOS_combined_antenna_gain = Get_UE_antenna_pattern(p, UE_NLOS_v_angle_theta, UE_NLOS_h_angle_pi, _ms_idx, 0, UE_NLOS_F_theta_GCS[0], UE_NLOS_F_pi_GCS[0], UE_NLOS_F_theta_GCS[1], UE_NLOS_F_pi_GCS[1]);

						r_rx.x = sin(UE_NLOS_v_angle_theta) * cos(UE_NLOS_h_angle_pi);
						r_rx.y = sin(UE_NLOS_v_angle_theta) * sin(UE_NLOS_h_angle_pi);
						r_rx.z = cos(UE_NLOS_v_angle_theta);

						weight = {0, 0};

						for (int k = 0; k < ue_K; k++) // vertical element
						{
							for (int l = 0; l < ue_L; l++) // horizental element
							{
								// w = ue_virtualization_weight_wv[links[_ms_idx].zenith_angle_idx_selected_for_interference][links[_ms_idx].azimuth_angle_idx_selected_for_interference][k][l];
								w = ue_virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].z]
															   [links[_ms_idx].analog_beam_selection[adj_sector].a][k][l];

								// d_rx = ms[_ms_idx].d_rx[m_idx*ue_K + k][n_idx*ue_L + l][p][0][links[_ms_idx].panel_idx_selected_for_interference];  //// [M][N][P][Mg][Ng] <panel1>
								d_rx = ms[_ms_idx].d_rx[m_idx * ue_K + k][n_idx * ue_L + l][p][0][links[_ms_idx].analog_beam_selection[adj_sector].p];
								weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_rx, d_rx));
							}
						}
						F_theta_temp = (UE_NLOS_F_theta_GCS[p] * weight);
						F_pi_temp = (UE_NLOS_F_pi_GCS[p] * weight);

						if (ue_antenna_element_gain == 0)
						{
							F_theta_temp = (UE_NLOS_F_theta_GCS[p]);
							F_pi_temp = (UE_NLOS_F_pi_GCS[p]);
						}

						RX_NLOS_gain_theta[m_idx * MS_Np * 2 + n_idx * 2 + p][i][j] = F_theta_temp;
						RX_NLOS_gain_pi[m_idx * MS_Np * 2 + n_idx * 2 + p][i][j] = F_pi_temp;
					}
				}
			}
		}
	}

	// 220809 jhnoh
	// for (int sector_num_idx = 0; sector_num_idx < 3; sector_num_idx++)
	//{
	// sector_num_idx 이미 위에서 선언됨

	////////////////// TX LOS
	Real h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOD_GCS * (pi / 180.);	//// pi in GCS
	Real v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOD_GCS * (pi / 180.); ///// theta in GCS

	Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_num_idx, F_theta_GCS[0], F_pi_GCS[0], F_theta_GCS[1], F_pi_GCS[1]);

	// for (int tx_port_idx = 0; tx_port_idx < (NUM_TX_Port / 2.); tx_port_idx++)  // 1/2 => polar
	//{
	for (int p = 0; p < 2; p++)
	// for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
	{
		for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
		// for (int n_idx = 0; n_idx < BS_Np; n_idx++)
		{
			for (int n_idx = 0; n_idx < BS_Np; n_idx++)
			// for (int p = 0; p < 2; p++)
			{
				r_tx.x = sin(v_angle_theta) * cos(h_angle_pi);
				r_tx.y = sin(v_angle_theta) * sin(h_angle_pi);
				r_tx.z = cos(v_angle_theta);

				weight = {0, 0};

				for (int k = 0; k < K; k++) // vertical element
				{
					for (int l = 0; l < L; l++) // horizental element  links[_ms_idx].analog_beam_selection[adj_sector].a
					{
						// w = virtualization_weight_wv[links[_ms_idx].sector_zenith_angle_idx][links[_ms_idx].sector_azimuth_angle_idx][k][l];
						w = virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].sector_z]
													[links[_ms_idx].analog_beam_selection[adj_sector].sector_a][k][l];
						d_tx = bs[_bs_idx].d_tx[sector_num_idx][m_idx * K + k][n_idx * L + l][p][0][0]; //// sector(rx) [M][N][P][Mg][Ng]
						weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_tx, d_tx));
					}
				}
				F_theta_temp = (F_theta_GCS[p] * weight);
				F_pi_temp = (F_pi_GCS[p] * weight);

				// TX_LOS_gain_theta[m_idx*BS_Np*2 + n_idx*2 + p]  = F_theta_temp;
				TX_LOS_gain_theta[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_theta_temp;
				TX_LOS_gain_pi[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_pi_temp;
				// TX_LOS_gain_pi[m_idx*BS_Np * 2 + n_idx * 2 + p] = F_pi_temp;
			}
		}
		//}
	}

	//////////////////////////////////// TX NLOS

	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
		{
			Real NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].AOD[i] + channel[_bs_idx][_ms_idx].cluster_ASD * channel[_bs_idx][_ms_idx].offset_angle_rand_coupling[i][j];
			Real NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ZOD[i] + (3. / 8.) * pow(10, channel[_bs_idx][_ms_idx].mu_ZSD) * channel[_bs_idx][_ms_idx].offset_angle_rand_coupling[i][j];

			Transform_angle_minus_180_to_plus_180(NLOS_h_angle_pi);
			Transform_angle_0_to_plus_180(NLOS_v_angle_theta);

			NLOS_h_angle_pi = NLOS_h_angle_pi * (pi / 180.);
			NLOS_v_angle_theta = NLOS_v_angle_theta * (pi / 180.);

			Real NLOS_combined_antenna_gain = Get_BS_antenna_pattern(NLOS_v_angle_theta, NLOS_h_angle_pi, _bs_idx, sector_num_idx, NLOS_F_theta_GCS[0], NLOS_F_pi_GCS[0], NLOS_F_theta_GCS[1], NLOS_F_pi_GCS[1]);

			for (int p = 0; p < 2; p++)
			// for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
			{
				for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
				// for (int n_idx = 0; n_idx < BS_Np; n_idx++)
				{
					for (int n_idx = 0; n_idx < BS_Np; n_idx++)
					// for (int p = 0; p < 2; p++)
					{
						r_tx.x = sin(NLOS_v_angle_theta) * cos(NLOS_h_angle_pi);
						r_tx.y = sin(NLOS_v_angle_theta) * sin(NLOS_h_angle_pi);
						r_tx.z = cos(NLOS_v_angle_theta);

						// 220810 jhnoh
						weight = {0, 0};

						for (int k = 0; k < K; k++) // vertical element
						{
							for (int l = 0; l < L; l++) // horizental element
							{
								// w = virtualization_weight_wv[links[_ms_idx].sector_zenith_angle_idx][links[_ms_idx].sector_azimuth_angle_idx][k][l];
								w = virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].sector_z]
															[links[_ms_idx].analog_beam_selection[adj_sector].sector_a][k][l];
								d_tx = bs[_bs_idx].d_tx[sector_num_idx][m_idx * K + k][n_idx * L + l][p][0][0]; //// sector(rx) [M][N][P][Mg][Ng]
								weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_tx, d_tx));
							}
						}
						F_theta_temp = (NLOS_F_theta_GCS[p] * weight);
						F_pi_temp = (NLOS_F_pi_GCS[p] * weight);

						TX_NLOS_gain_theta[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx][i][j] = F_theta_temp;
						TX_NLOS_gain_pi[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx][i][j] = F_pi_temp;
					}
				}
			}
		}
	}

	////////////// CHIR
	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			complex<Real> Big_PI_LOS(0., (360. * randnum.u() - 180.0) * (pi / 180.0));

			complex<Real> alpha_zero_temp = (RX_LOS_gain_theta[rp] * exp(Big_PI_LOS) * TX_LOS_gain_theta[tp] - RX_LOS_gain_pi[rp] * exp(Big_PI_LOS) * TX_LOS_gain_pi[tp]);

			CHIR_LOS[sector_num_idx][tp][rp] = alpha_zero_temp;

			for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
			{
				complex<Real> alpha_nmup(0, 0);
				for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
				{
					Real kappa = channel[_bs_idx][_ms_idx].kappa[i][j];
					Real _1_over_sqrt_K = 1.0 / sqrt(kappa);

					complex<Real> Big_pi_NLOS_thetatheta(0, (360. * randnum.u() - 180.0) * (pi / 180.0));
					complex<Real> Big_pi_NLOS_thetapi(0, (360. * randnum.u() - 180.0) * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pitheta(0, (360. * randnum.u() - 180.0) * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pipi(0, (360. * randnum.u() - 180.0) * (pi / 180.0));

					// panel 1
					complex<Real> alpha_nmup_temp =
					sqrt(channel[_bs_idx][_ms_idx].power[i] / (channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i])) *
					( ( RX_NLOS_gain_theta[rp][i][j] * exp(Big_pi_NLOS_thetatheta) +
					    RX_NLOS_gain_pi[rp][i][j] * _1_over_sqrt_K * exp(Big_pi_NLOS_pitheta) )
					  * TX_NLOS_gain_theta[tp][i][j] +
					  ( RX_NLOS_gain_theta[rp][i][j] * _1_over_sqrt_K * exp(Big_pi_NLOS_thetapi) +
					    RX_NLOS_gain_pi[rp][i][j] * exp(Big_pi_NLOS_pipi) )
					  * TX_NLOS_gain_pi[tp][i][j] );

					alpha_nmup = alpha_nmup + alpha_nmup_temp;

					// Store per-ray gain for ray-level precision mode
					if (USE_RAY_LEVEL_DOPPLER && ray_data_allocated)
					{
						ray_gain[sector_num_idx][tp][rp][i][j] = alpha_nmup_temp;

						// Store ray angles (AOA with offsets)
						Real ray_azimuth = channel[_bs_idx][_ms_idx].AOA[i] +
						                   channel[_bs_idx][_ms_idx].cluster_ASA *
						                   channel[_bs_idx][_ms_idx].offset_angle[i][j];
						Real ray_zenith = channel[_bs_idx][_ms_idx].ZOA[i] +
						                  channel[_bs_idx][_ms_idx].cluster_ZSA *
						                  channel[_bs_idx][_ms_idx].offset_angle[i][j];

						Transform_angle_minus_180_to_plus_180(ray_azimuth);
						Transform_angle_0_to_plus_180(ray_zenith);

						ray_AOA[i][j][0] = ray_azimuth;  // Azimuth
						ray_AOA[i][j][1] = ray_zenith;   // Zenith
					}
				}

				CHIR[sector_num_idx][tp][rp][i] = alpha_nmup;
				// COMMENTED OUT: CHIR_vec not allocated to save memory
				// if (sampled_delay[i] < 4096)
				// 	CHIR_vec[sector_num_idx][tp][rp][sampled_delay[i]] = alpha_nmup;
			}
		}
	}
	//}
}

void CHANNEL::Update_per_time(Real t, int adj_sector, int _ms_idx)
{
	//////////////////////////////////// RX NLOS

	complex<Real> jay(0, 1);

	LOCATION3D r_rx;
	LOCATION3D v_rx;

	// 3GPP standard: slot duration = 1ms / 2^μ
	Real slot_duration = 1.0e-3 / pow(2.0, numerology);  // seconds
	Real time = t * slot_duration;

	// 220809 jhnoh
	int _bs_idx;
	if (TYPE == 11 && num_Indoor_TRxP == 1) //// Get all BS-MS channel coefficient
		_bs_idx = adj_sector;
	else
		_bs_idx = (int)(adj_sector / 3);

	// CHIR 메모리 지연 할당 (필요한 섹터만, 안전장치)
	int sector_num_idx = adj_sector % 3;
	if (!sector_allocated[sector_num_idx]) {
		Allocate_CHIR_memory(sector_num_idx);
	}

	// 220809 jhnoh
	// for (int sector_num_idx = 0; sector_num_idx < 3; sector_num_idx++)
	//{
	// sector_num_idx 이미 위에서 선언됨
	//cout << "Memory Allocated 0" << endl;
	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			// FIX: Use cluster center angles for Doppler calculation (3GPP TR 38.901)
			// Previously: inner ray loop modified CHIR[i] repeatedly, only last ray's Doppler remained
			// Now: Apply Doppler once per cluster using cluster center AOA/ZOA
			// Use MS moving direction for velocity vector (independent of LOS)
			Real moving_elevation = ms[_ms_idx].moving_direction * (pi / 180.);        // elevation angle (theta)
			Real moving_azimuth = ms[_ms_idx].moving_direction_azimuth * (pi / 180.);  // azimuth angle (phi)

			// Velocity vector in GCS using MS's actual moving direction
			v_rx.x = ms[_ms_idx].speed * sin(moving_elevation) * cos(moving_azimuth);
			v_rx.y = ms[_ms_idx].speed * sin(moving_elevation) * sin(moving_azimuth);
			v_rx.z = ms[_ms_idx].speed * cos(moving_elevation);

			for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
			{
				// Use cluster center angles (not individual ray angles with offsets)
				Real UE_NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].AOA[i] * (pi / 180.);
				Real UE_NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ZOA[i] * (pi / 180.);

				// Ray arrival direction vector
				r_rx.x = sin(UE_NLOS_v_angle_theta) * cos(UE_NLOS_h_angle_pi);
				r_rx.y = sin(UE_NLOS_v_angle_theta) * sin(UE_NLOS_h_angle_pi);
				r_rx.z = cos(UE_NLOS_v_angle_theta);

				Real r_v = (r_rx.x * v_rx.x) + (r_rx.y * v_rx.y) + (r_rx.z * v_rx.z);

				r_v = r_v / Wavelength;

				// Apply Doppler shift once per cluster (not per ray)
				CHIR[sector_num_idx][tp][rp][i] *= exp(jay * (Real)2. * pi * r_v * time);
				// COMMENTED OUT: CHIR_vec not allocated to save memory
				// if (sampled_delay[i] < 4096)
				// 	CHIR_vec[sector_num_idx][tp][rp][sampled_delay[i]] *= exp(jay * (Real)2. * pi * r_v * time);
			}
		}
	}
	//cout << "CHIR is calculated " << endl;

	//}

	if (Propagation == LOS_propagation)
	{
		// Use MS moving direction for velocity vector (same as NLOS)
		Real moving_elevation = ms[_ms_idx].moving_direction * (pi / 180.);
		Real moving_azimuth = ms[_ms_idx].moving_direction_azimuth * (pi / 180.);

		v_rx.x = ms[_ms_idx].speed * sin(moving_elevation) * cos(moving_azimuth);
		v_rx.y = ms[_ms_idx].speed * sin(moving_elevation) * sin(moving_azimuth);
		v_rx.z = ms[_ms_idx].speed * cos(moving_elevation);

		for (int tp = 0; tp < NUM_TX_Port; tp++)
		{
			for (int rp = 0; rp < NUM_RX_Port; rp++)
			{
				Real UE_LOS_h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOA_GCS * (pi / 180.);	   //// pi in GCS, rad
				Real UE_LOS_v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOA_GCS * (pi / 180.); ///// theta in GCS, rad

				// LOS ray arrival direction
				r_rx.x = sin(UE_LOS_v_angle_theta) * cos(UE_LOS_h_angle_pi);
				r_rx.y = sin(UE_LOS_v_angle_theta) * sin(UE_LOS_h_angle_pi);
				r_rx.z = cos(UE_LOS_v_angle_theta);

				Real r_v = (r_rx.x * v_rx.x) + (r_rx.y * v_rx.y) + (r_rx.z * v_rx.z);

				r_v = r_v / Wavelength;

				CHIR_LOS[sector_num_idx][tp][rp] *= exp(jay * (Real)2. * pi * r_v * time);
				CHIR    [sector_num_idx][tp][rp][0] = sqrt(1 / (K_linear + 1)) * CHIR[sector_num_idx][tp][rp][0] + sqrt(K_linear / (K_linear + 1)) * CHIR_LOS[sector_num_idx][tp][rp];

				// COMMENTED OUT: CHIR_vec not allocated to save memory
				// CHIR_vec[sector_num_idx][tp][rp][0] = CHIR[sector_num_idx][tp][rp][0];
			}
		}
		//}
	}
	//cout << "LOS CHIR is calculated " << endl;
}

// Ray-level precision mode: Apply independent Doppler shift to each ray
void CHANNEL::Update_per_time_precise(Real t, int adj_sector, int _ms_idx)
{
	if (!USE_RAY_LEVEL_DOPPLER || !ray_data_allocated)
	{
		// Fall back to cluster-average mode
		Update_per_time(t, adj_sector, _ms_idx);
		return;
	}

	complex<Real> jay(0, 1);
	LOCATION3D r_rx;
	LOCATION3D v_rx;

	// 3GPP standard: slot duration = 1ms / 2^μ
	Real slot_duration = 1.0e-3 / pow(2.0, numerology);  // seconds
	Real time = t * slot_duration;

	// Determine BS index
	int _bs_idx;
	if (TYPE == 11 && num_Indoor_TRxP == 1)
		_bs_idx = adj_sector;
	else
		_bs_idx = (int)(adj_sector / 3);

	int sector_num_idx = adj_sector % 3;
	if (!sector_allocated[sector_num_idx]) {
		Allocate_CHIR_memory(sector_num_idx);
	}

	// Use MS moving direction for velocity vector (independent of LOS)
	Real moving_elevation = ms[_ms_idx].moving_direction * (pi / 180.);        // elevation angle (theta)
	Real moving_azimuth = ms[_ms_idx].moving_direction_azimuth * (pi / 180.);  // azimuth angle (phi)

	// Velocity vector in GCS using MS's actual moving direction
	v_rx.x = ms[_ms_idx].speed * sin(moving_elevation) * cos(moving_azimuth);
	v_rx.y = ms[_ms_idx].speed * sin(moving_elevation) * sin(moving_azimuth);
	v_rx.z = ms[_ms_idx].speed * cos(moving_elevation);

	//////////////////////////////////// NLOS: Ray-by-ray Doppler update
	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			// Process each cluster
			for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
			{
				// Sum over rays with independent Doppler shifts
				complex<Real> cluster_sum(0, 0);

				for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
				{
					// Get stored ray angles
					Real ray_azimuth = ray_AOA[i][j][0] * (pi / 180.);
					Real ray_zenith = ray_AOA[i][j][1] * (pi / 180.);

					// Compute ray direction vector
					r_rx.x = sin(ray_zenith) * cos(ray_azimuth);
					r_rx.y = sin(ray_zenith) * sin(ray_azimuth);
					r_rx.z = cos(ray_zenith);

					// Compute Doppler frequency: r̂·v̄ / λ
					Real r_v = (r_rx.x * v_rx.x + r_rx.y * v_rx.y + r_rx.z * v_rx.z) / Wavelength;

					// Apply Doppler phase shift to this ray's gain
					complex<Real> ray_contribution = ray_gain[sector_num_idx][tp][rp][i][j] *
					                                  exp(jay * (Real)2. * pi * r_v * time);

					cluster_sum += ray_contribution;
				}

				// Update cluster CHIR with sum of all ray contributions
				CHIR[sector_num_idx][tp][rp][i] = cluster_sum;
			}
		}
	}

	//////////////////////////////////// LOS component
	if (Propagation == LOS_propagation)
	{
		for (int tp = 0; tp < NUM_TX_Port; tp++)
		{
			for (int rp = 0; rp < NUM_RX_Port; rp++)
			{
				// LOS ray arrival direction
				Real UE_LOS_h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOA_GCS * (pi / 180.);
				Real UE_LOS_v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOA_GCS * (pi / 180.);

				r_rx.x = sin(UE_LOS_v_angle_theta) * cos(UE_LOS_h_angle_pi);
				r_rx.y = sin(UE_LOS_v_angle_theta) * sin(UE_LOS_h_angle_pi);
				r_rx.z = cos(UE_LOS_v_angle_theta);

				// LOS Doppler
				Real r_v = (r_rx.x * v_rx.x + r_rx.y * v_rx.y + r_rx.z * v_rx.z) / Wavelength;

				// Update LOS component
				CHIR_LOS[sector_num_idx][tp][rp] *= exp(jay * (Real)2. * pi * r_v * time);

				// Combine LOS and NLOS with K-factor (first cluster contains combined)
				CHIR[sector_num_idx][tp][rp][0] =
					sqrt(1 / (K_linear + 1)) * CHIR[sector_num_idx][tp][rp][0] +
					sqrt(K_linear / (K_linear + 1)) * CHIR_LOS[sector_num_idx][tp][rp];
			}
		}
	}
}

// Step 9 : Generate the cross polarization power ratios
void Generate_XPR(int ue_idx, int bs_idx, Real mu_XPR, Real sigma_XPR)
{
	for (int n = 0; n < MAX_NUM_CLUSTERS; n++)
	{
		for (int m = 0; m < MAX_NUM_RAYS; m++)
		{
			channel[bs_idx][ue_idx].kappa[n][m] = pow(10., (sigma_XPR * randnum.n() + mu_XPR) / 10.);
		}
	}
}

/*
Real Get_LCS_theta(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{
	return acos(cos(beta)*cos(gamma)*cos(GCS_theta) + (sin(beta)*cos(gamma)*cos(GCS_pi - alpha) - sin(gamma)*sin(GCS_pi - alpha))*sin(GCS_theta)) * 180. / pi;
}

Real Get_LCS_pi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{

	Real real;
	Real imag;

	real = (cos(beta)*sin(GCS_theta)*cos(GCS_pi - alpha) - sin(beta)*cos(GCS_theta));
	imag = (cos(beta)*sin(gamma)*cos(GCS_theta) + (sin(beta)*sin(gamma)*cos(GCS_pi - alpha) + cos(gamma)*sin(GCS_pi - alpha))*sin(GCS_theta));

	ComplexReal result = (real, imag);

	return arg(result) * 180. / pi;

}
*/
Real Get_cos_psi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{
	Real A = cos(beta) * cos(gamma) * sin(GCS_theta) - (sin(beta) * cos(gamma) * cos(GCS_pi - alpha) - sin(gamma) * sin(GCS_pi - alpha)) * cos(GCS_theta);

	Real B = cos(beta) * cos(gamma) * cos(GCS_theta) - (sin(beta) * cos(gamma) * cos(GCS_pi - alpha) - sin(gamma) * sin(GCS_pi - alpha)) * sin(GCS_theta);

	Real cos_psi = A / sqrt(1 - B * B);

	return cos_psi;
}

Real Get_sin_psi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{
	Real A = sin(beta) * cos(gamma) * sin(GCS_pi - alpha) + sin(gamma) * cos(GCS_pi - alpha);

	Real B = cos(beta) * cos(gamma) * cos(GCS_theta) - (sin(beta) * cos(gamma) * cos(GCS_pi - alpha) - sin(gamma) * sin(GCS_pi - alpha)) * sin(GCS_theta);

	Real sin_psi = A / sqrt(1 - B * B);

	return sin_psi;
}