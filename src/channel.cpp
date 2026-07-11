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

Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2, int port_idx = -1);

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
	circular_angle_spread_ZOA = 0;
	circular_angle_spread_ZOD = 0;
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
			offset_angle_perm_zoa[cluster_idx][ray_idx] = 0;
			offset_angle_perm_zod[cluster_idx][ray_idx] = 0;
		}
	}

	// ray_AOA/ray_AOD 초기화
	if (ray_AOA != NULL)
	{
		for (int cluster_idx = 0; cluster_idx < MAX_NUM_CLUSTERS; cluster_idx++)
		{
			for (int ray_idx = 0; ray_idx < MAX_NUM_RAYS; ray_idx++)
			{
				ray_AOA[cluster_idx][ray_idx][0] = 0.0;
				ray_AOA[cluster_idx][ray_idx][1] = 0.0;
				ray_AOD[cluster_idx][ray_idx][0] = 0.0;
				ray_AOD[cluster_idx][ray_idx][1] = 0.0;
			}
		}
	}
	ray_angles_precomputed = false;

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
						CHIR_init[sec_idx][u][s][i] = complex<Real>(0, 0);
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

	// SNS reset
	sns_Pr_sns = 0;
	sns_any_limited = false;
	sns_vr_los.limited = false;
	sns_rsrp_power_atten_los = REAL(1.0);
	for (int i = 0; i < MAX_NUM_CLUSTERS; i++) {
		sns_vr[i].limited = false;
		sns_rsrp_power_atten[i] = REAL(1.0);
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
	Generate_VisibilityRegion();
	Compute_SNS_RSRP_Attenuation();
	Find_Strong2Clusters();
	Set_AOAAOD(_bs_idx, _ue_idx);
	Set_ZOAZOD(_bs_idx, _ue_idx);

	Generate_XPR(_ue_idx, _bs_idx, mu_XPR, sigma_XPR);

	Set_InitialPhase();
	Set_SUBCLUSTER();

	Precompute_ray_angles();          // per-ray angles (after subcluster split)
	Set_circular_angle_spread();      // 3GPP TR 25.996 Annex A (per-ray level)
	Set_RMS_delay_spread();           // power-weighted RMS delay spread (cluster + sub-cluster level)

	Sampling_DelaySpread();
}

void CHANNEL::Allocate_memory()
{
	// jhnoh 240109
	sampled_delay =  new int[fft_size];

	delay = new Real[MAX_NUM_CLUSTERS];
	delay_LOS = new Real[MAX_NUM_CLUSTERS];
	power = new Real[MAX_NUM_CLUSTERS];
	powerForAngles = new Real[MAX_NUM_CLUSTERS];  // K-factor applied power for angle generation only
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
	offset_angle_perm_zoa = new Real *[MAX_NUM_CLUSTERS];
	offset_angle_perm_zod = new Real *[MAX_NUM_CLUSTERS];

	for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
	{
		random_phase_vv[i] = new Real[MAX_NUM_RAYS];
		random_phase_vh[i] = new Real[MAX_NUM_RAYS];
		random_phase_hv[i] = new Real[MAX_NUM_RAYS];
		random_phase_hh[i] = new Real[MAX_NUM_RAYS];
		offset_angle[i] = new Real[MAX_NUM_RAYS];
		offset_angle_rand_coupling[i] = new Real[MAX_NUM_RAYS];
		offset_angle_perm_zoa[i] = new Real[MAX_NUM_RAYS];
		offset_angle_perm_zod[i] = new Real[MAX_NUM_RAYS];
	}

	// ray_AOA/ray_AOD: per-ray angles [cluster][ray][2]
	ray_AOA = new Real **[MAX_NUM_CLUSTERS];
	ray_AOD = new Real **[MAX_NUM_CLUSTERS];
	for (int cluster = 0; cluster < MAX_NUM_CLUSTERS; cluster++)
	{
		ray_AOA[cluster] = new Real *[MAX_NUM_RAYS];
		ray_AOD[cluster] = new Real *[MAX_NUM_RAYS];
		for (int ray = 0; ray < MAX_NUM_RAYS; ray++)
		{
			ray_AOA[cluster][ray] = new Real[2];
			ray_AOD[cluster][ray] = new Real[2];
			ray_AOA[cluster][ray][0] = 0.0;
			ray_AOA[cluster][ray][1] = 0.0;
			ray_AOD[cluster][ray][0] = 0.0;
			ray_AOD[cluster][ray][1] = 0.0;
		}
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
		CHIR_init = new complex<Real> ***[3];

		for (int i = 0; i < 3; i++) {
			CHIR[i] = NULL;
			// CHIR_vec[i] = NULL;  // COMMENTED OUT: Not used in optimized version
			CHIR_LOS[i] = NULL;
			CHIR_init[i] = NULL;
		}

		CHIR_allocated = true;
	}

	// 이미 이 섹터가 할당되어 있으면 스킵
	if (sector_allocated[sector_idx]) return;

	// 특정 섹터만 할당
	CHIR[sector_idx] = new complex<Real> **[NUM_TX_Port];
	// CHIR_vec[sector_idx] = new complex<Real> **[NUM_TX_Port];  // COMMENTED OUT: Not used in optimized version
	CHIR_LOS[sector_idx] = new complex<Real> *[NUM_TX_Port];
	CHIR_init[sector_idx] = new complex<Real> **[NUM_TX_Port];

	for (int u = 0; u < NUM_TX_Port; u++)
	{
		CHIR[sector_idx][u] = new complex<Real> *[NUM_RX_Port];
		// CHIR_vec[sector_idx][u] = new complex<Real> *[NUM_RX_Port];  // COMMENTED OUT: Not used in optimized version
		CHIR_LOS[sector_idx][u] = new complex<Real>[NUM_RX_Port];
		CHIR_init[sector_idx][u] = new complex<Real> *[NUM_RX_Port];

		for (int s = 0; s < NUM_RX_Port; s++)
		{
			CHIR[sector_idx][u][s] = new complex<Real>[MAX_NUM_CLUSTERS];
			CHIR_init[sector_idx][u][s] = new complex<Real>[MAX_NUM_CLUSTERS];
			// CHIR_vec[sector_idx][u][s] = new complex<Real>[fft_size];  // COMMENTED OUT: Not used in optimized version

			for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
			{
				CHIR[sector_idx][u][s][i] = complex<Real>(0, 0);
				CHIR_init[sector_idx][u][s][i] = complex<Real>(0, 0);
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

		ray_data_allocated = true;
	}
}

void CHANNEL::Delete_memory()
{
	delete[] delay;
	delete[] delay_LOS;
	delete[] power;
	delete[] powerForAngles;
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
		delete[] offset_angle_rand_coupling[i];   // was leaked before
		delete[] offset_angle_perm_zoa[i];
		delete[] offset_angle_perm_zod[i];
	}

	delete[] random_phase_vv;
	delete[] random_phase_vh;
	delete[] random_phase_hv;
	delete[] random_phase_hh;

	delete[] offset_angle;
	delete[] offset_angle_rand_coupling;
	delete[] offset_angle_perm_zoa;
	delete[] offset_angle_perm_zod;

	// ray_AOA/ray_AOD deallocation
	if (ray_AOA != NULL)
	{
		for (int cluster = 0; cluster < MAX_NUM_CLUSTERS; cluster++)
		{
			for (int ray = 0; ray < MAX_NUM_RAYS; ray++)
			{
				delete[] ray_AOA[cluster][ray];
				delete[] ray_AOD[cluster][ray];
			}
			delete[] ray_AOA[cluster];
			delete[] ray_AOD[cluster];
		}
		delete[] ray_AOA;
		delete[] ray_AOD;
		ray_AOA = NULL;
		ray_AOD = NULL;
	}
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
				delete[] CHIR_init[sec_idx][u][s];
				// delete[] CHIR_vec[sec_idx][u][s];  // COMMENTED OUT: CHIR_vec not allocated
			}
			delete[] CHIR[sec_idx][u];
			delete[] CHIR_init[sec_idx][u];
			// delete[] CHIR_vec[sec_idx][u];  // COMMENTED OUT: CHIR_vec not allocated
			delete[] CHIR_LOS[sec_idx][u];
		}
		delete[] CHIR[sec_idx];
		delete[] CHIR_init[sec_idx];
		// delete[] CHIR_vec[sec_idx];  // COMMENTED OUT: CHIR_vec not allocated
		delete[] CHIR_LOS[sec_idx];

		sector_allocated[sec_idx] = false;
	}

	delete[] CHIR;
	delete[] CHIR_init;
	// delete[] CHIR_vec;  // COMMENTED OUT: CHIR_vec not allocated
	delete[] CHIR_LOS;

	CHIR = NULL;
	CHIR_init = NULL;
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
					LOS_prob = ((18. / distance) + exp(-1 * distance / 63.) * (1 - (18. / distance))) * (1 + C * (5.0 / 4.0) * pow((distance / 100), 3) * exp(-1 * distance / 150.));
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
					LOS_prob = ((18. / d_out) + exp(-1 * d_out / 63.) * (1 - (18. / d_out))) * (1 + C * (5.0 / 4.0) * pow(d_out / 100, 3) * exp(-1 * d_out / 150.));
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
	
	// PathLoss
	//------------------------------------------------- IMT-2020 --------------------------------------------------//
	// IMT 2020 EVAL Page 39 ~

	// InH
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
	// Dense_Urban
	else if (TYPE == 12)
	{
		// Real n = 0;
		// Real N = floor((5. * randnum.u()) + 4.); // 4 ~ 8
		// n = floor(((N)* randnum.u()) + 1.); // 1 ~ N

		// ms_height_in = 3 * (n - 1) + 1.5;
		int lc_Urban_env_Config_Type = Configuration_Type;


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
			/*
			if (distance > 18.)
			{
				g = (5.0 / 4.0) * pow((distance / 100), 3) * exp(-1 * distance / 150);
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
			*/

			hE = 1;

			Real d_BP = 4 * (ms_height_in_channel - hE) * (_bs_height - hE) * carrier_freq / light_speed;
			Real distance_3d;
			distance_3d = sqrt(distance * distance + (_bs_height - ms_height_in_channel) * (_bs_height - ms_height_in_channel)); // 3D distance

			/////////////////////////////////////////////// Channel Model A
			if (Channel_Model_Type == 0)
			{
				//////////////////////////////// LOS
				if (LOS == 1)
				{
					if (carrier_freq <= 7 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
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
					if (carrier_freq <= 7 * pow(10, 9)) ///// 0.5GHz <= fc <= 6GHz
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
						pathloss = 32.4 + 40 * log10(distance_3d)  + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
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
						pathloss = 32.4 + 40 * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.) - 9.5 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
					pathloss_2 = 22.4 + 35.3 * log10(distance_3d) + 21.3 * log10(carrier_freq / 1000000000.) - 0.3 * (ms_height_in_channel - 1.5);

					pathloss = MAX(pathloss, pathloss_2);
					sigma_SF = 7.82;
				}

				if (Propagation == OUT2IN_propagation)
				{
					sigma_SF = 7;
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

				if (Propagation == OUT2IN_propagation)
				{
					sigma_SF = 7;
				}

			}

		}
		//----------------------------------- UMa -------------------------------------------------------//
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
				g = (5.0 / 4.0) * pow((distance / 100), 3) * exp(-1 * distance / 150);
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

			// Channel Model A
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
						pathloss = 28 + 40. * log10(distance_3d) 
						              + 20. * log10(carrier_freq / 1000000000.) 
									  -  9. * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
				}
				// NLOS — TR 38.901 V19.1.0 (replaces M.2135 sub-6 formula)
				else if (LOS == 0)
				{
					sigma_SF = 6.;
					if (distance <= d_BP)
					{
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
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
			// Channel Model B — TR 38.901 V19.1.0 (same as Model A)
			else if (Channel_Model_Type == 1)
			{
				// LOS
				if (LOS == 1)
				{
					sigma_SF = 4.;
					if (distance <= d_BP)
					{
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
						pathloss = 40. * log10(distance_3d) + 28 + 20. * log10(carrier_freq / 1000000000.) - 9 * log10(pow(d_BP, 2) + pow((_bs_height - ms_height_in_channel), 2));
					}
				}
				// NLOS
				else if (LOS == 0)
				{
					sigma_SF = 6.;
					if (distance <= d_BP)
					{
						pathloss = 28. + 22. * log10(distance_3d) + 20. * log10(carrier_freq / 1000000000.);
					}
					else if (d_BP < distance)
					{
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
				XPR_mean[LOS_propagation ] = 11;
				XPR_mean[NLOS_propagation] = 10;
				XPR_std [LOS_propagation ] =  4;
				XPR_std [NLOS_propagation] =  4;

				if (LOS == 1) // LOS
				{
					num_path = 15; // Number of clusters
					r_tau = 3.6;

					// K_factor = normal(7, 4);
					mu_K_factor = 7.;
					sigma_K_factor = 4.;
					K_factor  = pow(10, normal(mu_K_factor, sigma_K_factor) / 10); // linear

					mu_DS     = -7.70;
					sigma_DS  = 0.18;
					DS        = pow(10, normal(-7.70, 0.18)); // = sigma_tau in M.2135

					mu_ASD    = 1.6;
					sigma_ASD = 0.18;
					mu_ASA    = 1.62;
					sigma_ASA = 0.22;

					mu_ZSA    = 1.22;
					sigma_ZSA = 0.23;
					mu_ZSD    = 1.02;
					sigma_ZSD = 0.41;
					mu_offset_ZOD = 0;

					cluster_DS  = -1; /// N/A
					cluster_ASD =  5;
					cluster_ASA =  8;
					cluster_ZSA =  9;

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
				{
					// ===== TR 38.901 Table 7.5-6 Part 1: UMi-Street Canyon =====
					Real fc = (carrier_freq / 1000000000.);  // fc [GHz]

					if (Propagation == 1) // LOS
					{
						num_path = 12;
						r_tau = 3;

						mu_K_factor = 9.;
						sigma_K_factor = 5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10);

						if (channel_param_legacy) {
							// Legacy (e00) frequency-dependent formulas
							mu_DS = -7.14 - 0.24 * log10(1 + fc);
							sigma_DS = 0.38;
							mu_ASD = 1.21 - 0.05 * log10(1 + fc);
							sigma_ASD = 0.41;
							mu_ASA = 1.73 - 0.08 * log10(1 + fc);
							sigma_ASA = 0.28 + 0.014 * log10(1 + fc);
							mu_ZSA = 0.73 - 0.1 * log10(1 + fc);
							sigma_ZSA = 0.34 - 0.04 * log10(1 + fc);
						} else {
							// Current spec (V19/j10) frequency-dependent formulas
							mu_DS = -0.18 * log10(1 + fc) - 7.28;
							sigma_DS = 0.39;
							mu_ASD = -0.05 * log10(1 + fc) + 1.21;
							sigma_ASD = 0.08 * log10(1 + fc) + 0.29;
							mu_ASA = -0.07 * log10(1 + fc) + 1.66;
							sigma_ASA = 0.021 * log10(1 + fc) + 0.26;
							mu_ZSA = -0.11 * log10(1 + fc) + 0.81;
							sigma_ZSA = -0.03 * log10(1 + fc) + 0.29;
						}
						DS = pow(10, normal(mu_DS, sigma_DS));

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
						num_path = 19;
						r_tau = 2.1;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						if (channel_param_legacy) {
							// Legacy (e00) frequency-dependent formulas
							mu_DS = -6.83 - 0.24 * log10(1 + fc);
							sigma_DS = 0.28 + 0.16 * log10(1 + fc);
							mu_ASD = 1.53 - 0.23 * log10(1 + fc);
							sigma_ASD = 0.33 + 0.11 * log10(1 + fc);
							mu_ASA = 1.81 - 0.08 * log10(1 + fc);
							sigma_ASA = 0.3 + 0.05 * log10(1 + fc);
							mu_ZSA = 0.92 - 0.04 * log10(1 + fc);
							sigma_ZSA = 0.41 - 0.07 * log10(1 + fc);
						} else {
							// Current spec (V19/j10) frequency-dependent formulas
							mu_DS = -0.22 * log10(1 + fc) - 6.87;
							sigma_DS = 0.19 * log10(1 + fc) + 0.22;
							mu_ASD = -0.24 * log10(1 + fc) + 1.54;
							sigma_ASD = 0.10 * log10(1 + fc) + 0.33;
							mu_ASA = -0.07 * log10(1 + fc) + 1.76;
							sigma_ASA = 0.05 * log10(1 + fc) + 0.27;
							mu_ZSA = -0.03 * log10(1 + fc) + 0.92;
							sigma_ZSA = -0.05 * log10(1 + fc) + 0.35;
						}
						DS = pow(10, normal(mu_DS, sigma_DS));

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
						num_path = 12;
						r_tau = 2.2;

						K_factor = -1; // N/A
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;
						DS = pow(10, normal(mu_DS, sigma_DS));

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
					// ===== TR 38.901 V19.1.0 Table 7.5-6 Part 1: UMa (sub-6 GHz) =====
					// (Replaces M.2412 Channel Model A sub-6 parameters)
					Real fc = (carrier_freq / 1000000000.);  // fc [GHz]
					// Note: Table 7.5-7 ZSD/ZOD uses fc=6 for fc<6 GHz
					Real fc_zsd = MAX(6.0, fc);

					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12;
						r_tau = 2.5;

						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10);

						mu_DS = -7.067 - 0.0794 * log10(fc);       // M.2412: -7.03
						sigma_DS = 0.57 + 0.026 * log10(fc);       // M.2412: 0.66

						mu_ASD = 0.92;                              // M.2412: 1.15
						sigma_ASD = 0.31;                           // M.2412: 0.28
						mu_ASA = 1.76;                              // M.2412: 1.81
						sigma_ASA = 0.19;                           // M.2412: 0.2

						mu_ZSA = 0.96;
						sigma_ZSA = 0.15;
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.955 - 0.0963 * log10(fc);
							sigma_DS = 0.66;
							mu_ASD = 1.06 + 0.1114 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 1.81;
							sigma_ASA = 0.2;
							mu_ZSA = 0.95;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622); // M.2412: N/A
						cluster_ASD = channel_param_legacy ? 5.0 : 3.58;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20;
						r_tau = 2.3;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.47 - 0.134 * log10(fc);         // M.2412: -6.44
						sigma_DS = 0.39;

						mu_ASD = 1.09;                              // M.2412: 1.41
						sigma_ASD = 0.44;                           // M.2412: 0.28
						mu_ASA = 2.04 - 0.25 * log10(fc);          // M.2412: 1.87
						sigma_ASA = 0.17 - 0.03 * log10(fc);       // M.2412: 0.11

						mu_ZSA = -0.2856 * log10(fc) + 1.445;
						sigma_ZSA = 0.17;

						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = 7.66 * log10(fc_zsd) - 5.96 - pow(10, (0.208 * log10(fc_zsd) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc_zsd) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.28 - 0.204 * log10(fc);
							mu_ASD = 1.5 - 0.1144 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 2.08 - 0.27 * log10(fc);
							sigma_ASA = 0.11;
							mu_ZSA = -0.3236 * log10(fc) + 1.512;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = channel_param_legacy ? 2.0 : 1.8;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12;
						r_tau = 2.2;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;

						mu_ASD = channel_param_legacy ? 1.25 : 0.58;
						sigma_ASD = channel_param_legacy ? 0.42 : 0.7;
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
							mu_offset_ZOD = 7.66 * log10(fc_zsd) - 5.96 - pow(10, (0.208 * log10(fc_zsd) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc_zsd) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = 11;
						cluster_ASD = channel_param_legacy ? 5.0 : 1.8;
						cluster_ASA = 8;
						cluster_ZSA = 3;

						XPR = normal(9, 5);
						cluster_shadowing = 4;
					}
				}
				else if (carrier_freq > 6000000000) /// 6GHz < fc <= 100GHz
				{
					// ===== TR 38.901 V19.1.0 Table 7.5-6 Part 1: UMa (above-6 GHz) =====
					// (Replaces M.2412 Channel Model A above-6 parameters)
					Real fc = (carrier_freq / 1000000000.);  // fc [GHz]

					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12;
						r_tau = 2.5;

						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10);

						mu_DS = -7.067 - 0.0794 * log10(fc);       // M.2412: -6.955 - 0.0963*log10(fc)
						sigma_DS = 0.57 + 0.026 * log10(fc);       // M.2412: 0.66

						mu_ASD = 0.92;                              // M.2412: 1.06 + 0.1114*log10(fc)
						sigma_ASD = 0.31;                           // M.2412: 0.28
						mu_ASA = 1.76;                              // M.2412: 1.81
						sigma_ASA = 0.19;                           // M.2412: 0.2

						mu_ZSA = 0.96;                              // M.2412: 0.95
						sigma_ZSA = 0.15;                           // M.2412: 0.16
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.955 - 0.0963 * log10(fc);
							sigma_DS = 0.66;
							mu_ASD = 1.06 + 0.1114 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 1.81;
							sigma_ASA = 0.2;
							mu_ZSA = 0.95;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = channel_param_legacy ? 5.0 : 3.58;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20;
						r_tau = 2.3;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.47 - 0.134 * log10(fc);         // M.2412: -6.28 - 0.204*log10(fc)
						sigma_DS = 0.39;

						mu_ASD = 1.09;                              // M.2412: 1.5 - 0.1144*log10(fc)
						sigma_ASD = 0.44;                           // M.2412: 0.28
						mu_ASA = 2.04 - 0.25 * log10(fc);          // M.2412: 2.08 - 0.27*log10(fc)
						sigma_ASA = 0.17 - 0.03 * log10(fc);       // M.2412: 0.11

						mu_ZSA = -0.2856 * log10(fc) + 1.445;      // M.2412: -0.3236*log10(fc) + 1.512
						sigma_ZSA = 0.17;                           // M.2412: 0.16
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = 7.66 * log10(fc) - 5.96 - pow(10, (0.208 * log10(fc) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.28 - 0.204 * log10(fc);
							mu_ASD = 1.5 - 0.1144 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 2.08 - 0.27 * log10(fc);
							sigma_ASA = 0.11;
							mu_ZSA = -0.3236 * log10(fc) + 1.512;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = channel_param_legacy ? 2.0 : 1.8;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12;
						r_tau = 2.2;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;                              // M.2412: -6.63
						sigma_DS = 0.32;

						mu_ASD = channel_param_legacy ? 1.25 : 0.58;
						sigma_ASD = channel_param_legacy ? 0.42 : 0.7;
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

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = 11;
						cluster_ASD = channel_param_legacy ? 5.0 : 1.8;
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
				{
					// ===== TR 38.901 V19.1.0 Table 7.5-6 Part 1: UMa =====
					// (Replaces M.2412 Channel Model B parameters)
					Real fc = (carrier_freq / 1000000000);  // V19: use actual fc [GHz]
					// M.2412 used: fc = max(6, carrier_freq_GHz)

					XPR_mean[LOS_propagation] = 8;
					XPR_mean[NLOS_propagation] = 7;
					XPR_mean[OUT2IN_propagation] = 9;

					XPR_std[LOS_propagation] = 4;
					XPR_std[NLOS_propagation] = 3;
					XPR_std[OUT2IN_propagation] = 5;

					if (Propagation == 1) // LOS
					{
						num_path = 12;
						r_tau = 2.5;

						mu_K_factor = 9.;
						sigma_K_factor = 3.5;
						K_factor = pow(10, normal(mu_K_factor, sigma_K_factor) / 10);

						mu_DS = -7.067 - 0.0794 * log10(fc);       // M.2412: -6.955 - 0.0963*log10(fc)
						sigma_DS = 0.57 + 0.026 * log10(fc);       // M.2412: 0.66

						mu_ASD = 0.92;                              // M.2412: 1.06 + 0.1114*log10(fc)
						sigma_ASD = 0.31;                           // M.2412: 0.28
						mu_ASA = 1.76;                              // M.2412: 1.81
						sigma_ASA = 0.19;                           // M.2412: 0.2

						mu_ZSA = 0.96;                              // M.2412: 0.95
						sigma_ZSA = 0.15;                           // M.2412: 0.16
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.75);
						sigma_ZSD = 0.4;
						mu_offset_ZOD = 0;

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.955 - 0.0963 * log10(fc);
							sigma_DS = 0.66;
							mu_ASD = 1.06 + 0.1114 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 1.81;
							sigma_ASA = 0.2;
							mu_ZSA = 0.95;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = channel_param_legacy ? 5.0 : 3.58;
						cluster_ASA = 11;
						cluster_ZSA = 7;

						XPR = normal(8, 4);
						cluster_shadowing = 3;
					}
					else if (Propagation == 0) //// NLOS
					{
						num_path = 20;
						r_tau = 2.3;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.47 - 0.134 * log10(fc);         // M.2412: -6.28 - 0.204*log10(fc)
						sigma_DS = 0.39;

						mu_ASD = 1.09;                              // M.2412: 1.5 - 0.1144*log10(fc)
						sigma_ASD = 0.44;                           // M.2412: 0.28
						mu_ASA = 2.04 - 0.25 * log10(fc);          // M.2412: 2.08 - 0.27*log10(fc)
						sigma_ASA = 0.17 - 0.03 * log10(fc);       // M.2412: 0.11

						mu_ZSA = -0.2856 * log10(fc) + 1.445;      // M.2412: -0.3236*log10(fc) + 1.512
						sigma_ZSA = 0.17;                           // M.2412: 0.16
						mu_ZSD = MAX(-0.5, -2.1 * (distance / 1000) - 0.01 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5) + 0.9);
						sigma_ZSD = 0.49;
						mu_offset_ZOD = 7.66 * log10(fc) - 5.96 - pow(10, (0.208 * log10(fc) - 0.782) * log10(MAX(25, distance)) - 0.13 * log10(fc) + 2.03 - 0.07 * (ms[self_ms_idx].MS_HEIGHT_FINAL - 1.5));

						// Legacy: old TR 38.901 Table 7.5-6 (pre-V19)
						if (channel_param_legacy) {
							mu_DS = -6.28 - 0.204 * log10(fc);
							mu_ASD = 1.5 - 0.1144 * log10(fc);
							sigma_ASD = 0.28;
							mu_ASA = 2.08 - 0.27 * log10(fc);
							sigma_ASA = 0.11;
							mu_ZSA = -0.3236 * log10(fc) + 1.512;
							sigma_ZSA = 0.16;
						}

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = MAX(0.25, -3.4084 * log10(fc) + 6.5622);
						cluster_ASD = channel_param_legacy ? 2.0 : 1.8;
						cluster_ASA = 15;
						cluster_ZSA = 7;

						XPR = normal(7, 3);
						cluster_shadowing = 3;
					}
					else if (Propagation == 2) ///// OUT2IN
					{
						num_path = 12;
						r_tau = 2.2;

						K_factor = -1;
						mu_K_factor = -1;
						sigma_K_factor = -1;

						mu_DS = -6.62;
						sigma_DS = 0.32;

						mu_ASD = channel_param_legacy ? 1.25 : 0.58;
						sigma_ASD = channel_param_legacy ? 0.42 : 0.7;
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

						DS = pow(10, normal(mu_DS, sigma_DS));

						cluster_DS = 11;
						cluster_ASD = channel_param_legacy ? 5.0 : 1.8;
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

	} 


	for (int i = 0; i < num_path; i++)
	{
		delay[i] = -r_tau * _DS * log(randnum.u());
	}


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
		powerForAngles[i] = 0.;
		power_LOS[i] = 0.;
		power_NLOS[i] = 0.;
	} //// initialize

	// Step 6: Generate cluster powers using UNSCALED delays (ns-3: GenerateClusterPowers)
	for (int i = 0; i < num_path; i++)
	{
		power[i]        = exp(-delay[i] * (r_tau - 1.) / (r_tau * _DS)) * pow(10., -1 * normal(0, cluster_shadowing) / 10.);

		power_LOS[i]    = power[i];
		power_NLOS[i]   = power[i];

		sum_power      += power[i];
		sum_power_LOS  += power_LOS[i];
		sum_power_NLOS += power_NLOS[i];
	}

	// Numerical stability check: detect near-zero denominator
	const Real EPSILON = REAL(1e-12);

	if (LOS == 1)
	{
		// Normalize cluster power WITHOUT K-factor (ns-3: m_clusterPower)
		// K-factor is only applied to powerForAngles for angle generation (Eq. 7.5-8)
		if (sum_power < EPSILON)
		{
			cout << "WARNING: Near-zero sum_power detected in LOS power normalization!" << endl;
			cout << "  sum_power = " << sum_power << ", k_R = " << k_R << endl;
			cout << "  BS_idx = " << self_bs_idx << ", MS_idx = " << self_ms_idx << endl;
		}
		for (int i = 0; i < num_path; i++)
		{
			power[i] = power[i] / sum_power;  // Simple normalization, NO K-factor
		}

		// Create powerForAngles with K-factor (Eq. 7.5-8) — used ONLY for angle generation
		// ns-3: clusterPowerForAngles in RemoveWeakClusters()
		for (int i = 0; i < num_path; i++)
		{
			if (i == 0)
			{
				// LOS cluster: P'_1 = P_1/(K_R+1) + K_R/(K_R+1)
				powerForAngles[i] = power[i] / (k_R + 1) + k_R / (k_R + 1);
			}
			else
			{
				// NLOS clusters: P'_n = P_n/(K_R+1)
				powerForAngles[i] = power[i] / (k_R + 1);
			}
		}
	}
	else if (LOS == 0)
	{
		if (sum_power < EPSILON)
		{
			cout << "WARNING: Near-zero sum_power detected in NLOS power normalization!" << endl;
			cout << "  sum_power = " << sum_power << endl;
			cout << "  BS_idx = " << self_bs_idx << ", MS_idx = " << self_ms_idx << endl;
		}
		for (int i = 0; i < num_path; i++)
		{
			power[i] = power[i] / sum_power;
			powerForAngles[i] = power[i];  // NLOS: same as power
		}
	}

	// Remove clusters with power more than 25 dB below the maximum cluster power
	// Use powerForAngles for threshold (ns-3: RemoveWeakClusters uses clusterPowerForAngles)
	Real max_power = powerForAngles[0];
	for (int i = 1; i < num_path; i++)
	{
		if (powerForAngles[i] > max_power)
			max_power = powerForAngles[i];
	}

	Real threshold = max_power * pow(10.0, -25.0 / 10.0); // -25 dB relative to max

	int new_num_path = 0;
	for (int i = 0; i < num_path; i++)
	{
		if (powerForAngles[i] >= threshold)
		{
			if (new_num_path != i)
			{
				power[new_num_path]           = power[i];
				powerForAngles[new_num_path]  = powerForAngles[i];
				power_LOS[new_num_path]       = power_LOS[i];
				power_NLOS[new_num_path]      = power_NLOS[i];
				delay[new_num_path]           = delay[i];
				delay_LOS[new_num_path]       = delay_LOS[i];
			}
			new_num_path++;
		}
	}

	// Zero out removed entries
	for (int i = new_num_path; i < num_path; i++)
	{
		power[i]           = 0.0;
		powerForAngles[i]  = 0.0;
		power_LOS[i]       = 0.0;
		power_NLOS[i]      = 0.0;
		delay[i]           = 0.0;
		delay_LOS[i]       = 0.0;
	}

	num_path = new_num_path;

	// Step 5 (resumed): Replace delays with LOS-scaled delays for channel generation
	// ns-3: AdjustClusterDelaysForLosCondition — applied AFTER power generation
	// The scaled delays (delay/D) are used in Fourier Transform for frequency-domain conversion
	if (LOS == 1)
	{
		for (int i = 0; i < num_path; i++)
		{
			delay[i] = delay_LOS[i];  // Replace with scaled delay (= original / D)
		}
	}
}

// ====================================================================
// BS-side Spatial Non-Stationarity: Generate Visibility Regions
// j20 Section 7.6.14.1.3 (Stochastic VR Model)
// Called after Set_POWER() — uses power[] (linear) per cluster
// ====================================================================
static void sns_setup_corner(ClusterVR& vr, Real a, Real b, Real W, Real H)
{
	// 50% left/right, 50% top/bottom (independent)
	bool right = (randnum.u() < REAL(0.5));
	bool top   = (randnum.u() < REAL(0.5));

	vr.a = a;
	vr.b = b;

	if (!right && !top) {       // bottom-left
		vr.x0 = REAL(0.0);  vr.y0 = REAL(0.0);
		vr.xa = a;           vr.yb = b;
		vr.xA = W;           vr.yB = H;
	} else if (right && !top) { // bottom-right
		vr.x0 = W;           vr.y0 = REAL(0.0);
		vr.xa = W - a;       vr.yb = b;
		vr.xA = REAL(0.0);   vr.yB = H;
	} else if (!right && top) { // top-left
		vr.x0 = REAL(0.0);  vr.y0 = H;
		vr.xa = a;           vr.yb = H - b;
		vr.xA = W;           vr.yB = REAL(0.0);
	} else {                    // top-right
		vr.x0 = W;           vr.y0 = H;
		vr.xa = W - a;       vr.yb = H - b;
		vr.xA = REAL(0.0);   vr.yB = REAL(0.0);
	}

	Real dx = vr.xA - vr.xa;
	Real dy = vr.yB - vr.yb;
	vr.D_n = sqrt(dx * dx + dy * dy);
}

void CHANNEL::Generate_VisibilityRegion()
{
	if (!g_sns_bs_enabled) {
		sns_any_limited = false;
		return;
	}

	// Step 1: Per-UT visibility probability — truncated normal clamped to [0,1]
	// Table 7.6.14.1.2-1: μ, σ per scenario
	Real Pr_raw = normal(g_sns_mu_P_vis, g_sns_sigma_P_vis);
	sns_Pr_sns = MAX(REAL(0.0), MIN(REAL(1.0), Pr_raw));

	// Step 2: Antenna array dimensions in wavelengths (including panel groups)
	Real W = (BS_N - 1) * BS_dH + (BS_Ng - 1) * (BS_N * BS_dH + BS_dgH);
	Real H = (BS_M - 1) * BS_dV + (BS_Mg - 1) * (BS_M * BS_dV + BS_dgV);

	// Step 3: Find max cluster power (dB)
	// Per spec (Eq. 7.6-57): "In the case of LOS condition, the LOS path is considered
	// as an additional cluster, and the power ratio of the LOS path to NLOS clusters
	// follows the Ricean K-factor generated in Clause 7.5."
	// Use powerForAngles[n] which includes K-factor scaling (Eq. 7.5-8):
	//   LOS: powerForAngles[0] = power[0]/(K+1) + K/(K+1)
	//   NLOS: powerForAngles[n] = power[n]/(K+1)
	// For NLOS case, powerForAngles[n] = power[n] (no K-factor).
	Real P_max_lin = REAL(0.0);
	for (int n = 0; n < num_path; n++) {
		if (powerForAngles[n] > P_max_lin) P_max_lin = powerForAngles[n];
	}
	if (P_max_lin <= REAL(0.0)) {
		sns_any_limited = false;
		return;
	}
	Real P_max_dB = REAL(10.0) * log10(P_max_lin);

	sns_any_limited = false;

	// NLOS clusters
	for (int n = 0; n < num_path; n++) {
		if (powerForAngles[n] <= REAL(0.0)) {
			sns_vr[n].limited = false;
			continue;
		}

		// Step 4: random draw to determine if cluster is visibility-limited
		Real x_n = randnum.u();
		if (x_n < sns_Pr_sns && W > REAL(0.0) && H > REAL(0.0)) {
			// Step 5: This cluster is visibility-limited
			sns_vr[n].limited = true;

			Real P_n_dB = REAL(10.0) * log10(powerForAngles[n]);
			Real delta = normal(REAL(0.0), g_sns_vr_delta);  // N(0, σ_δ)

			// V_n = A * exp(-(Pmax - Pn)/R) + B + δ, clamped to [0,1]  (Eq. 7.6-57)
			Real V_n = g_sns_vr_A * exp(-(P_max_dB - P_n_dB) / g_sns_vr_R) + g_sns_vr_B + delta;
			V_n = MAX(REAL(0.0), MIN(REAL(1.0), V_n));
			sns_vr[n].V_n = V_n;

			// a ~ unif(V_n*W, W), b = V_n*H*W/a (constant area = V_n*W*H)
			Real a_min = V_n * W;
			Real a = a_min + randnum.u() * (W - a_min);
			Real b = (a > REAL(0.0)) ? (V_n * H * W / a) : H;
			if (b > H) { b = H; a = (H > REAL(0.0)) ? (V_n * H * W / b) : W; }

			sns_setup_corner(sns_vr[n], a, b, W, H);
			sns_any_limited = true;
		} else {
			sns_vr[n].limited = false;
		}
	}

	// Zero out unused clusters
	for (int n = num_path; n < MAX_NUM_CLUSTERS; n++) {
		sns_vr[n].limited = false;
	}

	// LOS path: shares VR with cluster 0 (first cluster).
	// powerForAngles[0] = power[0]/(K+1) + K/(K+1) already includes LOS power,
	// so sns_vr[0]'s VR already reflects the LOS contribution.
	// The LOS path and cluster 0 share the same spatial direction,
	// so they should have the same visibility region.
	if (Propagation == LOS_propagation) {
		sns_vr_los = sns_vr[0];
	} else {
		sns_vr_los.limited = false;
	}
}

// Forward declaration (defined later in this file, before GetChannelImpulseResponse)
static inline Real compute_sns_attenuation(
	Real pos_h, Real pos_v, const ClusterVR& vr, Real rolloff_C);

// ====================================================================
// Compute average SNS power attenuation across all BS elements
// for each cluster. Used by Get_RSRP() to reflect SNS in coupling loss.
// Must be called after Generate_VisibilityRegion().
// ====================================================================
void CHANNEL::Compute_SNS_RSRP_Attenuation()
{
	// Default: no attenuation
	for (int n = 0; n < MAX_NUM_CLUSTERS; n++)
		sns_rsrp_power_atten[n] = REAL(1.0);
	sns_rsrp_power_atten_los = REAL(1.0);

	if (!sns_any_limited) return;

	int M = BS_M;
	int N = BS_N;
	Real total_elements = (Real)(M * N);

	// NLOS clusters
	for (int n = 0; n < num_path; n++) {
		if (!sns_vr[n].limited) {
			sns_rsrp_power_atten[n] = REAL(1.0);
			continue;
		}
		Real sum_gamma = REAL(0.0);
		for (int s_m = 0; s_m < M; s_m++) {
			Real pos_v = s_m * BS_dV;
			for (int s_n = 0; s_n < N; s_n++) {
				Real pos_h = s_n * BS_dH;
				Real atten = compute_sns_attenuation(pos_h, pos_v, sns_vr[n], g_sns_rolloff_C);
				sum_gamma += atten * atten;  // power domain (atten is amplitude)
			}
		}
		sns_rsrp_power_atten[n] = sum_gamma / total_elements;
	}

	// LOS path
	if (sns_vr_los.limited) {
		Real sum_gamma = REAL(0.0);
		for (int s_m = 0; s_m < M; s_m++) {
			Real pos_v = s_m * BS_dV;
			for (int s_n = 0; s_n < N; s_n++) {
				Real pos_h = s_n * BS_dH;
				Real atten = compute_sns_attenuation(pos_h, pos_v, sns_vr_los, g_sns_rolloff_C);
				sum_gamma += atten * atten;
			}
		}
		sns_rsrp_power_atten_los = sum_gamma / total_elements;
	}
}

void CHANNEL::Find_Strong2Clusters()
{
	/////////////////////////////////////////////////
	///
	/// Small Scale parameter --  find the 2 strongest clusters
	/// Uses powerForAngles[] (K-factor applied) for angle generation reference
	///
	/////////////////////////////////////////////////

	strongest_power = powerForAngles[0];
	strongest_power2 = powerForAngles[1];
	strongest_power_idx = 0;
	strongest_power_idx2 = 1;

	if (powerForAngles[0] < powerForAngles[1])
	{
		strongest_power = powerForAngles[1];
		strongest_power2 = powerForAngles[0];
		strongest_power_idx = 1;
		strongest_power_idx2 = 0;
	}

	for (int i = 2; i < num_path; i++)
	{
		if (powerForAngles[i] > strongest_power2)
		{
			if (powerForAngles[i] > strongest_power)
			{
				strongest_power2 = strongest_power;
				strongest_power_idx2 = strongest_power_idx;
				strongest_power = powerForAngles[i];
				strongest_power_idx = i;
			}
			else
			{
				strongest_power2 = powerForAngles[i];
				strongest_power_idx2 = i;
			}
		}
	}
}

void CHANNEL::Set_AOAAOD(int _bs_idx, int _ue_idx)
{
	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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
		// TR 38.901 Step 7 (Eq. 7.5-11/-16): X_n is an INDEPENDENT uniform{+1,-1}
		// per angle type. A single shared X rigidly mirror-couples departure and
		// arrival geometry (unphysical AOD/AOA correlation).
		Real X_AOD = (randnum.u() > 0.5) ? (Real)1 : (Real)-1;
		Real X_AOA = (randnum.u() > 0.5) ? (Real)1 : (Real)-1;

		Real Y_AOA = normal(0, ASA / 7.);
		Real Y_AOD = normal(0, ASD / 7.);

		// Real Y_AOA = normal(0, (ASA / 7. , 2.));
		// Real Y_AOD = normal(0, (ASD / 7. , 2.));

		if (i == 0)
		{
			X1_AOD = X_AOD;
			Y1_AOD = Y_AOD;

			X1_AOA = X_AOA;
			Y1_AOA = Y_AOA;
		}

		if (powerForAngles[i] == 0.)
		{
			AOD[i] = 0.;
			AOA[i] = 0.;
		}
		else
		{
			if ((TYPE == 11) &&  ((Channel_Model_Type == 0) || (Channel_Model_Type == 2)) && (carrier_freq >= 500000000 && carrier_freq <= 6000000000)) //// 0.5GHz <= fc <= 6GHz)   //// Laplacian (InH_A && F<=6GHz)
			{
				AOD[i] = -1 * ASD * log(powerForAngles[i] / strongest_power) / C;
				AOA[i] = -1 * ASA * log(powerForAngles[i] / strongest_power) / C;
			}
			else ///// Gaussian
			{
				AOD[i] = 2. * (ASD / 1.4) * sqrt(-1 * log(powerForAngles[i] / strongest_power)) / C;
				AOA[i] = 2. * (ASA / 1.4) * sqrt(-1 * log(powerForAngles[i] / strongest_power)) / C;
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

			AOD[i] = X_AOD * AOD[i] + Y_AOD - (X1_AOD * first_AOD + Y1_AOD - LOS_AOD_GCS);
			AOA[i] = X_AOA * AOA[i] + Y_AOA - (X1_AOA * first_AOA + Y1_AOA - LOS_AOA_GCS);
		}
		else // NLOS, O2I
		{
			AOD[i] = X_AOD * AOD[i] + Y_AOD + LOS_AOD_GCS;
			AOA[i] = X_AOA * AOA[i] + Y_AOA + LOS_AOA_GCS;
		}
	}
}

void CHANNEL::Set_ZOAZOD(int _bs_idx, int _ue_idx)
{
	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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
		// TR 38.901 Step 7 (Eq. 7.5-16 analog for zenith): independent uniform{+1,-1}
		// per angle type — a shared X mirror-coupled ZOD and ZOA (see Set_AOAAOD).
		Real X_ZOD = (randnum.u() > 0.5) ? (Real)1 : (Real)-1;
		Real X_ZOA = (randnum.u() > 0.5) ? (Real)1 : (Real)-1;

		Real Y_ZOA = normal(0, ZSA / 7.); //
		Real Y_ZOD = normal(0, ZSD / 7.);

		if (i == 0)
		{
			X1_ZOD = X_ZOD;
			Y1_ZOD = Y_ZOD;

			X1_ZOA = X_ZOA;
			Y1_ZOA = Y_ZOA;
		}

		if (powerForAngles[i] == 0.)
		{
			ZOD[i] = 0.;
			ZOA[i] = 0.;
		}
		else
		{
			ZOD[i] = -1 * ZSD * log(powerForAngles[i] / strongest_power) / C;
			ZOA[i] = -1 * ZSA * log(powerForAngles[i] / strongest_power) / C;
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

			ZOA[i] = (X_ZOA * ZOA[i] + Y_ZOA) - (X1_ZOA * first_ZOA + Y1_ZOA - LOS_ZOA_GCS);
			ZOD[i] = (X_ZOD * ZOD[i] + Y_ZOD) - (X1_ZOD * first_ZOD + Y1_ZOD - LOS_ZOD_GCS);
		}
		else // NLOS, O2I
		{
			if (Propagation == 2)
			{
				ZOA[i] = X_ZOA * ZOA[i] + Y_ZOA + 90.;
				ZOD[i] = X_ZOD * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
			}

			else if (Propagation == 0)
			{
				ZOA[i] = X_ZOA * ZOA[i] + Y_ZOA + LOS_ZOA_GCS;
				ZOD[i] = X_ZOD * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
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
	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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
			offset_angle_perm_zoa[path_idx][ray_idx] = 0.;
			offset_angle_perm_zod[path_idx][ray_idx] = 0.;
		}
	}

	NUM_PATH_for_channelcoeff = 0; //// find SSP per time, so initialize

	//////////////////////////////////////////////
	/////
	///// Set_subcluster
	////  strong 2 cluster -> 
	////  divide in 3 cluster and add offset, -> 
	////  In my code, not divide -> 
	////  add 4 cluster... the number of ray is different
	////
	/////////////////////////////////////////////

	// modified by jhnoh (multiplying 1e-9 for making nanosec)
	if (cluster_DS < 0) ///// cluster_DS == N/A
	{
		delay[num_path]     = delay[strongest_power_idx]  + 1.28 * 3.91 * 1e-9; 
		delay[num_path + 1] = delay[strongest_power_idx]  + 2.56 * 3.91 * 1e-9;
		delay[num_path + 2] = delay[strongest_power_idx2] + 1.28 * 3.91 * 1e-9;
		delay[num_path + 3] = delay[strongest_power_idx2] + 2.56 * 3.91 * 1e-9;
	}
	else
	{
		delay[num_path]     = delay[strongest_power_idx]  + 1.28 * cluster_DS * 1e-9;
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

	// SNS: Subclusters inherit parent cluster's Visibility Region
	// Subclusters are split from the 2 strongest clusters, so they share
	// the same spatial visibility characteristics.
	sns_vr[num_path]     = sns_vr[strongest_power_idx];
	sns_vr[num_path + 1] = sns_vr[strongest_power_idx];
	sns_vr[num_path + 2] = sns_vr[strongest_power_idx2];
	sns_vr[num_path + 3] = sns_vr[strongest_power_idx2];

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
		if (path_idx == strongest_power_idx || path_idx == strongest_power_idx2) 
		/// Subcluster No.1 (1,2,3,4,5,6,7,8,19,20 rays)
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

		else if (path_idx == num_path || path_idx == num_path + 2) 
		/// Subcluster No.2 (9,10,11,12,17,18 rays)
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

		else if (path_idx == num_path + 1 || path_idx == num_path + 3) 
		/// Subcluster No.3 (13,14,15,16 rays)
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
	// +4: strongest 2 clusters each split into 3 sub-clusters (Step 11), 
	// adding 2 extra entries per strong cluster

	// TR 38.901 Step 8 random ray coupling: relative to the AOA (identity) ray
	// ordering, AOD, ZOA and ZOD each need an INDEPENDENT random permutation of
	// the ray offsets. The fills above produced only one shuffled set (used for
	// both tx-azimuth and tx-zenith) and left the rx pair identically ordered —
	// regenerate all three here with independent Fisher-Yates draws (serial SSP
	// phase: the global randnum stream keeps this deterministic).
	for (int path_idx = 0; path_idx < NUM_PATH_for_channelcoeff; path_idx++)
	{
		int nr = (int)NUM_RAY_per_ClusterNUM[path_idx];
		if (nr <= 0) continue;
		int idx_aod[MAX_NUM_RAYS], idx_zoa[MAX_NUM_RAYS], idx_zod[MAX_NUM_RAYS];
		for (int k = 0; k < nr; k++) { idx_aod[k] = k; idx_zoa[k] = k; idx_zod[k] = k; }
		for (int k = nr - 1; k > 0; k--)
		{
			int j = (int)(randnum.u() * (k + 1)); if (j > k) j = k;
			int tmp = idx_aod[k]; idx_aod[k] = idx_aod[j]; idx_aod[j] = tmp;
		}
		for (int k = nr - 1; k > 0; k--)
		{
			int j = (int)(randnum.u() * (k + 1)); if (j > k) j = k;
			int tmp = idx_zoa[k]; idx_zoa[k] = idx_zoa[j]; idx_zoa[j] = tmp;
		}
		for (int k = nr - 1; k > 0; k--)
		{
			int j = (int)(randnum.u() * (k + 1)); if (j > k) j = k;
			int tmp = idx_zod[k]; idx_zod[k] = idx_zod[j]; idx_zod[j] = tmp;
		}
		for (int k = 0; k < nr; k++)
		{
			offset_angle_rand_coupling[path_idx][k] = offset_angle[path_idx][idx_aod[k]];
			offset_angle_perm_zoa[path_idx][k]      = offset_angle[path_idx][idx_zoa[k]];
			offset_angle_perm_zod[path_idx][k]      = offset_angle[path_idx][idx_zod[k]];
		}
	}
}

void CHANNEL::Precompute_ray_angles()
{
	for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
	{
		int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
		for (int j = 0; j < num_rays; j++)
		{
			// RX (Arrival) angles — AOA keeps the identity ray ordering; ZOA uses
			// its own independent permutation (TR 38.901 Step 8 random coupling).
			Real rx_azimuth = AOA[i] + cluster_ASA * offset_angle[i][j];
			Real rx_zenith  = ZOA[i] + cluster_ZSA * offset_angle_perm_zoa[i][j];

			// TX (Departure) angles — independent permutations for AOD and ZOD.
			Real tx_azimuth = AOD[i] + cluster_ASD * offset_angle_rand_coupling[i][j];
			Real tx_zenith  = ZOD[i] + (3.0 / 8.0) * pow(10, mu_ZSD) * offset_angle_perm_zod[i][j];

			// Transform to proper ranges
			rx_azimuth = Transform_angle_minus_180_to_plus_180(rx_azimuth);
			rx_zenith  = Transform_angle_0_to_plus_180(rx_zenith);
			tx_azimuth = Transform_angle_minus_180_to_plus_180(tx_azimuth);
			tx_zenith  = Transform_angle_0_to_plus_180(tx_zenith);

			// Store
			ray_AOA[i][j][0] = rx_azimuth;  // azimuth
			ray_AOA[i][j][1] = rx_zenith;   // zenith
			ray_AOD[i][j][0] = tx_azimuth;  // azimuth
			ray_AOD[i][j][1] = tx_zenith;   // zenith
		}
	}
	ray_angles_precomputed = true;
}

void CHANNEL::Set_RMS_delay_spread()
{
	// Power-weighted RMS delay spread from all clusters + sub-clusters
	// Note: delay[i] already contains LOS-scaled delays for LOS case (line 3193)
	// NUM_PATH_for_channelcoeff = num_path + 4 (includes sub-cluster splits)
	int N = NUM_PATH_for_channelcoeff;
	if (N <= 0) { RMS_delay_spread = 0.0; return; }

	// Power-weighted mean delay
	Real power_sum = 0.0;
	Real delay_mean = 0.0;
	for (int i = 0; i < N; i++)
	{
		delay_mean += delay[i] * power[i];
		power_sum += power[i];
	}
	if (power_sum <= 0.0) { RMS_delay_spread = 0.0; return; }
	delay_mean /= power_sum;

	// Power-weighted RMS
	Real rms_sum = 0.0;
	for (int i = 0; i < N; i++)
	{
		Real diff = delay[i] - delay_mean;
		rms_sum += diff * diff * power[i];
	}
	RMS_delay_spread = sqrt(rms_sum / power_sum);
}

// Static helper: 3GPP TR 25.996 Annex A circular angular spread with Delta optimization
static Real compute_circular_angle_spread(Real* angles, Real* powers, int N)
{
	Real power_sum = 0.0;
	for (int i = 0; i < N; i++) power_sum += powers[i];
	if (power_sum <= 0.0) return 0.0;

	Real min_spread = 1e30;
	for (int delta = 0; delta < 360; delta++)
	{
		Real d = (Real)delta;
		// Eq. A-6: weighted mean of shifted angles
		Real mu = 0.0;
		for (int i = 0; i < N; i++)
		{
			Real shifted = fmod(angles[i] + d + 180.0, 360.0);
			if (shifted < 0.0) shifted += 360.0;
			shifted -= 180.0;
			mu += shifted * powers[i];
		}
		mu /= power_sum;

		// Eq. A-4,5: spread
		Real variance = 0.0;
		for (int i = 0; i < N; i++)
		{
			Real shifted = fmod(angles[i] + d + 180.0, 360.0);
			if (shifted < 0.0) shifted += 360.0;
			shifted -= 180.0;
			Real dev = fmod(shifted - mu + 180.0, 360.0);
			if (dev < 0.0) dev += 360.0;
			dev -= 180.0;
			variance += dev * dev * powers[i];
		}
		Real spread = sqrt(variance / power_sum);
		if (spread < min_spread) min_spread = spread;
	}
	return min_spread;
}

void CHANNEL::Set_circular_angle_spread()
{
	// Count total rays across all clusters/subclusters
	int total_rays = 0;
	for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
		total_rays += (int)NUM_RAY_per_ClusterNUM[i];

	if (total_rays <= 0) return;

	// Flatten per-ray angles and powers into temporary arrays
	Real* flat_aoa = new Real[total_rays];
	Real* flat_aod = new Real[total_rays];
	Real* flat_zoa = new Real[total_rays];
	Real* flat_zod = new Real[total_rays];
	Real* flat_pow = new Real[total_rays];

	int idx = 0;
	for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
	{
		int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
		Real ray_power = power[i] / num_rays;
		for (int j = 0; j < num_rays; j++)
		{
			flat_aoa[idx] = ray_AOA[i][j][0];  // azimuth
			flat_zoa[idx] = ray_AOA[i][j][1];  // zenith
			flat_aod[idx] = ray_AOD[i][j][0];  // azimuth
			flat_zod[idx] = ray_AOD[i][j][1];  // zenith
			flat_pow[idx] = ray_power;
			idx++;
		}
	}

	circular_angle_spread_AOA = compute_circular_angle_spread(flat_aoa, flat_pow, total_rays);
	circular_angle_spread_AOD = compute_circular_angle_spread(flat_aod, flat_pow, total_rays);
	circular_angle_spread_ZOA = compute_circular_angle_spread(flat_zoa, flat_pow, total_rays);
	circular_angle_spread_ZOD = compute_circular_angle_spread(flat_zod, flat_pow, total_rays);

	delete[] flat_aoa;
	delete[] flat_aod;
	delete[] flat_zoa;
	delete[] flat_zod;
	delete[] flat_pow;
}

#if 0  // unused functions: Load_precalculate, Update, Update_per_time, Update_per_time_precise — replaced by Update_v2/Update_per_time_v2
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

	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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
			for (int p = 0; p < MS_P; p++)
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

				RX_LOS_gain_theta[m_idx * MS_Np * MS_P + n_idx * MS_P + p] = F_theta_temp;
				RX_LOS_gain_pi[m_idx * MS_Np * MS_P + n_idx * MS_P + p] = F_pi_temp;
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
					for (int p = 0; p < MS_P; p++)
					{
						Real UE_NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].ray_AOA[i][j][0] * (pi / 180.);
						Real UE_NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ray_AOA[i][j][1] * (pi / 180.);

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

						RX_NLOS_gain_theta[m_idx * MS_Np * MS_P + n_idx * MS_P + p][i][j] = F_theta_temp;
						RX_NLOS_gain_pi[m_idx * MS_Np * MS_P + n_idx * MS_P + p][i][j] = F_pi_temp;
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
	for (int p = 0; p < BS_P; p++)
	{
		for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
		{
			for (int n_idx = 0; n_idx < BS_Np; n_idx++)
			{
				r_tx.x = sin(v_angle_theta) * cos(h_angle_pi);
				r_tx.y = sin(v_angle_theta) * sin(h_angle_pi);
				r_tx.z = cos(v_angle_theta);

				weight = {0, 0};

				for (int k = 0; k < K; k++) // vertical element
				{
					for (int l = 0; l < L; l++) // horizental element
					{
						w = virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].sector_z]
													[links[_ms_idx].analog_beam_selection[adj_sector].sector_a][k][l];
						d_tx = bs[_bs_idx].d_tx[sector_num_idx][m_idx * K + k][n_idx * L + l][p][0][0]; //// sector(rx) [M][N][P][Mg][Ng]
						weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_tx, d_tx));
					}
				}
				F_theta_temp = (F_theta_GCS[p] * weight);
				F_pi_temp = (F_pi_GCS[p] * weight);

				TX_LOS_gain_theta[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_theta_temp;
				TX_LOS_gain_pi[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_pi_temp;
			}
		}
	}

	//////////////////////////////////// TX NLOS

	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
		{
			Real NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].ray_AOD[i][j][0] * (pi / 180.);
			Real NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ray_AOD[i][j][1] * (pi / 180.);

			Real NLOS_combined_antenna_gain = Get_BS_antenna_pattern(NLOS_v_angle_theta, NLOS_h_angle_pi, _bs_idx, sector_num_idx, NLOS_F_theta_GCS[0], NLOS_F_pi_GCS[0], NLOS_F_theta_GCS[1], NLOS_F_pi_GCS[1]);

			for (int p = 0; p < BS_P; p++)
			{
				for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
				{
					for (int n_idx = 0; n_idx < BS_Np; n_idx++)
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
	// NOTE: power[] no longer includes K-factor.
	// v0 Update stores raw NLOS CHIR (without K-factor scaling).
	// K-factor scaling is applied in Update_per_time (v1) for cluster 0 only,
	// and for all clusters in Update_per_time_v2 via CHIR_init with nlos_scale.
	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			complex<Real> Big_PI_LOS(0., channel[_bs_idx][_ms_idx].random_phase_vv_LOS * (pi / 180.0));   // use stored initial phase

			complex<Real> alpha_zero_temp = (RX_LOS_gain_theta[rp] * exp(Big_PI_LOS) * TX_LOS_gain_theta[tp] - RX_LOS_gain_pi[rp] * exp(Big_PI_LOS) * TX_LOS_gain_pi[tp]);

			CHIR_LOS[sector_num_idx][tp][rp] = alpha_zero_temp;

			for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
			{
				complex<Real> alpha_nmup(0, 0);
				for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
				{
					Real kappa = channel[_bs_idx][_ms_idx].kappa[i][j];
					Real _1_over_sqrt_K = 1.0 / sqrt(kappa);

					complex<Real> Big_pi_NLOS_thetatheta(0, channel[_bs_idx][_ms_idx].random_phase_vv[i][j] * (pi / 180.0));   // use stored initial phase
					complex<Real> Big_pi_NLOS_thetapi(0, channel[_bs_idx][_ms_idx].random_phase_vh[i][j] * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pitheta(0, channel[_bs_idx][_ms_idx].random_phase_hv[i][j] * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pipi(0, channel[_bs_idx][_ms_idx].random_phase_hh[i][j] * (pi / 180.0));

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
	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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

		Real v1_nlos_scale = sqrt(1.0 / (K_linear + 1.0));
		Real v1_los_scale  = sqrt(K_linear / (K_linear + 1.0));

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

				// Apply K-factor scaling to ALL clusters (power[] no longer includes K-factor)
				for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++) {
					CHIR[sector_num_idx][tp][rp][i] *= v1_nlos_scale;
				}

				// Add LOS component to cluster 0
				CHIR[sector_num_idx][tp][rp][0] += v1_los_scale * CHIR_LOS[sector_num_idx][tp][rp];

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
	// Thread-safe RNG: Use thread-local instance in OpenMP parallel regions
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

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
	Real moving_azimuth   = ms[_ms_idx].moving_direction_azimuth * (pi / 180.);  // azimuth angle (phi)

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
					Real ray_zenith  = ray_AOA[i][j][1] * (pi / 180.);

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

				// Apply K-factor scaling to ALL clusters (power[] no longer includes K-factor)
				Real precise_nlos_scale = sqrt(1.0 / (K_linear + 1.0));
				Real precise_los_scale  = sqrt(K_linear / (K_linear + 1.0));
				for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++) {
					CHIR[sector_num_idx][tp][rp][i] *= precise_nlos_scale;
				}

				// Add LOS component to cluster 0
				CHIR[sector_num_idx][tp][rp][0] += precise_los_scale * CHIR_LOS[sector_num_idx][tp][rp];
			}
		}
	}
}
#endif  // unused: Load_precalculate, Update, Update_per_time, Update_per_time_precise

// ====================================================================
// v2: Channel coefficient generation with K-factor correction and
//     Doppler frequency precomputation (3GPP TR 38.901 Eq. 7.5-22/29/30)
// ====================================================================
void CHANNEL::Update_v2(Real t, int _ms_idx, int adj_sector)
{
	#ifdef ENABLE_MULTITHREADING
	Rand& randnum = get_thread_local_rng();
	#endif

	int _bs_idx;
	if (TYPE == 11 && num_Indoor_TRxP == 1)
		_bs_idx = adj_sector;
	else
		_bs_idx = (int)(adj_sector / 3);

	int sector_num_idx = adj_sector % 3;
	Allocate_CHIR_memory(sector_num_idx);

	Real F_theta_GCS[2];
	Real F_pi_GCS[2];
	Real NLOS_F_theta_GCS[2];
	Real NLOS_F_pi_GCS[2];

	Real UE_F_theta_GCS[2];
	Real UE_F_pi_GCS[2];
	Real UE_NLOS_F_theta_GCS[2];
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

	int K = BS_M / BS_Mp;
	int L = BS_N / BS_Np;
	int ue_K = MS_M / MS_Mp;
	int ue_L = MS_N / MS_Np;

	// ---- RX LOS antenna pattern (same as Update) ----
	for (int m_idx = 0; m_idx < MS_Mp; m_idx++)
	{
		for (int n_idx = 0; n_idx < MS_Np; n_idx++)
		{
			for (int p = 0; p < MS_P; p++)
			{
				Real UE_LOS_h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOA_GCS * (pi / 180.);
				Real UE_LOS_v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOA_GCS * (pi / 180.);

				Real UE_LOS_combined_antenna_gain = Get_UE_antenna_pattern(p, UE_LOS_v_angle_theta, UE_LOS_h_angle_pi, _ms_idx, 0, UE_F_theta_GCS[0], UE_F_pi_GCS[0], UE_F_theta_GCS[1], UE_F_pi_GCS[1]);

				r_rx.x = sin(UE_LOS_v_angle_theta) * cos(UE_LOS_h_angle_pi);
				r_rx.y = sin(UE_LOS_v_angle_theta) * sin(UE_LOS_h_angle_pi);
				r_rx.z = cos(UE_LOS_v_angle_theta);

				weight = {0, 0};
				for (int k = 0; k < ue_K; k++)
				{
					for (int l = 0; l < ue_L; l++)
					{
						w = ue_virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].z]
													   [links[_ms_idx].analog_beam_selection[adj_sector].a][k][l];
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

				if (std::isnan(F_theta_temp.real()) || std::isnan(F_pi_temp.real()))
				{
					cout << "NLOS_F_theta_GCS[p] =" << NLOS_F_theta_GCS[p] << endl;
					cout << "NLOS_F_pi_GCS[p] = " << NLOS_F_pi_GCS[p] << endl;
					cout << "weight = " << weight << endl; 
					cout <<"==============================" << endl;
				}				

				RX_LOS_gain_theta[m_idx * MS_Np * MS_P + n_idx * MS_P + p] = F_theta_temp;
				RX_LOS_gain_pi[m_idx * MS_Np * MS_P + n_idx * MS_P + p] = F_pi_temp;
			}
		}
	}

	// ---- RX NLOS antenna pattern (same as Update) ----
	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
		{
			for (int m_idx = 0; m_idx < MS_Mp; m_idx++)
			{
				for (int n_idx = 0; n_idx < MS_Np; n_idx++)
				{
					for (int p = 0; p < MS_P; p++)
					{
						Real UE_NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].ray_AOA[i][j][0] * (pi / 180.);
						Real UE_NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ray_AOA[i][j][1] * (pi / 180.);

						Real UE_NLOS_combined_antenna_gain = Get_UE_antenna_pattern(p, UE_NLOS_v_angle_theta, UE_NLOS_h_angle_pi, _ms_idx, 0, UE_NLOS_F_theta_GCS[0], UE_NLOS_F_pi_GCS[0], UE_NLOS_F_theta_GCS[1], UE_NLOS_F_pi_GCS[1]);

						r_rx.x = sin(UE_NLOS_v_angle_theta) * cos(UE_NLOS_h_angle_pi);
						r_rx.y = sin(UE_NLOS_v_angle_theta) * sin(UE_NLOS_h_angle_pi);
						r_rx.z = cos(UE_NLOS_v_angle_theta);

						weight = {0, 0};
						for (int k = 0; k < ue_K; k++)
						{
							for (int l = 0; l < ue_L; l++)
							{
								w = ue_virtualization_weight_wv[links[_ms_idx].analog_beam_selection[adj_sector].z]
															   [links[_ms_idx].analog_beam_selection[adj_sector].a][k][l];
								d_rx = ms[_ms_idx].d_rx[m_idx * ue_K + k][n_idx * ue_L + l][p][0][links[_ms_idx].analog_beam_selection[adj_sector].p];
								weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_rx, d_rx));
							}
						}
						F_theta_temp = (UE_NLOS_F_theta_GCS[p] * weight);
						F_pi_temp    = (UE_NLOS_F_pi_GCS[p] * weight);

						if (ue_antenna_element_gain == 0)
						{
							F_theta_temp = (UE_NLOS_F_theta_GCS[p]);
							F_pi_temp = (UE_NLOS_F_pi_GCS[p]);
						}

						if (std::isnan(F_theta_temp.real()) || std::isnan(F_pi_temp.real()))
						{
							cout << "NLOS_F_theta_GCS[p] =" << NLOS_F_theta_GCS[p] << endl;
							cout << "NLOS_F_pi_GCS[p] = " << NLOS_F_pi_GCS[p] << endl;
							cout << "weight = " << weight << endl; 
							cout <<"==============================" << endl;
						}						

						RX_NLOS_gain_theta[m_idx * MS_Np * MS_P + n_idx * MS_P + p][i][j] = F_theta_temp;
						RX_NLOS_gain_pi[m_idx * MS_Np * MS_P + n_idx * MS_P + p][i][j] = F_pi_temp;
					}
				}
			}
		}
	}

	// ---- TX LOS antenna pattern (same as Update) ----
	Real h_angle_pi = channel[_bs_idx][_ms_idx].LOS_AOD_GCS * (pi / 180.);
	Real v_angle_theta = channel[_bs_idx][_ms_idx].LOS_ZOD_GCS * (pi / 180.);

	Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_num_idx, F_theta_GCS[0], F_pi_GCS[0], F_theta_GCS[1], F_pi_GCS[1]);

	for (int p = 0; p < BS_P; p++)
	{
		for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
		{
			// Row-beam override kept consistent with the live W_tx build in
			// Fourier_Transform_WithBF (Update_v2 currently has no callers).
			int rb_z = row_beam_enable ? sector[adj_sector].row_beam_z[m_idx]
			                           : links[_ms_idx].analog_beam_selection[adj_sector].sector_z;
			for (int n_idx = 0; n_idx < BS_Np; n_idx++)
			{
				int rb_a = links[_ms_idx].analog_beam_selection[adj_sector].sector_a;
				if (row_beam_enable && row_beam_az_mode == 1) rb_a = row_beam_boresight_a;
				if (row_beam_enable && row_beam_az_mode == 2) rb_a = sector[adj_sector].col_beam_a[n_idx];

				r_tx.x = sin(v_angle_theta) * cos(h_angle_pi);
				r_tx.y = sin(v_angle_theta) * sin(h_angle_pi);
				r_tx.z = cos(v_angle_theta);

				weight = {0, 0};
				for (int k = 0; k < K; k++)
				{
					for (int l = 0; l < L; l++)
					{
						w = virtualization_weight_wv[rb_z][rb_a][k][l];
						d_tx = bs[_bs_idx].d_tx[sector_num_idx][m_idx * K + k][n_idx * L + l][p][0][0];
						weight += w * exp(jay * (Real)2.0 * pi / Wavelength * dot(r_tx, d_tx));
					}
				}
				F_theta_temp = (F_theta_GCS[p] * weight);
				F_pi_temp = (F_pi_GCS[p] * weight);

				TX_LOS_gain_theta[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_theta_temp;
				TX_LOS_gain_pi[p * BS_Mp * BS_Np + m_idx * BS_Np + n_idx] = F_pi_temp;
			}
		}
	}

	// ---- TX NLOS antenna pattern ----
	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
		{
			Real NLOS_h_angle_pi = channel[_bs_idx][_ms_idx].ray_AOD[i][j][0] * (pi / 180.);
			Real NLOS_v_angle_theta = channel[_bs_idx][_ms_idx].ray_AOD[i][j][1] * (pi / 180.);

			Real NLOS_combined_antenna_gain = Get_BS_antenna_pattern(NLOS_v_angle_theta, NLOS_h_angle_pi, _bs_idx, sector_num_idx, NLOS_F_theta_GCS[0], NLOS_F_pi_GCS[0], NLOS_F_theta_GCS[1], NLOS_F_pi_GCS[1]);

			for (int p = 0; p < BS_P; p++)
			{
				for (int m_idx = 0; m_idx < BS_Mp; m_idx++)
				{
					// Row-beam override kept consistent with the live W_tx build in
					// Fourier_Transform_WithBF (Update_v2 currently has no callers).
					int rb_z = row_beam_enable ? sector[adj_sector].row_beam_z[m_idx]
					                           : links[_ms_idx].analog_beam_selection[adj_sector].sector_z;
					for (int n_idx = 0; n_idx < BS_Np; n_idx++)
					{
						int rb_a = links[_ms_idx].analog_beam_selection[adj_sector].sector_a;
						if (row_beam_enable && row_beam_az_mode == 1) rb_a = row_beam_boresight_a;
						if (row_beam_enable && row_beam_az_mode == 2) rb_a = sector[adj_sector].col_beam_a[n_idx];

						r_tx.x = sin(NLOS_v_angle_theta) * cos(NLOS_h_angle_pi);
						r_tx.y = sin(NLOS_v_angle_theta) * sin(NLOS_h_angle_pi);
						r_tx.z = cos(NLOS_v_angle_theta);

						weight = {0, 0};
						for (int k = 0; k < K; k++)
						{
							for (int l = 0; l < L; l++)
							{
								w = virtualization_weight_wv[rb_z][rb_a][k][l];
								d_tx = bs[_bs_idx].d_tx[sector_num_idx][m_idx * K + k][n_idx * L + l][p][0][0];
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

	// ====================================================================
	// v2 CHIR assembly per 3GPP TR 38.901 Eq. 7.5-22/29/30
	// NOTE: power[i] from Set_POWER() is normalized WITHOUT K-factor.
	//       For LOS case, NLOS clusters are scaled by sqrt(1/(K+1)) * sqrt(P_n/M_n)
	//       and LOS component is scaled by sqrt(K/(K+1)).
	// ====================================================================
	Real K_R = channel[_bs_idx][_ms_idx].K_linear;  // Ricean K-factor (linear)
	Real nlos_scale = (Propagation == LOS_propagation) ? sqrt(1.0 / (K_R + 1.0)) : 1.0;
	Real los_scale  = (Propagation == LOS_propagation) ? sqrt(K_R / (K_R + 1.0))  : 0.0;

	static int nan_report_count = 0;  // 전역 NaN 리포트 횟수 제한

	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			// LOS component — Eq. 7.5-29: [1,0;0,-1] polarization matrix
			// Scaled by sqrt(K/(K+1)) since H_LOS has no cluster power
			complex<Real> Big_PI_LOS(0., channel[_bs_idx][_ms_idx].random_phase_vv_LOS * (pi / 180.0));
			complex<Real> alpha_LOS = los_scale *
				(RX_LOS_gain_theta[rp] * exp(Big_PI_LOS) * TX_LOS_gain_theta[tp]
			   - RX_LOS_gain_pi[rp]    * exp(Big_PI_LOS) * TX_LOS_gain_pi[tp]);

			CHIR_LOS[sector_num_idx][tp][rp] = alpha_LOS;

			// NLOS components — Eq. 7.5-22
			// power[i] is normalized WITHOUT K-factor; nlos_scale applies sqrt(1/(K+1)) for LOS
			for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
			{
				Real cluster_power = channel[_bs_idx][_ms_idx].power[i];

				complex<Real> alpha_nmup(0, 0);
				for (int j = 0; j < channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]; j++)
				{
					Real kappa = channel[_bs_idx][_ms_idx].kappa[i][j];
					Real _1_over_sqrt_K = 1.0 / sqrt(kappa);

					complex<Real> Big_pi_NLOS_thetatheta(0, channel[_bs_idx][_ms_idx].random_phase_vv[i][j] * (pi / 180.0));
					complex<Real> Big_pi_NLOS_thetapi(0, channel[_bs_idx][_ms_idx].random_phase_vh[i][j] * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pitheta(0, channel[_bs_idx][_ms_idx].random_phase_hv[i][j] * (pi / 180.0));
					complex<Real> Big_pi_NLOS_pipi(0, channel[_bs_idx][_ms_idx].random_phase_hh[i][j] * (pi / 180.0));

					// nlos_scale * sqrt(P_n / M_n) — nlos_scale = sqrt(1/(K+1)) for LOS, 1.0 for NLOS
					complex<Real> alpha_nmup_temp =
					nlos_scale * sqrt(cluster_power / (channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i])) *
					( ( RX_NLOS_gain_theta[rp][i][j] * exp(Big_pi_NLOS_thetatheta) +
					    RX_NLOS_gain_pi[rp][i][j] * _1_over_sqrt_K * exp(Big_pi_NLOS_pitheta) )
					  * TX_NLOS_gain_theta[tp][i][j] +
					  ( RX_NLOS_gain_theta[rp][i][j] * _1_over_sqrt_K * exp(Big_pi_NLOS_thetapi) +
					    RX_NLOS_gain_pi[rp][i][j] * exp(Big_pi_NLOS_pipi) )
					  * TX_NLOS_gain_pi[tp][i][j] );

					// NaN 검출: ray 단위
					if (std::isnan(alpha_nmup_temp.real()) || std::isnan(alpha_nmup_temp.imag())) {
						if (nan_report_count < 10) {
							nan_report_count++;
							cout << "\n*** NaN in Update_v2 NLOS ray ***" << endl;
							cout << "  ms=" << _ms_idx << " bs=" << _bs_idx
							     << " sec=" << sector_num_idx
							     << " tp=" << tp << " rp=" << rp
							     << " cluster=" << i << " ray=" << j << endl;
							cout << "  cluster_power=" << cluster_power
							     << " power[" << i << "]=" << channel[_bs_idx][_ms_idx].power[i]
							     << " NUM_RAY=" << channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i] << endl;
							cout << "  kappa=" << kappa
							     << " 1/sqrt(kappa)=" << _1_over_sqrt_K << endl;
							cout << "  sqrt(Pn/Mn)=" << sqrt(cluster_power / channel[_bs_idx][_ms_idx].NUM_RAY_per_ClusterNUM[i]) << endl;
							cout << "  RX_theta=(" << RX_NLOS_gain_theta[rp][i][j].real() << "," << RX_NLOS_gain_theta[rp][i][j].imag() << ")"
							     << " RX_pi=(" << RX_NLOS_gain_pi[rp][i][j].real() << "," << RX_NLOS_gain_pi[rp][i][j].imag() << ")" << endl;
							cout << "  TX_theta=(" << TX_NLOS_gain_theta[tp][i][j].real() << "," << TX_NLOS_gain_theta[tp][i][j].imag() << ")"
							     << " TX_pi=(" << TX_NLOS_gain_pi[tp][i][j].real() << "," << TX_NLOS_gain_pi[tp][i][j].imag() << ")" << endl;
							cout << "  phases vv=" << channel[_bs_idx][_ms_idx].random_phase_vv[i][j]
							     << " vh=" << channel[_bs_idx][_ms_idx].random_phase_vh[i][j]
							     << " hv=" << channel[_bs_idx][_ms_idx].random_phase_hv[i][j]
							     << " hh=" << channel[_bs_idx][_ms_idx].random_phase_hh[i][j] << endl;
							cout << "  Prop=" << Propagation << " LOS=" << channel[_bs_idx][_ms_idx].LOS
							     << " K_R=" << K_R << endl;
						}
					}

					alpha_nmup += alpha_nmup_temp;
				}

				// NaN 검출: cluster 합산 후
				if (std::isnan(alpha_nmup.real()) || std::isnan(alpha_nmup.imag())) {
					if (nan_report_count < 10) {
						nan_report_count++;
						cout << "\n*** NaN in Update_v2 NLOS cluster sum ***" << endl;
						cout << "  ms=" << _ms_idx << " bs=" << _bs_idx
						     << " cluster=" << i << " tp=" << tp << " rp=" << rp
						     << " alpha_nmup=(" << alpha_nmup.real() << "," << alpha_nmup.imag() << ")" << endl;
					}
				}

				// Store NLOS-only initial value
				CHIR_init[sector_num_idx][tp][rp][i] = alpha_nmup;

				// CHIR = NLOS (Doppler added in Update_per_time_v2)
				CHIR[sector_num_idx][tp][rp][i] = alpha_nmup;
			}

			// Add LOS to cluster 0 (Eq. 7.5-30)
			if (Propagation == LOS_propagation)
			{
				CHIR[sector_num_idx][tp][rp][0] += CHIR_LOS[sector_num_idx][tp][rp];

				// NaN 검출: LOS 추가 후
				if (std::isnan(CHIR[sector_num_idx][tp][rp][0].real()) || std::isnan(CHIR[sector_num_idx][tp][rp][0].imag())) {
					if (nan_report_count < 10) {
						nan_report_count++;
						cout << "\n*** NaN in Update_v2 after LOS add ***" << endl;
						cout << "  ms=" << _ms_idx << " bs=" << _bs_idx
						     << " tp=" << tp << " rp=" << rp
						     << " CHIR_LOS=(" << CHIR_LOS[sector_num_idx][tp][rp].real() << "," << CHIR_LOS[sector_num_idx][tp][rp].imag() << ")"
						     << " los_scale=" << los_scale << " K_R=" << K_R << endl;
					}
				}
			}
		}
	}

	// ====================================================================
	// v2: Pre-compute Doppler frequencies per cluster and LOS
	// ====================================================================
	LOCATION3D v_rx;
	Real moving_elevation = ms[_ms_idx].moving_direction * (pi / 180.);
	Real moving_azimuth = ms[_ms_idx].moving_direction_azimuth * (pi / 180.);

	v_rx.x = ms[_ms_idx].speed * sin(moving_elevation) * cos(moving_azimuth);
	v_rx.y = ms[_ms_idx].speed * sin(moving_elevation) * sin(moving_azimuth);
	v_rx.z = ms[_ms_idx].speed * cos(moving_elevation);

	// NLOS cluster Doppler (using cluster center AOA/ZOA)
	for (int i = 0; i < channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff; i++)
	{
		Real aoa_rad = channel[_bs_idx][_ms_idx].AOA[i] * (pi / 180.);
		Real zoa_rad = channel[_bs_idx][_ms_idx].ZOA[i] * (pi / 180.);

		LOCATION3D r_hat;
		r_hat.x = sin(zoa_rad) * cos(aoa_rad);
		r_hat.y = sin(zoa_rad) * sin(aoa_rad);
		r_hat.z = cos(zoa_rad);

		doppler_freq_per_cluster[i] = (r_hat.x * v_rx.x + r_hat.y * v_rx.y + r_hat.z * v_rx.z) / Wavelength;
	}

	// LOS Doppler
	{
		Real los_aoa_rad = channel[_bs_idx][_ms_idx].LOS_AOA_GCS * (pi / 180.);
		Real los_zoa_rad = channel[_bs_idx][_ms_idx].LOS_ZOA_GCS * (pi / 180.);

		LOCATION3D r_hat;
		r_hat.x = sin(los_zoa_rad) * cos(los_aoa_rad);
		r_hat.y = sin(los_zoa_rad) * sin(los_aoa_rad);
		r_hat.z = cos(los_zoa_rad);

		doppler_freq_LOS_val = (r_hat.x * v_rx.x + r_hat.y * v_rx.y + r_hat.z * v_rx.z) / Wavelength;
	}
}

// ====================================================================
// v2: Time-domain update with linear-phase Doppler (no cumulative error)
// CHIR[i] = CHIR_init[i] * exp(j*2π*fd[i]*t) + LOS
// ====================================================================
void CHANNEL::Update_per_time_v2(Real t, int adj_sector, int _ms_idx)
{
	complex<Real> jay(0, 1);

	Real slot_duration = 1.0e-3 / pow(2.0, numerology);
	Real time = t * slot_duration;

	int _bs_idx;
	if (TYPE == 11 && num_Indoor_TRxP == 1)
		_bs_idx = adj_sector;
	else
		_bs_idx = (int)(adj_sector / 3);

	int sector_num_idx = adj_sector % 3;
	if (!sector_allocated[sector_num_idx])
	{
		Allocate_CHIR_memory(sector_num_idx);
	}

	int num_path = channel[_bs_idx][_ms_idx].NUM_PATH_for_channelcoeff;

	// Pre-compute Doppler phasors (outside tp/rp loops)
	complex<Real> doppler_phase[24];
	for (int i = 0; i < num_path; i++)
		doppler_phase[i] = exp(jay * (Real)2.0 * pi * doppler_freq_per_cluster[i] * time);

	complex<Real> doppler_phase_LOS = exp(jay * (Real)2.0 * pi * doppler_freq_LOS_val * time);

	for (int tp = 0; tp < NUM_TX_Port; tp++)
	{
		for (int rp = 0; rp < NUM_RX_Port; rp++)
		{
			// NLOS: initial value × Doppler (clean linear phase)
			for (int i = 0; i < num_path; i++) {
				CHIR[sector_num_idx][tp][rp][i] = CHIR_init[sector_num_idx][tp][rp][i] * doppler_phase[i];
				//if (_ms_idx == 30 )
				//{
				//	cout << "CHIR["<<sector_num_idx<<"]["<<tp<<"]["<<rp<<"]["<<i<<"] = " 
				//	<< CHIR_init[sector_num_idx][tp][rp][i] * doppler_phase[i] << endl;
				//}
			}

			// LOS: add to cluster 0 with separate Doppler
			if (Propagation == LOS_propagation)
				CHIR[sector_num_idx][tp][rp][0] += CHIR_LOS[sector_num_idx][tp][rp] * doppler_phase_LOS;
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
#if 0  // unused: Get_cos_psi, Get_sin_psi — never called
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
#endif  // unused: Get_cos_psi, Get_sin_psi

// ====================================================================
// Element-Level Channel Matrix Generation (ns-3 GetChannelImpulseResponse style)
// ====================================================================
// Generates H_usn[cluster](rxElement, txElement) without beamforming weights.
// This is equivalent to ns-3's ThreeGppChannelModel::GetChannelImpulseResponse().
// ====================================================================

void CHANNEL::Allocate_H_usn_memory(int N)
{
	// Allocate once with MAX sizes (MAX_NUM_CLUSTERS × MAX_NUM_RAYS).
	// Subsequent calls skip allocation — only update H_usn_num_clusters.
	if (H_usn_allocated) {
		H_usn_num_clusters = N;
		return;
	}

	// raysPreComp[polTx][polRx][cluster][ray] — fixed MAX sizes
	elem_raysPreComp = new ComplexReal***[BS_P];
	for (int pt = 0; pt < BS_P; pt++) {
		elem_raysPreComp[pt] = new ComplexReal**[MS_P];
		for (int pr = 0; pr < MS_P; pr++) {
			elem_raysPreComp[pt][pr] = new ComplexReal*[MAX_NUM_CLUSTERS];
			for (int n = 0; n < MAX_NUM_CLUSTERS; n++) {
				elem_raysPreComp[pt][pr][n] = new ComplexReal[MAX_NUM_RAYS]{};
			}
		}
	}

	// Spherical unit vector cache [cluster][ray] — fixed MAX sizes
	elem_sinZoA_cosAoA = new Real*[MAX_NUM_CLUSTERS];
	elem_sinZoA_sinAoA = new Real*[MAX_NUM_CLUSTERS];
	elem_cosZoA        = new Real*[MAX_NUM_CLUSTERS];
	elem_sinZoD_cosAoD = new Real*[MAX_NUM_CLUSTERS];
	elem_sinZoD_sinAoD = new Real*[MAX_NUM_CLUSTERS];
	elem_cosZoD        = new Real*[MAX_NUM_CLUSTERS];
	for (int n = 0; n < MAX_NUM_CLUSTERS; n++) {
		elem_sinZoA_cosAoA[n] = new Real[MAX_NUM_RAYS]{};
		elem_sinZoA_sinAoA[n] = new Real[MAX_NUM_RAYS]{};
		elem_cosZoA[n]        = new Real[MAX_NUM_RAYS]{};
		elem_sinZoD_cosAoD[n] = new Real[MAX_NUM_RAYS]{};
		elem_sinZoD_sinAoD[n] = new Real[MAX_NUM_RAYS]{};
		elem_cosZoD[n]        = new Real[MAX_NUM_RAYS]{};
	}

	// H_usn: vector of (totalRx × totalTx) matrices — fixed MAX_NUM_CLUSTERS
	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;
	H_usn.resize(MAX_NUM_CLUSTERS);
	for (int n = 0; n < MAX_NUM_CLUSTERS; n++) {
		H_usn[n] = MatrixXcReal::Zero(totalRx, totalTx);
	}

	H_usn_allocated = true;
	H_usn_num_clusters = N;
}

#if 0  // unused: Reset_H_usn_memory — never called
void CHANNEL::Reset_H_usn_memory()
{
	// Zero-fill all element-level buffers for reuse. No deallocation.
	if (!H_usn_allocated) return;

	int N = H_usn_num_clusters;

	// Zero-fill raysPreComp (only used range)
	for (int pt = 0; pt < BS_P; pt++) {
		for (int pr = 0; pr < MS_P; pr++) {
			for (int n = 0; n < N; n++) {
				int M = (int)NUM_RAY_per_ClusterNUM[n];
				for (int m = 0; m < M; m++) {
					elem_raysPreComp[pt][pr][n][m] = ComplexReal(REAL(0.0), REAL(0.0));
				}
			}
		}
	}

	// Zero-fill spherical caches (only used range)
	for (int n = 0; n < N; n++) {
		int M = (int)NUM_RAY_per_ClusterNUM[n];
		for (int m = 0; m < M; m++) {
			elem_sinZoA_cosAoA[n][m] = REAL(0.0);
			elem_sinZoA_sinAoA[n][m] = REAL(0.0);
			elem_cosZoA[n][m]        = REAL(0.0);
			elem_sinZoD_cosAoD[n][m] = REAL(0.0);
			elem_sinZoD_sinAoD[n][m] = REAL(0.0);
			elem_cosZoD[n][m]        = REAL(0.0);
		}
	}

	// Zero-fill H_usn matrices (only used range)
	for (int n = 0; n < N; n++) {
		H_usn[n].setZero();
	}
}
#endif  // unused: Reset_H_usn_memory

// ====================================================================
// SNS helper: compute amplitude attenuation for a TX element position
// relative to a cluster's Visibility Region boundary.
// j20 Eq. 7.6-58: γ = exp(-d / C)  [power domain]
// C = roll-off parameter (wavelength units), Table 7.6.14.1.2-3
// Returns sqrt(γ) for amplitude scaling.
// pos_h, pos_v: element position in wavelength units (LCS)
// ====================================================================
// compute_sns_attenuation() is now defined in h/channel.h (shared with Link.cpp)

void CHANNEL::GetChannelImpulseResponse(int bs_idx, int ms_idx, int sector_idx)
{
	// ====================================================================
	// ns-3 style element-level channel coefficient generation
	// Output: H_usn[cluster](rxElement, txElement) without beamforming
	// Reference: 3GPP TR 38.901 Eq. 7.5-22 (NLOS), 7.5-29/30 (LOS)
	// ====================================================================

	int N = NUM_PATH_for_channelcoeff;
	if (N <= 0) return;

	// Phase A: Allocate memory
	Allocate_H_usn_memory(N);

	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;
	Real k_2pi = REAL(2.0) * pi / Wavelength;
	complex<Real> jay(REAL(0.0), REAL(1.0));
	const Real deg2rad = pi / REAL(180.0);

	// Phase B: Cache spherical unit vectors per ray
	for (int n = 0; n < N; n++) {
		int M = (int)NUM_RAY_per_ClusterNUM[n];
		for (int m = 0; m < M; m++) {
			// RX arrival direction (AOA, ZOA) — ray_AOA[n][m][0]=azimuth, [1]=zenith
			Real aoa_rad = ray_AOA[n][m][0] * deg2rad;
			Real zoa_rad = ray_AOA[n][m][1] * deg2rad;
			elem_sinZoA_cosAoA[n][m] = sin(zoa_rad) * cos(aoa_rad);
			elem_sinZoA_sinAoA[n][m] = sin(zoa_rad) * sin(aoa_rad);
			elem_cosZoA[n][m]        = cos(zoa_rad);

			// TX departure direction (AOD, ZOD) — ray_AOD[n][m][0]=azimuth, [1]=zenith
			Real aod_rad = ray_AOD[n][m][0] * deg2rad;
			Real zod_rad = ray_AOD[n][m][1] * deg2rad;
			elem_sinZoD_cosAoD[n][m] = sin(zod_rad) * cos(aod_rad);
			elem_sinZoD_sinAoD[n][m] = sin(zod_rad) * sin(aod_rad);
			elem_cosZoD[n][m]        = cos(zod_rad);
		}
	}

	// Phase C: Pre-compute raysPreComp[polTx][polRx][cluster][ray]
	// = antenna_pattern × polarization_matrix × initial_phases
	// (independent of element position u, s)
	Real tx_F_theta_GCS[2], tx_F_pi_GCS[2];
	Real rx_F_theta_GCS[2], rx_F_pi_GCS[2];

	// For handheld: store TX × polarization intermediate per ray
	// raysPreComp = rx_Ft · txPol_theta_col + rx_Fp · txPol_phi_col
	// (RX pattern computed per-port in Phase D)
	ComplexReal txPol_theta_col[2][MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	ComplexReal txPol_phi_col[2][MAX_NUM_CLUSTERS][MAX_NUM_RAYS];

	for (int n = 0; n < N; n++) {
		int M = (int)NUM_RAY_per_ClusterNUM[n];
		for (int m = 0; m < M; m++) {
			// Initial random phases (degrees → radians)
			Real phase_vv = random_phase_vv[n][m] * deg2rad;
			Real phase_vh = random_phase_vh[n][m] * deg2rad;
			Real phase_hv = random_phase_hv[n][m] * deg2rad;
			Real phase_hh = random_phase_hh[n][m] * deg2rad;
			Real kappa_nm = kappa[n][m];
			Real inv_sqrt_kappa = REAL(1.0) / sqrt(kappa_nm);

			// TX antenna pattern (BS) — returns both polarization patterns
			Real aod_rad = ray_AOD[n][m][0] * deg2rad;
			Real zod_rad = ray_AOD[n][m][1] * deg2rad;
			Get_BS_antenna_pattern(zod_rad, aod_rad, bs_idx, sector_idx,
				tx_F_theta_GCS[0], tx_F_pi_GCS[0], tx_F_theta_GCS[1], tx_F_pi_GCS[1]);

			// Pre-compute exp(j*phase)
			ComplexReal exp_vv = exp(jay * phase_vv);
			ComplexReal exp_vh = exp(jay * phase_vh);
			ComplexReal exp_hv = exp(jay * phase_hv);
			ComplexReal exp_hh = exp(jay * phase_hh);

			if (handheld_mode) {
				// Handheld: store TX × polarization matrix intermediate
				// Per-port RX field will be combined in Phase D
				for (int polTx = 0; polTx < BS_P; polTx++) {
					Real tx_Ft = tx_F_theta_GCS[polTx];
					Real tx_Fp = tx_F_pi_GCS[polTx];
					// θ-row: [e^jΦ_θθ, √(1/κ)·e^jΦ_θφ] · [F_tx_θ; F_tx_φ]
					txPol_theta_col[polTx][n][m] = exp_vv * tx_Ft + inv_sqrt_kappa * exp_vh * tx_Fp;
					// φ-row: [√(1/κ)·e^jΦ_φθ, e^jΦ_φφ] · [F_tx_θ; F_tx_φ]
					txPol_phi_col[polTx][n][m] = inv_sqrt_kappa * exp_hv * tx_Ft + exp_hh * tx_Fp;
				}
			} else {
				// Non-handheld: pre-compute full raysPreComp (RX pattern same for all ports)
				Real aoa_rad = ray_AOA[n][m][0] * deg2rad;
				Real zoa_rad = ray_AOA[n][m][1] * deg2rad;
				if (ue_antenna_element_gain == 0) {
					rx_F_theta_GCS[0] = REAL(1.0); rx_F_pi_GCS[0] = REAL(0.0);
					rx_F_theta_GCS[1] = REAL(0.0); rx_F_pi_GCS[1] = REAL(1.0);
				} else {
					Get_UE_antenna_pattern(0, zoa_rad, aoa_rad, ms_idx, 0,
						rx_F_theta_GCS[0], rx_F_pi_GCS[0], rx_F_theta_GCS[1], rx_F_pi_GCS[1]);
				}

				for (int polTx = 0; polTx < BS_P; polTx++) {
					Real tx_Ft = tx_F_theta_GCS[polTx];
					Real tx_Fp = tx_F_pi_GCS[polTx];
					for (int polRx = 0; polRx < MS_P; polRx++) {
						Real rx_Ft = rx_F_theta_GCS[polRx];
						Real rx_Fp = rx_F_pi_GCS[polRx];
						// Eq. 7.5-22 polarization matrix
						elem_raysPreComp[polTx][polRx][n][m] =
							rx_Ft * exp_vv * tx_Ft +
							rx_Ft * inv_sqrt_kappa * exp_vh * tx_Fp +
							rx_Fp * inv_sqrt_kappa * exp_hv * tx_Ft +
							rx_Fp * exp_hh * tx_Fp;
					}
				}
			}
		}
	}

	// (SNS array dimensions now stored in ClusterVR struct)

	// Phase D: Compute H_usn[n](u, s) = Σ_m { √(P_n/M) · raysPreComp · exp(j·rxPhase) · exp(j·txPhase) }
	for (int n = 0; n < N; n++) {
		H_usn[n].setZero();
		int M = (int)NUM_RAY_per_ClusterNUM[n];
		Real norm = sqrt(power[n] / REAL(M));

		for (int u = 0; u < totalRx; u++) {
			// RX element decomposition: u = u_m*(MS_N*MS_P) + u_n*MS_P + u_p
			int u_p = u % MS_P;                    // polRx index
			int u_n = (u / MS_P) % MS_N;
			int u_m = u / (MS_N * MS_P);
			LOCATION3D rxLoc = ms[ms_idx].d_rx[u_m][u_n][u_p][0][0];

			// Handheld: pre-compute per-port RX field pattern for all rays in this cluster
			Real port_rx_Ft[MAX_NUM_RAYS], port_rx_Fp[MAX_NUM_RAYS];
			if (handheld_mode) {
				for (int m = 0; m < M; m++) {
					Real aoa_rad = ray_AOA[n][m][0] * deg2rad;
					Real zoa_rad = ray_AOA[n][m][1] * deg2rad;
					Real dummy1, dummy2;
					Get_UE_antenna_pattern(0, zoa_rad, aoa_rad, ms_idx, 0,
						port_rx_Ft[m], port_rx_Fp[m], dummy1, dummy2, u_m);
				}
			}

			for (int s = 0; s < totalTx; s++) {
				// TX element decomposition: s = s_m*(BS_N*BS_P) + s_n*BS_P + s_p
				int s_p = s % BS_P;                // polTx index
				int s_n = (s / BS_P) % BS_N;
				int s_m = s / (BS_N * BS_P);
				LOCATION3D txLoc = bs[bs_idx].d_tx[sector_idx][s_m][s_n][s_p][0][0];

				ComplexReal h_val(REAL(0.0), REAL(0.0));
				for (int m = 0; m < M; m++) {
					// RX spatial phase: 2π/λ · r̂_rx · d_rx
					Real rxPhase = k_2pi * (
						elem_sinZoA_cosAoA[n][m] * rxLoc.x +
						elem_sinZoA_sinAoA[n][m] * rxLoc.y +
						elem_cosZoA[n][m]        * rxLoc.z);

					// TX spatial phase: 2π/λ · r̂_tx · d_tx
					Real txPhase = k_2pi * (
						elem_sinZoD_cosAoD[n][m] * txLoc.x +
						elem_sinZoD_sinAoD[n][m] * txLoc.y +
						elem_cosZoD[n][m]        * txLoc.z);

					// Handheld: combine per-port RX field with TX+polarization intermediate
					// Non-handheld: use pre-computed raysPreComp
					ComplexReal rayComp;
					if (handheld_mode) {
						rayComp = port_rx_Ft[m] * txPol_theta_col[s_p][n][m] +
						          port_rx_Fp[m] * txPol_phi_col[s_p][n][m];
					} else {
						rayComp = elem_raysPreComp[s_p][u_p][n][m];
					}

					h_val += rayComp *
						ComplexReal(cos(rxPhase), sin(rxPhase)) *
						ComplexReal(cos(txPhase), sin(txPhase));
				}
				// SNS: BS-side per-element per-cluster attenuation
				// Element LCS position: horizontal = s_n * dH, vertical = s_m * dV
				// (For multi-panel, add ng*dgH / mg*dgV when Phase D supports Mg/Ng)
				Real sns_atten = REAL(1.0);
				if (sns_any_limited) {
					Real pos_h = s_n * BS_dH;
					Real pos_v = s_m * BS_dV;
					sns_atten = compute_sns_attenuation(pos_h, pos_v, sns_vr[n], g_sns_rolloff_C);
				}
				H_usn[n](u, s) = norm * sns_atten * h_val;
			}
		}
	}

	// Phase E: LOS component (Eq. 7.5-29, 7.5-30) + Store H_usn_init / H_usn_LOS for Doppler
	Allocate_H_usn_Init_memory(N);

	if (Propagation == LOS_propagation) {
		Real K_R = K_linear;

		// LOS direction
		Real los_aoa_rad = LOS_AOA_GCS * deg2rad;
		Real los_zoa_rad = LOS_ZOA_GCS * deg2rad;
		Real los_aod_rad = LOS_AOD_GCS * deg2rad;
		Real los_zod_rad = LOS_ZOD_GCS * deg2rad;

		// LOS propagation phase (ns-3: phaseDiffDueToDistance)
		ComplexReal losPhase = exp(-jay * REAL(2.0) * pi * distance / Wavelength);

		// LOS initial phase (from random_phase_vv_LOS)
		Real los_init_phase_rad = random_phase_vv_LOS * deg2rad;
		ComplexReal losInitPhase = exp(jay * los_init_phase_rad);

		// LOS antenna patterns (TX)
		Real los_tx_Ft[2], los_tx_Fp[2];
		Real los_rx_Ft[2], los_rx_Fp[2];
		Get_BS_antenna_pattern(los_zod_rad, los_aod_rad, bs_idx, sector_idx,
			los_tx_Ft[0], los_tx_Fp[0], los_tx_Ft[1], los_tx_Fp[1]);
		// RX pattern: for handheld, computed per-port in the u loop below
		if (!handheld_mode) {
			if (ue_antenna_element_gain == 0) {
				los_rx_Ft[0] = REAL(1.0); los_rx_Fp[0] = REAL(0.0);
				los_rx_Ft[1] = REAL(0.0); los_rx_Fp[1] = REAL(1.0);
			} else {
				Get_UE_antenna_pattern(0, los_zoa_rad, los_aoa_rad, ms_idx, 0,
					los_rx_Ft[0], los_rx_Fp[0], los_rx_Ft[1], los_rx_Fp[1]);
			}
		}

		// LOS unit vectors
		Real los_sinZoA_cosAoA = sin(los_zoa_rad) * cos(los_aoa_rad);
		Real los_sinZoA_sinAoA = sin(los_zoa_rad) * sin(los_aoa_rad);
		Real los_cosZoA        = cos(los_zoa_rad);
		Real los_sinZoD_cosAoD = sin(los_zod_rad) * cos(los_aod_rad);
		Real los_sinZoD_sinAoD = sin(los_zod_rad) * sin(los_aod_rad);
		Real los_cosZoD        = cos(los_zod_rad);

		Real scale_nlos = sqrt(REAL(1.0) / (K_R + REAL(1.0)));
		Real scale_los  = sqrt(K_R / (K_R + REAL(1.0)));

		// Build H_usn_LOS matrix (before applying K-factor to H_usn)
		H_usn_LOS.setZero();
		for (int u = 0; u < totalRx; u++) {
			int u_p = u % MS_P;
			int u_n = (u / MS_P) % MS_N;
			int u_m = u / (MS_N * MS_P);
			LOCATION3D rxLoc = ms[ms_idx].d_rx[u_m][u_n][u_p][0][0];

			// Handheld: compute per-port LOS RX field pattern
			Real port_los_rx_Ft, port_los_rx_Fp;
			if (handheld_mode) {
				Real dummy1, dummy2;
				Get_UE_antenna_pattern(0, los_zoa_rad, los_aoa_rad, ms_idx, 0,
					port_los_rx_Ft, port_los_rx_Fp, dummy1, dummy2, u_m);
			} else {
				port_los_rx_Ft = los_rx_Ft[u_p];
				port_los_rx_Fp = los_rx_Fp[u_p];
			}

			Real rxLosPhase = k_2pi * (
				los_sinZoA_cosAoA * rxLoc.x +
				los_sinZoA_sinAoA * rxLoc.y +
				los_cosZoA        * rxLoc.z);

			for (int s = 0; s < totalTx; s++) {
				int s_p = s % BS_P;
				int s_n = (s / BS_P) % BS_N;
				int s_m = s / (BS_N * BS_P);
				LOCATION3D txLoc = bs[bs_idx].d_tx[sector_idx][s_m][s_n][s_p][0][0];

				Real txLosPhase = k_2pi * (
					los_sinZoD_cosAoD * txLoc.x +
					los_sinZoD_sinAoD * txLoc.y +
					los_cosZoD        * txLoc.z);

				// LOS polarization (Eq. 7.5-29): [F_rx_θ·F_tx_θ − F_rx_φ·F_tx_φ]
				ComplexReal losRay =
					(port_los_rx_Ft * los_tx_Ft[s_p] - port_los_rx_Fp * los_tx_Fp[s_p]) *
					losInitPhase * losPhase *
					ComplexReal(cos(rxLosPhase), sin(rxLosPhase)) *
					ComplexReal(cos(txLosPhase), sin(txLosPhase));

				// SNS: BS-side LOS attenuation
				// (LOS path is almost never limited since it dominates total power,
				//  but included for completeness per spec Section 7.6.14.1.3)
				Real sns_los_atten = REAL(1.0);
				if (sns_any_limited) {
					Real los_pos_h = s_n * BS_dH;
					Real los_pos_v = s_m * BS_dV;
					sns_los_atten = compute_sns_attenuation(los_pos_h, los_pos_v, sns_vr_los, g_sns_rolloff_C);
				}
				H_usn_LOS(u, s) = scale_los * sns_los_atten * losRay;
			}
		}

		// Store H_usn_init = NLOS scaled by sqrt(1/(K+1))
		for (int n = 0; n < N; n++) {
			H_usn_init[n] = H_usn[n] * ComplexReal(scale_nlos, REAL(0.0));
		}

		// Combine for t=0: H_usn[0] = H_usn_init[0] + H_usn_LOS
		H_usn[0] = H_usn_init[0] + H_usn_LOS;
		for (int n = 1; n < N; n++) {
			H_usn[n] = H_usn_init[n];
		}
	}
	else {
		// NLOS/O2I: H_usn_init = H_usn (no K-factor scaling)
		for (int n = 0; n < N; n++) {
			H_usn_init[n] = H_usn[n];
		}
		H_usn_LOS.setZero();
	}

	// Phase F: Pre-compute Doppler frequencies
	PrecomputeDoppler(bs_idx, ms_idx);
}

// ====================================================================
// Allocate H_usn_init and H_usn_LOS for Doppler time evolution
// ====================================================================
void CHANNEL::Allocate_H_usn_Init_memory(int N)
{
	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;

	if (!H_usn_init_allocated) {
		H_usn_init.resize(MAX_NUM_CLUSTERS);
		for (int n = 0; n < MAX_NUM_CLUSTERS; n++) {
			H_usn_init[n] = MatrixXcReal::Zero(totalRx, totalTx);
		}
		H_usn_LOS = MatrixXcReal::Zero(totalRx, totalTx);
		H_usn_init_allocated = true;
	} else {
		// Zero-fill for reuse
		for (int n = 0; n < N; n++) {
			H_usn_init[n].setZero();
		}
		H_usn_LOS.setZero();
	}
}

// ====================================================================
// Pre-compute Doppler frequencies per cluster and LOS
// ====================================================================
void CHANNEL::PrecomputeDoppler(int bs_idx, int ms_idx)
{
	LOCATION3D v_rx;
	Real moving_elevation = ms[ms_idx].moving_direction * (pi / 180.);
	Real moving_azimuth = ms[ms_idx].moving_direction_azimuth * (pi / 180.);

	v_rx.x = ms[ms_idx].speed * sin(moving_elevation) * cos(moving_azimuth);
	v_rx.y = ms[ms_idx].speed * sin(moving_elevation) * sin(moving_azimuth);
	v_rx.z = ms[ms_idx].speed * cos(moving_elevation);

	// NLOS cluster Doppler (using cluster center AOA/ZOA)
	for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
	{
		Real aoa_rad = AOA[i] * (pi / 180.);
		Real zoa_rad = ZOA[i] * (pi / 180.);

		LOCATION3D r_hat;
		r_hat.x = sin(zoa_rad) * cos(aoa_rad);
		r_hat.y = sin(zoa_rad) * sin(aoa_rad);
		r_hat.z = cos(zoa_rad);

		doppler_freq_per_cluster[i] = (r_hat.x * v_rx.x + r_hat.y * v_rx.y + r_hat.z * v_rx.z) / Wavelength;
	}

	// LOS Doppler
	{
		Real los_aoa_rad = LOS_AOA_GCS * (pi / 180.);
		Real los_zoa_rad = LOS_ZOA_GCS * (pi / 180.);

		LOCATION3D r_hat;
		r_hat.x = sin(los_zoa_rad) * cos(los_aoa_rad);
		r_hat.y = sin(los_zoa_rad) * sin(los_aoa_rad);
		r_hat.z = cos(los_zoa_rad);

		doppler_freq_LOS_val = (r_hat.x * v_rx.x + r_hat.y * v_rx.y + r_hat.z * v_rx.z) / Wavelength;
	}
}

// ====================================================================
// Element-level time-domain update with linear-phase Doppler
// H_usn[n] = H_usn_init[n] * exp(j*2π*fd[n]*t)
// H_usn[0] += H_usn_LOS * exp(j*2π*fd_LOS*t)
// ====================================================================
void CHANNEL::Update_H_usn_per_time(Real t_idx, int ms_idx, int sector_idx)
{
	complex<Real> jay(0, 1);

	Real slot_duration = 1.0e-3 / pow(2.0, numerology);
	Real time = t_idx * slot_duration;

	int num_path = NUM_PATH_for_channelcoeff;
	if (num_path <= 0) return;
	if (!H_usn_init_allocated) return;

	// Pre-compute Doppler phasors
	complex<Real> doppler_phase[24];
	for (int i = 0; i < num_path; i++)
		doppler_phase[i] = exp(jay * (Real)2.0 * pi * doppler_freq_per_cluster[i] * time);

	complex<Real> doppler_phase_LOS = exp(jay * (Real)2.0 * pi * doppler_freq_LOS_val * time);

	// NLOS: H_usn[n] = H_usn_init[n] * doppler_phase[n]
	for (int n = 0; n < num_path; n++) {
		H_usn[n] = H_usn_init[n] * doppler_phase[n];
	}

	// LOS: add to cluster 0 with separate Doppler
	if (Propagation == LOS_propagation) {
		H_usn[0] += H_usn_LOS * doppler_phase_LOS;
	}
}