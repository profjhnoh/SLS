#include "const.h"

#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

#ifndef _CHANNEL_H_
#define _CHANNEL_H_

class CHANNEL
{
public:

	void Set_channel(Real, bool, int, int, LOCATION, LOCATION);
	void Set_SmallScaleParameter(int, int);
	void Allocate_memory();
	void Delete_memory();
	void Allocate_CHIR_memory(int sector_idx);
	void Delete_CHIR_memory();
	void Update(Real, int, int);
	void Update_per_time(Real,int, int);
	void Update_per_time_precise(Real,int, int);  // Ray-level precision mode
	void Update_v2(Real, int, int);               // v2: K-factor corrected channel coefficients
	void Update_per_time_v2(Real, int, int);      // v2: Linear-phase Doppler update
	void Precalculate(int);
	void Load_precalculate(Real, int, int);
	void Set_LOS_Prob();
	void Set_PATHLOSS();
	void Set_Channel_Parameters();
	void Set_DELAY();
	void Set_POWER();
	void Find_Strong2Clusters();
	void Set_AOAAOD(int, int);
	void Set_ZOAZOD(int, int);
	void Set_InitialPhase();
	void Set_SUBCLUSTER();
	void Set_W_matrix();
	void Set_RMS_delay_spread();
	void Set_circular_angle_spread();
	void Precompute_ray_angles();

	void GetNewChannel(int bs_idx, int ms_idx, int sector_idx);
	void Allocate_H_usn_memory(int N);
	void Allocate_H_usn_Init_memory(int N);
	void Reset_H_usn_memory();
	void Update_H_usn_per_time(Real t, int ms_idx, int sector_idx);
	void PrecomputeDoppler(int bs_idx, int ms_idx);

	void Set_AOAAOD_LOS();
	void Set_AOAAOD_NLOS();
	void Set_ZOAZOD_LOS();
	void Set_ZOAZOD_NLOS();
	//void Kill_Memory();
	void Reset2Default();
	void Sampling_DelaySpread();

	Real            _bs_height                        = 0;

	LOCATION          self_bs                           = {0,0};
	LOCATION          self_ms                           = {0,0};
	LOCATION          bs_location                       = {0,0};
	LOCATION          ms_location                       = {0,0};

	Real            LOS_prob                          = 0;
	Real            pathloss                          = 0;
	Real            pathloss_2                        = 0;
	Real            pathloss_final                    = 0;
	bool              LOS                               = false;
	int               Propagation                       = 0;
	bool              Indoor                            = false;
	Real            RMS_delay_spread                  = 0;
	Real            circular_angle_spread_AOA         = 0;
	Real            circular_angle_spread_AOD         = 0;
	Real            circular_angle_spread_ZOA         = 0;
	Real            circular_angle_spread_ZOD         = 0;
	int               self_bs_idx                       = 0;
	int               self_ms_idx                       = 0;
	Real            LOS_AOA_GCS                       = 0;
	Real            LOS_AOD_GCS                       = 0;
	Real            LOS_ZOA_GCS                       = 0;
	Real            LOS_ZOD_GCS                       = 0;
	Real            strongest_power                   = 0;
	Real            strongest_power2                  = 0;
	int               strongest_power_idx               = 0;
	int               strongest_power_idx2              = 0;
	Real            d_rx                              = 0;
	Real            d_tx                              = 0;
    int               ms_idx                            = 0;
	Real            random_phase_vv_LOS               = 0;
	Real            random_phase_hh_LOS               = 0;
	int               num_ray                           = 0;
	int               num_path                          = 0;
	int               NUM_PATH_for_channelcoeff         = 0; 
	Real            distance                          = 0;
	Real            D                                 = 0;
	Real            r_tau                             = 0;
	Real            K_factor                          = 0;
	Real            DS                                = 0; // = sigma_tau in M.2135
	Real            mu_K_factor                       = 0;
	Real            sigma_K_factor                    = 0;
	Real            mu_DS                             = 0;
	Real            sigma_DS                          = 0;
	Real            mu_ASD                            = 0;
	Real            sigma_ASD                         = 0;
	Real            mu_ASA                            = 0;
	Real            sigma_ASA                         = 0;
	Real            mu_ZSA                            = 0;
	Real            sigma_ZSA                         = 0;
	Real            mu_ZSD                            = 0;
	Real            sigma_ZSD                         = 0;
	Real            _DS                               = 0;
	Real            K_dB                              = 0;
	Real            K_linear                          = 0;
	Real            ASA                               = 0;
	Real            ASD                               = 0;
	Real            ZSA                               = 0;
	Real            ZSD                               = 0;
	Real            P_1LOS                            = 0;
	Real            mu_XPR                            = 0;
	Real            sigma_XPR                         = 0;
	Real            XPR_mean[3]                       = {0,};
	Real            XPR_std[3]                        = {0,};
	Real            kappa[24][20]                     = {0,};
	Real            mu_offset_ZOD                     = 0;
	Real            cluster_DS                        = 0;
	Real            cluster_ASD                       = 0;
	Real            cluster_ASA                       = 0;
	Real            cluster_ZSA                       = 0;
	Real            XPR                               = 0;
	Real            cluster_shadowing                 = 0;
	Real            sigma_SF                          = 0;
	Real            BS_d_H                            = 0;
	Real            BS_d_V                            = 0;
	Real            BS_d_gH                           = 0;
	Real            BS_d_gV                           = 0;
	Real            MS_d_H                            = 0;
	Real            MS_d_V                            = 0;
	Real            MS_d_gH                           = 0;
	Real            MS_d_gV                           = 0;
	Real            alpha                             = 0;
	Real            beta                              = 0;
	Real            gamma                             = 0;
	Real            correlated_randnum_SF             = 0;
	Real            correlated_randnum_DS             = 0;
	Real            correlated_randnum_ASD            = 0;
	Real            correlated_randnum_ASA            = 0;
	Real            correlated_randnum_K              = 0;
	Real            final_ant_gain                    = 0;
	Real            signal_RSRP_gain[3][4][8][2][4][2] = {0,}; //[sec_idx] [sec_zenith][sec_azimuth][ue_zenith_idx][ue_azimuth_idx][pannel];  
	
	ComplexReal  F_tx_v                            = ComplexReal (0,0);
	ComplexReal  F_tx_h                            = ComplexReal (0,0);
	ComplexReal  F_rx_v                            = ComplexReal (0,0);
	ComplexReal  F_rx_h                            = ComplexReal (0,0);
	ComplexReal  F_vv                              = ComplexReal (0,0);
	ComplexReal  F_vh                              = ComplexReal (0,0);
	ComplexReal  F_hv                              = ComplexReal (0,0);
	ComplexReal  F_hh                              = ComplexReal (0,0);

	int *             sampled_delay                     = NULL;
	Real*           delay                             = NULL;
	Real*           delay_LOS                         = NULL;
	Real*           power                             = NULL;
	Real*           powerForAngles                    = NULL;  // K-factor applied power for angle generation only (Eq. 7.5-8)
	Real*           power_LOS                         = NULL;
	Real*           power_NLOS                        = NULL;

	Real*           sync_delay                        = NULL;
	Real*           sync_delay_LOS                    = NULL;
	Real*           sync_power                        = NULL;
	Real*           sync_power_LOS                    = NULL;
	Real*           sync_power_NLOS                   = NULL;

	Real*           AOA                               = NULL;
	Real*           AOD                               = NULL;
	Real*           ZOA                               = NULL;
	Real*           ZOD                               = NULL;
	//Real*           AOA_LOS                           = NULL;
	//Real*           AOD_LOS                           = NULL;
	//Real*           ZOA_LOS                           = NULL;
	//Real*           ZOD_LOS                           = NULL;
	//Real*           AOA_NLOS                          = NULL;
	//Real*           AOD_NLOS                          = NULL;
	//Real*           ZOA_NLOS                          = NULL;
	//Real*           ZOD_NLOS                          = NULL;
	Real*           delay_power                       = NULL;
	Real*           power_ray_sum                     = NULL;

	Real*           RMS_delay_power                   = NULL;

	Real*           power_ray_sum_circular            = NULL;
	Real*           AOA_power                         = NULL;
	Real*           AOD_power                         = NULL;
	Real*           theta_n_m_mu_AOA                  = NULL;
	Real*           theta_n_m_mu_AOD                  = NULL;
	Real*           theta_power_AOA                   = NULL;
	Real*           theta_power_AOD                   = NULL;
	Real**          random_phase_vv                   = NULL;
	Real**          random_phase_vh                   = NULL;
	Real**          random_phase_hv                   = NULL;
	Real**          random_phase_hh                   = NULL;
	Real **         offset_angle                      = NULL;
	Real **         offset_angle_rand_coupling        = NULL;
	Real *          NUM_RAY_per_ClusterNUM            = NULL;

	ComplexReal **** CHIR                          = NULL;
	ComplexReal *** CHIR_LOS                       = NULL;

	// v2: Initial CHIR storage (Doppler-free, K-factor applied)
	ComplexReal **** CHIR_init                     = NULL;  // [sector][TX_port][RX_port][cluster]

	// v2: Pre-computed Doppler frequencies (Hz)
	Real doppler_freq_per_cluster[24]              = {0};   // per NLOS cluster
	Real doppler_freq_LOS_val                      = 0;     // LOS path

	// ====================================================================
	// CHIR_vec: Channel Impulse Response in time domain (FFT input/output)
	// ====================================================================
	// MEMORY OPTIMIZATION: Currently DISABLED to reduce memory usage.
	// This array is not allocated or used when Fourier_Transform_of_Channel_Optimized()
	// is employed, which performs continuous-time DFT directly from delay/gain.
	//
	// Memory saved per channel: 3 sectors × NUM_TX_Port × NUM_RX_Port × fft_size × 16 bytes
	// Example: 3 × 8 × 4 × 4096 × 16 = ~6.3 MB per BS-MS channel pair
	//
	// To re-enable for FFT-based approach:
	// 1. Uncomment allocation/deallocation in channel.cpp (Allocate_CHIR_memory, Delete_CHIR_memory)
	// 2. Uncomment usage in channel coefficient calculation (lines ~5093, 5163, 5199)
	// 3. Uncomment clearing in Clear_channel_per_drop() (line ~190)
	// ====================================================================
	ComplexReal **** CHIR_vec                      = NULL;  // NOT ALLOCATED - See note above

	// ====================================================================
	// Ray-level precision mode arrays
	// ====================================================================
	// These arrays store per-ray initial complex gains and angle information
	// to enable precise time evolution with independent Doppler per ray.
	// Only allocated when USE_RAY_LEVEL_DOPPLER is enabled in config.
	// ====================================================================
	ComplexReal ***** ray_gain                     = NULL;  // [sector][tx][rx][cluster][ray]
	Real ***          ray_AOA                      = NULL;  // [cluster][ray][2] - [0]=azimuth, [1]=zenith
	Real ***          ray_AOD                      = NULL;  // [cluster][ray][2] - [0]=azimuth, [1]=zenith
	bool ray_data_allocated                             = false;
	bool ray_angles_precomputed                         = false;

	bool CHIR_allocated                                 = false;
	bool sector_allocated[3]                            = {false, false, false};

	// ====================================================================
	// Element-Level Channel Matrix (ns-3 style)
	// H_usn[cluster](rxElement, txElement) — no beamforming weights
	// ====================================================================

	// raysPreComp: pre-computed polarization matrix per (polTx, polRx, cluster, ray)
	// Value: F_rx_θ·e^(jΦ_θθ)·F_tx_θ + F_rx_θ·√(1/κ)·e^(jΦ_θφ)·F_tx_φ + ...
	ComplexReal **** elem_raysPreComp              = NULL;  // [polTx][polRx][cluster][ray]

	// Spherical unit vector cache (per ray)
	Real ** elem_sinZoA_cosAoA                     = NULL;  // [cluster][ray]
	Real ** elem_sinZoA_sinAoA                     = NULL;
	Real ** elem_cosZoA                            = NULL;
	Real ** elem_sinZoD_cosAoD                     = NULL;
	Real ** elem_sinZoD_sinAoD                     = NULL;
	Real ** elem_cosZoD                            = NULL;

	// Element-level channel matrix: cluster-indexed vector of (totalRx × totalTx) matrices
	std::vector<MatrixXcReal> H_usn;

	// Element-level Doppler storage (Doppler-free initial values)
	std::vector<MatrixXcReal> H_usn_init;  // [cluster] Doppler-free NLOS
	MatrixXcReal H_usn_LOS;                // LOS component (separate Doppler)
	bool H_usn_init_allocated                           = false;

	bool H_usn_allocated                                = false;
	int  H_usn_num_clusters                             = 0;

private:
	
};

#endif