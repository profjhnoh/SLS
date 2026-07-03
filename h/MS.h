#ifndef _MS_H_
#define _MS_H_
#include "const.h"
#include <complex>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/FFT>

using namespace Eigen;


LL_CURVE_LOOK eesm_methology(INPUT_TO_LLS_MAPPING input, int _self_idx);
LL_CURVE_LOOK rbir_methology(INPUT_TO_LLS_MAPPING input, int _self_idx);
LL_CURVE_LOOK mib_methology(INPUT_TO_LLS_MAPPING input, int _self_idx);

class MS
{
    public:
		void   Configuration(int,int);
		void   Set_MS_Location();
		void   Channel_Update_MIMO(int);
		void   Quantization_of_Ch(void);
		void   Quantization_of_Ch_CSIRS_2(void);
		void   Quantization_of_Ch_Type2(void);    // TS 38.214 §5.2.2.2.3 Type II (Rel-15, rank 1)
		void   Quantization_of_Ch_EType2(void);   // TS 38.214 §5.2.2.2.5 eType II (Rel-16, rank 1..4)
		void   Compute_RI(void);                  // Rank Indicator selection (SVD capacity)
		void   Store_CSI_for_TDD(void);  // TDD: Store full channel matrix for reciprocity-based precoding
		void   CQI_Update(void);
		void   Received_SINR(void);
		void   Fourier_Transform_of_Channel(int);
		void   Fourier_Transform_of_Channel_Optimized(int);  // Optimized version using recursive phasor updates
		void   Fourier_Transform_H_usn(int);                  // H_usn DFT (no beamforming)
		void   Fourier_Transform_WithBF(int);                // Element-level DFT + BF weights → H_m
		void   Declare_ch_matrix(void);
		void   Receive_DL(int);
		void   Find_Allocated_Rbs_Mcs(int);
		void   Rate_matching(void);
		void   Compute_Num_Bin_Symbol(void);
		void   Compute_Effective_SINR(void);
		void   Compute_lls_mapping(void);
		void   AddThroughput(Real prob);
		Real Return_initial_packet_error(void) const;
		void   Kill_Memory(void);
		void   Receive_DL_mTRP( void );
		MatrixXcReal MMSE_Receiver (int, int);
		Real compute_tone_SINR_MRC(int rb_idx);
		Real compute_tone_SINR(int rb_idx);
		void   compute_CQI_offset(Real Received_SINR, int rb_idx);
		Real Return_Throughput(int) const;
		void   Allocate_memory(void);
		void   Delete_memory(void);
		void   Find_Allocated_Rbs_Mcs_mTRP  ( void );
		Real Compute_Tone_SINR_NCJT       ( int, int );
		Real Compute_Tone_SINR_CJT        ( int, int );
		void   Compute_RBs_SINR             ( void );
		void   Compute_Transport_Block_Size ( void );
		Real BLER_Mapping                 ( void );
		INPUT_TO_LLS_MAPPING Return_mapping_info( void );
		void   Reset2Default ( void );
		void   CQI_Update_CJT( void );
		VectorXcReal PowerIteration_DominantEigenvector(const MatrixXcReal& M, int max_iterations = 5);

		int              unscheduled_stack = 0;
		int              self_idx          = 0;
		int              self_sector_idx   = 0;
		int              nearest_bs_idx    = 0;
		LOCATION         loc               = {0,0};
		LOCATION3D       d_rx[4][8][2][2][2];//  antenna element location vector [M][N][P][Mg][Ng] (increased for safety)

		// array antenna rotate angle
		Real           alpha             = 0;
		Real           beta              = 0;
		Real           gamma             = 0;

		Real           speed             = 0;
		Real           moving_direction  = 0;
		Real           moving_direction_azimuth = 0;

		bool           Indoor            = false;

		//receive_downlink
		int               receive_self_ms_idx = 0;
		Real            received_sinr_avg   = 0;
		Real            ESINR_linear        = 0;
		long unsigned int sum_cqi             = 0;
		long unsigned int sum_mcs_type        = 0;
		long double       total_esinr         = 0;
		long double       total_estimate_sinr = 0;
		int               _info_bits          = 0;
		int               _mcs_idx            = 0;
		int               _cqi_idx            = 0;
		// Per-layer MCS/CQI/TBS (used when g_per_layer_mcs==1). Index = layer (0..num_layers_actual-1),
		// matching rbs_rx_per_layer ordering. Legacy path uses _mcs_idx/_cqi_idx (= layer 0).
		int               _mcs_per_layer[4]       = {0,0,0,0};
		int               _cqi_per_layer[4]       = {0,0,0,0};
		int               _info_bits_per_layer[4] = {0,0,0,0};
		Real             _avr_sinr           = 0;
		int               _mod_type           = 0;
		int               num_rx_rb           = 0;
		int               _num_traffic        = 0;
		Real            _code_rate          = 0;
		Real            _snr_scale_factor   = 0;
 		int               num_scheduled       = 0;
		int               num_comp_tx         = 0;
			
		// cqi_offset
		Real            CQI_offset          = 0;
		int               num_added           = 0;
		Real            ue_BLER_Value       = 0;

		//int *propagation_condition;
		int               floor_idx           = 0;
		Real            ms_d_in             = 0;
		Real            MS_HEIGHT_FINAL     = 0;
		LOCATION          pos_indoor          = {0,0};

		bool              comp_mode_rx_flag   = 0;
		int               num_rx_rbs          = 0;
		int               re_tx_failed_cnt    = 0;
		int               cqi_offset          = 0;

		// OLLA (Outer Loop Link Adaptation) offset in dB
		Real              olla_offset         = 0.0;
		// Moving window BLER tracking for OLLA
		std::vector<bool> olla_history;               // Circular buffer: true=ACK, false=NACK
		int               olla_history_idx    = 0;    // Current index in circular buffer
		int               olla_nack_count     = 0;    // Number of NACKs in current window

		int               pmi_l               = 0;
		int               pmi_m               = 0;
		int               pmi_n               = 0;

		// Rank Indicator (RI): UE-selected rank for eType II when g_rank_adaptive=1.
		// Default 1. Updated in Compute_RI() each time slot.
		int               self_RI             = 1;
		Real              su_capacity         = 0.0;   // Best SU Shannon capacity for scheduler SU-vs-MU

		// Actual number of layers allocated to this UE on the most recent scheduled RB.
		// Set inside compute_tone_SINR() as my_streams.size(). Used by Rate_matching to
		// compute the correct N_info (transport-block size scales linearly with layers).
		int               num_layers_actual   = 1;

		// Per-layer reception buffers (used when num_layers_actual > 1).
		// rbs_rx_per_layer[l] = vector of post-equalizer SINR (linear) for layer l across scheduled RBs.
		// per_layer_ESINR[l]  = EESM-aggregated per-layer effective SINR.
		// per_layer_BLER[l]   = looked-up BLER per layer.
		// TB ACK = AND of all layer ACKs (Bernoulli trials).
		std::vector< std::vector<Real> > rbs_rx_per_layer;   // [layer][rb_in_alloc]
		Real              per_layer_ESINR[4]    = {0,0,0,0};
		Real              per_layer_BLER[4]     = {0,0,0,0};

		// PMI quality / channel-rank statistics, accumulated in Quantization_of_Ch_EType2().
		// Used to diagnose RI selection (σ-spread) and PMI codebook quantization loss.
		Real              pmi_sigma2_sum[5]     = {0,0,0,0,0};  // σ_l² (index 1..4)
		Real              pmi_overlap_sum[5]    = {0,0,0,0,0};  // |w_l^H v_l|² per column (legacy)
		Real              pmi_cos2_sum[5]       = {0,0,0,0,0};  // cos²(θ_l): principal angles
		Real              pmi_chordal_sum       = 0.0;          // R − Σ cos²(θ_l) per sample
		int               pmi_quality_count     = 0;            // # (subband × time) samples

		// Per-layer SINR comparison (predicted vs received)
		Real              pred_perlayer_se_sum[5]  = {0,0,0,0,0};  // Compute_RI's estimated log2(1+sv²·ρ_norm) (bits/s/Hz)
		Real              recv_perlayer_sinr_sum[5]= {0,0,0,0,0};  // Receiver's actual per-layer SINR (linear)
		int               recv_perlayer_count[5]   = {0,0,0,0,0};

		Real ***            CQI             = NULL;
		Real ***            CQI_comp        = NULL;

		MatrixXcReal **          H_m             = NULL;
		MatrixXcReal *           H_m_elem        = NULL;  // Element-level freq-domain channel [rb_idx](totalRx, totalTx)
		PMI_FEEDBACK ***      PMI             = NULL;
		MatrixXcReal ***      PMI_vector      = NULL;   // N_tx × R (rank 1 default). col 0 = dominant layer.
		MatrixXcReal ***      CSI_matrix      = NULL;  // TDD: Full channel matrix (NUM_RX_Port × NUM_TX_Port) for reciprocity-based precoding
		SCHEDULE_DECISION **  ppSchedulerRead = NULL;
		//LOCATION *            pos             = NULL;
		ArrayXReal  *            TLSPs           = NULL;
		ArrayXReal  *            LSPs            = NULL;

		std::vector <Real> rbs_rx; 
		std::vector <int> rb_indices;

		Real geometry = 0;
		Real coupling_loss = 0;
		Real max_RSRP = 0;
		Real noise_linear_cached = 0;  // Cached noise value in linear scale

		Eigen::FFT<Real> fft;
};

#endif

Real ComputeBLER(LL_CURVE_LOOK);
void Get_Geometry(Real geometry_ue, int ue_idx);
void Get_Wideband_SIR(Real sir_ue, int ue_idx);
void Coupling_Loss_Get(Real cp_ue, int ue_idx);
void Get_RMS_delay_spread(Real RMS_delay_spread_ue, int ue_idx);
void Get_circular_angle_spread_AOA(Real angle_spread_AOA_ue, int ue_idx);
void Get_circular_angle_spread_AOD(Real angle_spread_AOD_ue, int ue_idx);