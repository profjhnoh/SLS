#include "const.h"
#include "channel.h"
#include <Eigen/Dense>
#include <vector>

// 220809 jhnoh ---------------------------------------------------
typedef std::pair<Real,int> Real_int_pair;
bool comparator ( const Real_int_pair& l, const Real_int_pair& r);

typedef struct {
    int a;
    int z;
	int p;
	int sector_a;
	int sector_z;
} beam_selection;

// Signal/interference computation helper structs
struct ScenarioConfig {
    int  num_candidate_bs;                // number of candidate BSs to evaluate
    int  sectors_per_bs;                  // 1 for InH 1TRxP, 3 for InH 3TRxP / DU / Rural
    bool use_wraparound;                  // DU/Rural: wraparound topology
    bool apply_penetration_loss;          // DU/Rural: incar + O2I loss in pathloss
    bool use_single_cell_mode;            // single_cell_mode constraint for serving cell selection
    bool use_rsrp_table_for_interference; // DU/Rural: signal_RSRP_gain table lookup for interference
    bool store_los_indoor_info;           // DU/Rural: store link_los, link_indoor
    int  configuration_type;              // 0=ConfigA, 1=ConfigB (InH interference method)
    int  ms_to_bs_wrap_idx;               // wraparound mapping index
};

struct BeamSearchResult {
    Real max_rsrp_dB  = -1e30;   // maximum RSRP (dB)
    int  ue_azimuth_idx  = 0;    // selected UE beam azimuth index
    int  ue_zenith_idx   = 0;    // selected UE beam zenith index
    int  ue_panel_idx    = 0;    // selected UE panel index
    int  sec_azimuth_idx = 0;    // selected BS beam azimuth index
    int  sec_zenith_idx  = 0;    // selected BS beam zenith index
};

struct CandidateCell {
    int   bs_idx       = 0;
    int   sector_idx   = 0;
    int   flat_idx     = 0;       // bs_idx * sectors_per_bs + sector_idx
    Real  signal_dBm   = -1e30;   // tx_power + max_rsrp_dB - pathloss_final
    Real  pathloss_final = 0;
    Real  max_rsrp_dB  = -1e30;
    Real  distance     = 0;
    Real  rms_delay_spread = 0;
    Real  aoa_spread   = 0;
    Real  aod_spread   = 0;
    bool  is_los       = false;
    bool  is_indoor    = false;
    BeamSearchResult beam;
};

class LINK
{
   public:
	    void   Configuration(int);
	    void   Delete_link_memory();
	    // void   Get_signal_interference_old();  // deprecated: original monolithic version
	    void   Get_adj_SECTORS();
	    void   Get_interference();
	    Real Get_angle(LOCATION, LOCATION, LOCATION);
	    Real Get_antgain(CHANNEL *, int , int , int , int );
	    Real Get_MS_antgain(int, int, int, CHANNEL * , int, int, int, int );
	    void   Get_TX_SmallScale_antgain(CHANNEL * , int , int , int , int );
	    void   Get_RX_SmallScale_antgain(int, int, int, CHANNEL * , int, int, int, int);
	    Real Get_RSRP(CHANNEL * , int, int, int, int);
	    void   UE_Initial_Setting(void);
		void   Reset2Default(void);

	    // Refactored signal/interference computation
	    void            Get_signal_interference();
	    ScenarioConfig  build_scenario_config() const;
	    Real            compute_pathloss_final(const ScenarioConfig& cfg, int bs_idx);
	    Real            compute_tx_antenna_gains(CHANNEL* ch, int bs_idx, int sector_idx, int tilt_a, int tilt_z);
	    void            compute_tx_smallscale_gains(CHANNEL* ch, int bs_idx, int sector_idx, int tilt_a, int tilt_z);
	    Real            compute_rsrp(CHANNEL* ch, int sec_number, int sec_z_idx, int sec_a_idx, int mode);
	    BeamSearchResult find_best_tx_beam(CHANNEL* ch, int bs_idx, int sector_idx);
	    CandidateCell   evaluate_candidate_cell(const ScenarioConfig& cfg, int bs_idx, int sector_idx);
	    std::vector<CandidateCell> compute_all_candidate_cells(const ScenarioConfig& cfg);
	    void            select_serving_cell(const ScenarioConfig& cfg, const std::vector<CandidateCell>& candidates);
	    void            compute_interference(const ScenarioConfig& cfg, const std::vector<CandidateCell>& candidates);

	    // Real combined_antgain;
	    ComplexReal TransmitterAntennaGainXLOS_theta       = ComplexReal(0, 0);
	    ComplexReal TransmitterAntennaGainXLOS_pi          = ComplexReal(0, 0);
	    ComplexReal ReceiverAntennaGainXLOS_theta          = ComplexReal(0, 0);
	    ComplexReal ReceiverAntennaGainXLOS_pi             = ComplexReal(0, 0);
	    ComplexReal ReceiverAntennaGainXLOS_theta_panel_2  = ComplexReal(0, 0);
	    ComplexReal ReceiverAntennaGainXLOS_pi_panel_2     = ComplexReal(0, 0);
	    ComplexReal MAX_ReceiverAntennaGainXLOS_theta      = ComplexReal(0, 0);
	    ComplexReal MAX_ReceiverAntennaGainXLOS_pi         = ComplexReal(0, 0);

	    //Real RSRP_antgain;
	    int                 self_bs_idx           = 0;
	    int                 self_ms_idx           = 0;
	    int                 self_sector_idx       = 0;		
	    Real              link_distance         = 0;
	    Real              link_pathloss         = 0;
	    Real              link_coupling_loss    = 0;
	    Real              link_prev_MAX_RSRP    = 0;
	    Real              link_prev_Random_RSRP = 0;
	    Real              link_prev_UE_MAX_RSRP = 0;
	    Real              link_antgain          = 0;
	    Real              link_array_factor[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	    Real              link_RMS_delay_spread = 0;
	    Real              link_AOA_spread       = 0;
	    Real              link_AOD_spread       = 0;
	    Real              distance              = 0;
	    Real              signal                = 0;
	    Real              only_pathloss_signal  = 0;
	    Real              interference          = 0;
	    Real              SINR                  = 0;
	    Real              str_signal            = 0;
	    LOCATION            serving_bs            = {0,0};
 
	    //UE_Initial_Setting
	    Real              incar_loss            = 0;
	    Real              Otoi_loss             = 0;

		bool                high_loss_flag        = 0;


	    unsigned long int   _throughput           = 0;
	    Real              CQI_offset            = 0;
	    
        int                 num_added             = 0;
        int                 num_re_tx             = 0;
        int                 total_num_tx          = 0;
        int                 total_num_re_tx       = 0;
        int                 re_tx_failed          = 0;
        int                 failed                = 0;
        int                 total_num_scheduled   = 0;

        // Per-layer HARQ state (used when g_per_layer_mcs==1). Index = layer 0..num_layers_actual-1.
        int                 per_layer_num_re_tx[4] = {0,0,0,0};
        bool                per_layer_done[4]      = {false,false,false,false};
        Real                per_layer_accum_esinr[4] = {0,0,0,0};  // IR/Chase: accumulated effective SINR (linear) across (re)tx
        // Consecutive scheduled slots in which this pending layer received no new
        // energy (retx allocated fewer streams than the in-flight TB). Reset on
        // service; at 8 the layer is declared starved and force-dropped so a rank
        // that never recovers cannot hold the HARQ process open forever.
        int                 per_layer_starve[4]    = {0,0,0,0};

	    bool                ACK                   = true;
        Real              geometry              = 0;
        int                 selected_a            = 0;
        int                 selected_z            = 0;
        int                 selected_p            = 0;

        int                 azimuth_angle_idx_selected                  = 0;
        int                 zenith_angle_idx_selected                   = 0;
        int                 panel_idx_selected                          = 0;
        int                 azimuth_angle_idx_selected_for_interference = 0;  // selected angle
        int                 zenith_angle_idx_selected_for_interference  = 0;
        int                 panel_idx_selected_for_interference         = 0;
	   
        int                 sector_a              = 0;
        int                 sector_z              = 0;

        int sector_azimuth_angle_idx              = 0;
        int sector_zenith_angle_idx               = 0;

        int    third_adj_sector                   = 0;
        int    second_adj_sector                  = 0;
        int    _sector_in_control                 = 0;

	    MatrixXcReal SmallScale_TX_AntennaGainXLOS_theta;
	    MatrixXcReal SmallScale_TX_AntennaGainXLOS_pi;
	    MatrixXcReal SmallScale_RX_AntennaGainXLOS_theta;
	    MatrixXcReal SmallScale_RX_AntennaGainXLOS_pi;
	    MatrixXcReal SmallScale_RX_AntennaGainXLOS_theta_panel_2;
	    MatrixXcReal SmallScale_RX_AntennaGainXLOS_pi_panel_2;
	    MatrixXcReal MAX_SmallScale_RX_AntennaGainXLOS_theta;
	    MatrixXcReal MAX_SmallScale_RX_AntennaGainXLOS_pi;

	    Real_int_pair * static_gain             = NULL;
	    beam_selection  * analog_beam_selection   = NULL;
	    Real *          intf_w_rnd_RSRP         = NULL;
        Real *          LS_gain                 = NULL;
        int *             adj_sector              = NULL;
        int *             rand_sec_a              = NULL;
        int *             rand_sec_z              = NULL;

		// For comp implementation by jhnoh 230109
		int    comp_sector_idx                 = -1;

        Real comp_signal_strength[2]            = {-INFINITY,-INFINITY};
        Real comp_interf_strength[2]            = {0,0};		
		Real comp_max_RSRP                   = 0;
		Real comp_inf_RSRP                   = 0;
		int    comp_azimuth_angle_idx_selected = 0;
		int    comp_zenith_angle_idx_selected  = 0;
		int    comp_panel_idx_selected         = 0;
		int    comp_sector_a                   = 0;
		int    comp_sector_z                   = 0;

		Real link_RSRP = 0;

		bool link_los = 0;
		bool link_indoor = 0;
};