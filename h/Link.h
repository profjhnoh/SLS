#include "const.h"
#include "channel.h"
#include <Eigen/Dense>

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

class LINK
{
   public:
	    void   Configuration(int);
	    void   Delete_link_memory();
	    void   Get_signal_interference();
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
        //Real Get_antgain(LOCATION, LOCATION, int, int, int, int);
	    //void Get_TX_SmallScale_antgain(LOCATION, int, int, int, int);
	    //Real Get_RSRP(int, int, int, int, int, int);

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