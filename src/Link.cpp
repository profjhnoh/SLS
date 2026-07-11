#include "common.h"
#include <chrono>  // For performance timing

// Real precision에 맞는 Eigen 매트릭스 타입 정의
//#ifdef USE_FLOAT
//    typedef Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> MatrixXcReal;
//#else
//    typedef Eigen::MatrixXcd MatrixXcReal;
//#endif

// ============================================================================
// Complex Exponential Optimization
// ============================================================================
// Set to 1 to use optimized cos/sin instead of exp() for complex exponentials
// This can provide 2-4x speedup for exp(i*theta) calculations
#define USE_FAST_COMPLEX_EXP 1

// Set to 1 to enable detailed profiling of exp() optimization
#define PROFILE_FAST_EXP 0

#if USE_FAST_COMPLEX_EXP
// Optimized version: exp(i*theta) = cos(theta) + i*sin(theta)
// For pure imaginary exponent: exp(complex(0, theta))
inline ComplexReal fast_exp_imag(Real theta) {
    return ComplexReal(cos(theta), sin(theta));
}

// For complex exponent where real part is 0: exp(complex(0, imag))
inline ComplexReal fast_exp(const ComplexReal& z) {
    if (z.real() == 0.0) {
        // Pure imaginary: exp(i*imag) = cos(imag) + i*sin(imag)
        return ComplexReal(cos(z.imag()), sin(z.imag()));
    } else {
        // Has real part: exp(real + i*imag) = exp(real) * [cos(imag) + i*sin(imag)]
        Real exp_real = std::exp(z.real());
        return ComplexReal(exp_real * cos(z.imag()), exp_real * sin(z.imag()));
    }
}

#if PROFILE_FAST_EXP
// Profiling wrapper
static long long fast_exp_call_count = 0;
static double fast_exp_total_time = 0.0;
inline ComplexReal fast_exp_profiled(const ComplexReal& z) {
    auto start = std::chrono::high_resolution_clock::now();
    ComplexReal result = fast_exp(z);
    auto end = std::chrono::high_resolution_clock::now();
    fast_exp_total_time += std::chrono::duration<double, std::micro>(end - start).count();
    fast_exp_call_count++;
    return result;
}
#define FAST_EXP(z) fast_exp_profiled(z)
#else
#define FAST_EXP(z) fast_exp(z)
#endif

#else
// Original version: use standard exp()
#if PROFILE_FAST_EXP
// Profiling wrapper for original exp()
static long long exp_call_count = 0;
static double exp_total_time = 0.0;
inline ComplexReal exp_profiled(const ComplexReal& z) {
    auto start = std::chrono::high_resolution_clock::now();
    ComplexReal result = exp(z);
    auto end = std::chrono::high_resolution_clock::now();
    exp_total_time += std::chrono::duration<double, std::micro>(end - start).count();
    exp_call_count++;
    return result;
}
#define FAST_EXP(z) exp_profiled(z)
#else
#define FAST_EXP(z) exp(z)
#endif
#endif
// ============================================================================

//Real Get_distance(LOCATION, LOCATION);

Real Get_LCS_theta(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);
Real Get_LCS_pi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi);

Real Get_BS_antenna_pattern(Real theta_GCS, Real pi_GCS, int bs_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2);
Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2, int port_idx = -1);

ComplexReal  Get_BS_antenna_field_pattern(LOCATION interferer, int _bs_idx, int _sec_idx, int tilt_z_idx, int tilt_a_idx, Real v_angle_theta, Real h_angle_pi, Real F_theta_GCS_P1, Real F_pi_GCS_P1, Real F_theta_GCS_P2, Real F_pi_GCS_P2, ComplexReal  &F_tx_theta, ComplexReal  &F_tx_pi, const ClusterVR* sns_vr_ptr = nullptr, Real sns_rolloff_C = 0.0);
ComplexReal  Get_UE_antenna_field_pattern(LOCATION UE , int M, int N, int P, int ms_idx, int sector_idx, int tilt_z_idx, int tilt_a_idx, Real v_angle_theta, Real h_angle_pi, Real F_theta_GCS_P1, Real F_pi_GCS_P1, Real F_theta_GCS_P2, Real F_pi_GCS_P2, ComplexReal  &F_rx_theta_panel_1, ComplexReal  &F_rx_pi_panel_1, ComplexReal  &F_rx_theta_panel_2, ComplexReal  &F_rx_pi_panel_2);

MatrixXReal Get_distance_angular(Real a, Real b, Real c, Real x, Real y, Real z);
void LCS_Antenna_field_to_GCS_antenna_pattern(Real alpha, Real beta, Real gamma, Real GCS_angle_theta, Real GCS_angle_pi, Real LCS_field_theta, Real LCS_field_pi, Real &GCS_field_theta, Real &GCS_field_pi); //GCS_angle degree
void PIot_sector_antenna_gain();
void PIot_ue_antenna_gain();
ANTENNA_FIELD Antenna_field_in_LCS(Real LCS_antenna_gain_in_dB, Real zeta);
Real dot(LOCATION3D a, LOCATION3D b);
ANTENNA_FIELD LCS_Antenna_field_to_GCS_antenna_pattern_rev(Real alpha, Real beta, Real gamma, SPHERICAL_ANGLE GCS_angle, ANTENNA_FIELD LCS_field);
LOCATION3D Transform_angle_to_spheical_vector(Real phi, Real theta);
SPHERICAL_ANGLE GCS_Angle_to_LCS_Angle(Real alpha, Real beta, Real gamma, SPHERICAL_ANGLE GCS_angle);
Real Sector_Antenna_gain_in_dB(SPHERICAL_ANGLE LCS_angle);
Real UE_Antenna_gain_in_dB(SPHERICAL_ANGLE LCS_angle); //LCS_angle : degree

ComplexReal  Generate_vertical_virtualization_weight_link(Real theta, Real dv, Real k, Real K);
ComplexReal  Generate_horizontal_virtualization_weight_link(Real theta, Real phi, Real dH, Real l, Real L);

Real Transform_angle_minus_180_to_plus_180(Real x);
Real Transform_angle_0_to_plus_180(Real x);


// 220809 jhnoh ---------------------------------------------------
bool comparator ( const Real_int_pair& l, const Real_int_pair& r) 
{
    return l.first > r.first;
}
//-----------------------------------------------------------------

void LINK::Reset2Default( void )
{
	    //Real RSRP_antgain;
	    self_bs_idx           = 0;
	    self_ms_idx           = 0;
	    self_sector_idx       = 0;		
	    link_distance         = 0;
	    link_pathloss         = 0;
	    link_coupling_loss    = 0;
	    link_prev_MAX_RSRP    = 0;
	    link_prev_Random_RSRP = 0;
	    link_prev_UE_MAX_RSRP = 0;
	    link_antgain          = 0;
	    link_array_factor[0]  = 0;
		link_array_factor[1]  = 0;
		link_array_factor[2]  = 0;
		link_array_factor[3]  = 0;
		link_array_factor[4]  = 0;
		link_array_factor[5]  = 0;
		link_array_factor[6]  = 0;
		link_array_factor[7]  = 0;
		link_array_factor[8]  = 0;
		link_array_factor[9]  = 0;
		link_array_factor[10]  = 0;
		link_array_factor[11]  = 0;
	    link_RMS_delay_spread = 0;
	    link_AOA_spread       = 0;
	    link_AOD_spread       = 0;
	    distance              = 0;
	    signal                = 0;
	    only_pathloss_signal  = 0;
	    interference          = 0;
	    SINR                  = 0;
	    str_signal            = 0;
	    serving_bs            = {0,0};
 
	    
	    incar_loss            = 0;
	    Otoi_loss             = 0;
	    _throughput           = 0;
	    CQI_offset            = 0;
	    num_added             = 0;
        num_re_tx             = 0;
        total_num_tx          = 0;
        total_num_re_tx       = 0;
        re_tx_failed          = 0;
        failed                = 0;
        total_num_scheduled   = 0;
        total_first_tx_layers  = 0;
        failed_first_tx_layers = 0;
        for (int l = 0; l < 4; l++)
        {
            per_layer_num_re_tx[l]   = 0;
            per_layer_done[l]        = false;
            per_layer_accum_esinr[l] = 0.0;
            per_layer_starve[l]      = 0;
        }

	    ACK                   = true;
        geometry              = 0;
        selected_a            = 0;
        selected_z            = 0;
        selected_p            = 0;

        azimuth_angle_idx_selected                  = 0;
        zenith_angle_idx_selected                   = 0;
        panel_idx_selected                          = 0;
        azimuth_angle_idx_selected_for_interference = 0;  // selected angle
        zenith_angle_idx_selected_for_interference  = 0;
        panel_idx_selected_for_interference         = 0;
	   
        sector_a                           = 0;
        sector_z                           = 0;

        sector_azimuth_angle_idx           = 0;
        sector_zenith_angle_idx            = 0;

        third_adj_sector                   = 0;
        second_adj_sector                  = 0;
        _sector_in_control                 = 0;

	    SmallScale_TX_AntennaGainXLOS_theta        .resize(0,0);
	    SmallScale_TX_AntennaGainXLOS_pi           .resize(0,0);
	    SmallScale_RX_AntennaGainXLOS_theta        .resize(0,0);
	    SmallScale_RX_AntennaGainXLOS_pi           .resize(0,0);
	    SmallScale_RX_AntennaGainXLOS_theta_panel_2.resize(0,0);
	    SmallScale_RX_AntennaGainXLOS_pi_panel_2   .resize(0,0);
	    MAX_SmallScale_RX_AntennaGainXLOS_theta    .resize(0,0);
	    MAX_SmallScale_RX_AntennaGainXLOS_pi       .resize(0,0);

	    //Real_int_pair * static_gain             = NULL;
	    //beam_selection  * analog_beam_selection   = NULL;
	    //Real *          intf_w_rnd_RSRP         = NULL;
        //Real *          LS_gain                 = NULL;
        //int *             adj_sector              = NULL;
        //int *             rand_sec_a              = NULL;
        //int *             rand_sec_z              = NULL;	
		for (int sec_idx = 0; sec_idx < num_SECTORS; sec_idx++ )
		{
			LS_gain              [sec_idx] = 0;
			adj_sector           [sec_idx] = 0;
			static_gain          [sec_idx].first  = 0;
			static_gain          [sec_idx].second = 0;
			analog_beam_selection[sec_idx].a = 0;
			analog_beam_selection[sec_idx].z = 0;
			analog_beam_selection[sec_idx].p = 0;
			analog_beam_selection[sec_idx].sector_a = 0;
			analog_beam_selection[sec_idx].sector_z = 0;
			intf_w_rnd_RSRP	     [sec_idx] = 0;
			rand_sec_a           [sec_idx] = 0;
			rand_sec_z           [sec_idx] = 0;				
		}
}


void LINK::Configuration(int _self_ms_idx)
{
	self_ms_idx = _self_ms_idx;
	UE_Initial_Setting();
	Get_signal_interference();
	Get_adj_SECTORS();

	SINR = 10. * log10(dBm2linear(signal) / (dBm2linear(interference) + noise));
}



#if 0  // deprecated: original monolithic Get_signal_interference — replaced by refactored version below
void LINK::Get_signal_interference_old()
{

	Real strongest_signal;
	Real strongest_only_pathloss_signal;
	Real pathloss;
	Real antgain;
	Real random_antgain;
	Real max_array_factor_gain = 0;
	Real random_array_factor_gain = 0;  /// for interference
	Real UE_max_array_factor_gain = 0;
	Real UE_random_array_factor_gain = 0;
	Real RMS_delay;
	Real AOA_spread;
	Real AOD_spread;
	Real distance;
	Real UE_antenna_gain = 0.;
	interference = 0.;
	Real RSRP_antgain = 0;
	Real RSRP_signal = 0;
	Real MAX_RSRP_antgain = 0.;
	//Real MAX_RSRP_signal = 0.;
	Real random_RSRP_antgain = 0;
	//link_array_factor[] = 0;

	//jhnoh 220810
	CHANNEL * channel_of_interest;

	if (TYPE == 11 && (num_Indoor_TRxP == 1)) //// 5G InH, 12 TRxP
	{
		///////////////////////////////////////////////////////////////// Get Signal ////////////////
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP ; bs_idx++)
		{
			channel_of_interest = &channel[bs_idx][self_ms_idx];
			pathloss   = channel_of_interest->pathloss;
			RMS_delay  = channel_of_interest->RMS_delay_spread;
			AOA_spread = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread = channel_of_interest->circular_angle_spread_AOD;
			distance   = channel_of_interest->distance;

			pathloss = pathloss + ms[self_ms_idx].LSPs[bs_idx](0);   ////// Add Shadow Fading (use stored LSP SF)
			channel_of_interest->pathloss_final = pathloss;

			///*
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////   Antenna gain - NEW version //////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
			{
				for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
				{
					antgain = Get_antgain(channel_of_interest, bs_idx, 0, a, z);   ///// MAX antgain - combined gain

					Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, 0, a, z);

					RSRP_antgain = Get_RSRP(channel_of_interest,0,z,a, 0);

					if (a == 0 && z == 0)  // first loop
					{
						MAX_RSRP_antgain = RSRP_antgain;

						azimuth_angle_idx_selected = selected_a;
						zenith_angle_idx_selected  = selected_z;
						panel_idx_selected         = selected_p;

						sector_a = a;
						sector_z = z;


					}
					else if (MAX_RSRP_antgain < RSRP_antgain)   // find max 
					{
						MAX_RSRP_antgain = RSRP_antgain;

						azimuth_angle_idx_selected = selected_a;
						zenith_angle_idx_selected  = selected_z;
						panel_idx_selected         = selected_p;

						sector_a = a;
						sector_z = z;
					}
					else
					{
						MAX_RSRP_antgain = MAX_RSRP_antgain;
					}
				}
			}

			signal = bs_maxpower + MAX_RSRP_antgain - pathloss;

			// 220826 jhnoh
			static_gain          [bs_idx].first  = signal;
			static_gain          [bs_idx].second = bs_idx;

			analog_beam_selection[bs_idx].a = azimuth_angle_idx_selected;
			analog_beam_selection[bs_idx].z = zenith_angle_idx_selected;
			analog_beam_selection[bs_idx].p = panel_idx_selected;

			analog_beam_selection[bs_idx].sector_a = sector_a;
			analog_beam_selection[bs_idx].sector_a = sector_z;


			//////////////////////////////////////////////////////////////////////////////////////////////////////////

			if (bs_idx == 0)
			{
				strongest_signal = signal;
				link_prev_MAX_RSRP = MAX_RSRP_antgain;
				link_prev_Random_RSRP = random_RSRP_antgain;

				self_bs_idx = bs_idx;
				self_sector_idx = 0;
				link_pathloss = pathloss;
				link_antgain = RSRP_antgain;
				link_RMS_delay_spread = RMS_delay;
				link_AOA_spread = AOA_spread;
				link_AOD_spread = AOD_spread;
				link_distance = distance;

				azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
				zenith_angle_idx_selected_for_interference = zenith_angle_idx_selected;
				panel_idx_selected_for_interference = panel_idx_selected;

				sector_azimuth_angle_idx = sector_a;
				sector_zenith_angle_idx = sector_z;

				// indoor���� num_BS = num_sector , bs = sector 

				_sector_in_control = self_bs_idx;

			}
			else if (signal > strongest_signal)
			{

				if (Configuration_Type == 0)  // config A
				{
					interference += dBm2linear(strongest_signal); //( pow(10., (strongest_signal) / 10.);
				}
				else
				{
					//interference += pow(10., (strongest_signal - link_prev_MAX_RSRP + link_prev_Random_RSRP) / 10.);
				}

				strongest_signal = signal;
				//link_prev_MAX_RSRP = MAX_RSRP_antgain;
				//link_prev_Random_RSRP = random_RSRP_antgain;

				self_bs_idx = bs_idx;
				self_sector_idx = 0;
				link_pathloss = pathloss;
				link_antgain = RSRP_antgain;
				link_RMS_delay_spread = RMS_delay;
				link_AOA_spread = AOA_spread;
				link_AOD_spread = AOD_spread;
				link_distance = distance;

				azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
				zenith_angle_idx_selected_for_interference = zenith_angle_idx_selected;
				panel_idx_selected_for_interference = panel_idx_selected;

				sector_azimuth_angle_idx = sector_a;
				sector_zenith_angle_idx = sector_z;


				_sector_in_control = self_bs_idx;

			}
			else
			{
				if (Configuration_Type == 0)  // config A
				{
					interference += dBm2linear(signal); //pow(10., (signal) / 10.);
				}
				else
				{
					//interference += pow(10., (signal - MAX_RSRP_antgain + random_RSRP_antgain) / 10.);
				}
			}
		}

		/////////////////////// Get interference ///////////////////////////////////////////////////////////////////////////////////////////////////////
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			channel_of_interest = &channel[bs_idx][self_ms_idx];
			pathloss = channel_of_interest->pathloss_final;

			if (Configuration_Type == 0)
			{
			}
			else
			{
				///////////////////////////////////////////////////////////////// random beam interference
				int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
				int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

				random_antgain = Get_antgain(channel_of_interest, bs_idx, 0, aa, zz);   ///// MAX antgain - combined gain
				Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, 0, aa, zz);
				random_RSRP_antgain = Get_RSRP(channel_of_interest, 0, zz, aa, 1);
				/////////////////////////////////////////////////////////////////

				signal = bs_maxpower + random_RSRP_antgain - pathloss;
			}

			if (bs_idx != self_bs_idx)
			{
				if (Configuration_Type == 0)  // config A
				{

				}
				else
				{
					interference += dBm2linear(signal); //pow(10., (signal) / 10.);
				}
			}
			intf_w_rnd_RSRP[bs_idx] = signal;
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}
	else if (TYPE == 11 && (num_Indoor_TRxP == 3)) //// 5G InH, 36 TRxP
	{
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP ; bs_idx++)
		{
			// 220811 jhnoh
			channel_of_interest                 = &channel[bs_idx][self_ms_idx];
			pathloss                            = channel_of_interest->pathloss;
			RMS_delay                           = channel_of_interest->RMS_delay_spread;
			AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
			distance                            = channel_of_interest->distance;
			pathloss                            = pathloss + ms[self_ms_idx].LSPs[bs_idx](0);   ////// Add Shadow Fading (use stored LSP SF)
			channel_of_interest->pathloss_final = pathloss;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				///*
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////   Antenna gain - NEW version //////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
				{
					for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
					{
						//antgain = Get_antgain(channel_of_interest->self_bs, channel_of_interest->self_ms, bs_idx, sector_idx-1, a, z);   ///// MAX antgain - combined gain
						antgain = Get_antgain(channel_of_interest, bs_idx, sector_idx, a, z);   ///// MAX antgain - combined gain

						//Get_TX_SmallScale_antgain(channel_of_interest->self_bs, bs_idx, sector_idx - 1, a, z);
						Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, sector_idx, a, z);

						RSRP_antgain = Get_RSRP(channel_of_interest,sector_idx,z,a, 0);

						if (a == 0 && z == 0)  // first loop
						{
							MAX_RSRP_antgain = RSRP_antgain;

							azimuth_angle_idx_selected = selected_a;
							zenith_angle_idx_selected  = selected_z;
							panel_idx_selected         = selected_p;

							sector_a = a;
							sector_z = z;
						}
						else if (MAX_RSRP_antgain < RSRP_antgain)   // find max 
						{
							MAX_RSRP_antgain = RSRP_antgain;

							azimuth_angle_idx_selected = selected_a;
							zenith_angle_idx_selected  = selected_z;
							panel_idx_selected         = selected_p;

							sector_a = a;
							sector_z = z;
						}
						else
						{
							MAX_RSRP_antgain = MAX_RSRP_antgain;
						}
					}
				}

				signal = bs_maxpower + MAX_RSRP_antgain - pathloss;

				// 220808 jhnoh
				static_gain  [bs_idx * 3 + (sector_idx)].first     = signal;
				static_gain  [bs_idx * 3 + (sector_idx)].second    = bs_idx * 3 + (sector_idx);

				analog_beam_selection[bs_idx * 3 + (sector_idx)].a = azimuth_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].z = zenith_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].p = panel_idx_selected;

				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_a = sector_a;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_z = sector_z;


				//sort(A.begin(),A.end(), comparator);

				//*/
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////

				// 220826 jhnoh 확인필요!!
				//if (bs_idx == 0)
				if (bs_idx == 0 && sector_idx == 0 )
				{
					strongest_signal      = signal;
					//link_prev_MAX_RSRP    = MAX_RSRP_antgain;
					//link_prev_Random_RSRP = random_RSRP_antgain;
					link_RSRP             = MAX_RSRP_antgain;

					self_bs_idx           = bs_idx;
					self_sector_idx       = sector_idx;
					link_pathloss         = pathloss;
					link_antgain          = RSRP_antgain;
					link_RMS_delay_spread = RMS_delay;
					link_AOA_spread       = AOA_spread;
					link_AOD_spread       = AOD_spread;
					link_distance         = distance;

					azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
					zenith_angle_idx_selected_for_interference  = zenith_angle_idx_selected;
					panel_idx_selected_for_interference         = panel_idx_selected;

					sector_azimuth_angle_idx = sector_a;
					sector_zenith_angle_idx  = sector_z;

					// 
					//adj_sector2[0] = _sector_in_control;
					_sector_in_control =  self_bs_idx * 3 + (self_sector_idx);

				}
				else if (signal > strongest_signal)
				{
					if (Configuration_Type == 0)  // config A
					{
						interference += dBm2linear(strongest_signal);//pow(10., (strongest_signal) / 10.);
					}
					else
					{
						//interference += pow(10., (strongest_signal - link_prev_MAX_RSRP + link_prev_Random_RSRP) / 10.);
					}

					strongest_signal      = signal;
					//link_prev_MAX_RSRP    = MAX_RSRP_antgain;
					//link_prev_Random_RSRP = random_RSRP_antgain;

					self_bs_idx           = bs_idx;
					self_sector_idx       = sector_idx;
					link_pathloss         = pathloss;
					link_RSRP             = MAX_RSRP_antgain;
					
					link_antgain          = RSRP_antgain;
					link_RMS_delay_spread = RMS_delay;
					link_AOA_spread       = AOA_spread;
					link_AOD_spread       = AOD_spread;
					link_distance         = distance;

					azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
					zenith_angle_idx_selected_for_interference  = zenith_angle_idx_selected;
					panel_idx_selected_for_interference         = panel_idx_selected;

					sector_azimuth_angle_idx = sector_a;
					sector_zenith_angle_idx  = sector_z;

					_sector_in_control = self_bs_idx * 3 + (self_sector_idx);

				}
				else
				{
					if (Configuration_Type == 0)  // config A
					{
						interference += dBm2linear(signal);  //pow(10., (signal) / 10.);
					}
					else
					{
						//interference += pow(10., (signal - MAX_RSRP_antgain + random_RSRP_antgain) / 10.);
					}
				}
			}
			
		}

		#if 0
		// Select Comp Sector
		int    comp_Rx_panel_idx            = 0;
		Real temp_signal_strength           = 0;

		comp_sector_idx                     = -1;

		comp_signal_strength[1]             = -INFINITY;
		comp_interf_strength[1]             = 0;

        comp_max_RSRP                       = 0;
        comp_inf_RSRP                       = 0;
        comp_azimuth_angle_idx_selected     = 0;
        comp_zenith_angle_idx_selected      = 0;
        comp_panel_idx_selected             = 0;
        comp_sector_a                       = 0;
        comp_sector_z                       = 0;

		if ( panel_idx_selected_for_interference == 0 )
		{
			comp_Rx_panel_idx = 1;
		}
		else
		{
			comp_Rx_panel_idx = 0;
		}

		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			// 220811 jhnoh
			channel_of_interest                 = &channel[bs_idx][self_ms_idx];
			pathloss                            = channel_of_interest->pathloss;
			RMS_delay                           = channel_of_interest->RMS_delay_spread;
			AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
			distance                            = channel_of_interest->distance;
			pathloss                            = channel_of_interest->pathloss_final;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				if ((bs_idx != self_bs_idx) || ( sector_idx != self_sector_idx))
				{			
					// Analog Beamforming 
					Real Tmp_RSRP_gain = 0;
					Real Max_RSRP_gain = 0;
					Real Inf_RSRP_gain = 0;

					int ue_a_idx_selected  = -1;
					int ue_z_idx_selected  = -1;
					int ue_p_idx_selected  = -1;
					int sec_a_idx_selected = -1;
					int sec_z_idx_selected = -1;

					// Sector beam index loop
					for (int sec_a_idx = 0; sec_a_idx < tilt_azimuth_angle_LCS_size; sec_a_idx++)
					{
						for (int sec_z_idx = 0; sec_z_idx < tilt_zenith_angle_LCS_size; sec_z_idx++)
						{
							// UE beam index loop
							for (int ue_z_idx = 0; ue_z_idx < ue_tilt_zenith_angle_LCS_size; ue_z_idx++)
							{
								for (int ue_a_idx = 0; ue_a_idx < ue_tilt_azimuth_angle_LCS_size; ue_a_idx++)
								{						
									Tmp_RSRP_gain = 10*log10(channel_of_interest->signal_RSRP_gain
															[sector_idx][sec_z_idx][sec_a_idx]
															[ue_z_idx][ue_a_idx][comp_Rx_panel_idx]);


									if ( Tmp_RSRP_gain > Max_RSRP_gain  || 
									( sec_a_idx == 0 && sec_z_idx == 0 && ue_z_idx == 0 && ue_a_idx == 0 ))
									{
										Max_RSRP_gain       = Tmp_RSRP_gain;
										ue_a_idx_selected   = ue_a_idx;
										ue_z_idx_selected   = ue_z_idx;
										ue_p_idx_selected   = comp_Rx_panel_idx;
										sec_a_idx_selected  = sec_a_idx;
										sec_z_idx_selected  = sec_z_idx;
									}
								}
							} // UE beam index loop end
						}
					} // Sector beam index loop end

					temp_signal_strength = bs_maxpower + Max_RSRP_gain - pathloss;
					if ( comp_sector_idx == -1 )
					{
						comp_signal_strength[1] = temp_signal_strength;

						// The interference induced from the comp sector to the signal of the main sector  
						comp_interf_strength[1] = bs_maxpower + 10*log10(channel_of_interest->signal_RSRP_gain
						[sector_idx][sec_z_idx_selected][sec_a_idx_selected]
						[zenith_angle_idx_selected_for_interference][azimuth_angle_idx_selected_for_interference][panel_idx_selected_for_interference])
						- pathloss;

						comp_azimuth_angle_idx_selected   = ue_a_idx_selected;
						comp_zenith_angle_idx_selected    = ue_z_idx_selected;
						comp_panel_idx_selected           = ue_p_idx_selected;
						comp_sector_a                     = sec_a_idx_selected;
						comp_sector_z                     = sec_z_idx_selected;
						comp_sector_idx                   = bs_idx * 3 + (sector_idx);						
					}
					if ( temp_signal_strength > comp_signal_strength[1] )
					{						
						Real temp_interf_strength = bs_maxpower + 10*log10(channel_of_interest->signal_RSRP_gain
						[sector_idx][sec_z_idx_selected][sec_a_idx_selected]
						[ue_z_idx_selected][ue_a_idx_selected][panel_idx_selected_for_interference])
						- pathloss;

						if (temp_signal_strength > temp_interf_strength)
						{
							comp_signal_strength[1] = temp_signal_strength;
							comp_interf_strength[1] = temp_interf_strength;
							comp_azimuth_angle_idx_selected   = ue_a_idx_selected;
							comp_zenith_angle_idx_selected    = ue_z_idx_selected;
							comp_panel_idx_selected           = ue_p_idx_selected;
							comp_sector_a                     = sec_a_idx_selected;
							comp_sector_z                     = sec_z_idx_selected;
							comp_sector_idx                   = bs_idx * 3 + (sector_idx);
						}
					}
				} // Exclude Serving BS & Sector 
			} // Sector index loop end
		} // BS index loop end
		#endif

		// The interference induced from the main sector to the signal of the comp sector
		comp_interf_strength[0] = bs_maxpower + 10*log10(
							channel[self_bs_idx][self_ms_idx].signal_RSRP_gain[self_sector_idx]
							[sector_zenith_angle_idx][sector_azimuth_angle_idx]
							[comp_azimuth_angle_idx_selected][comp_zenith_angle_idx_selected][comp_panel_idx_selected])
						  - channel[self_bs_idx][self_ms_idx].pathloss_final;
		//
		/////////////////////// Get interference ///////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			// 220811 jhnoh
			channel_of_interest = &channel[bs_idx][self_ms_idx];
			pathloss            = channel_of_interest->pathloss_final;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{

				if (Configuration_Type == 0)
				{
					// 
				}
				else
				{
					///////////////////////////////////////////////////////////////// random beam interference
					int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
					int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

					rand_sec_a[3*bs_idx+sector_idx] = aa;
					rand_sec_z[3*bs_idx+sector_idx] = zz;

					random_antgain = Get_antgain(channel_of_interest, bs_idx, sector_idx, aa, zz);   ///// MAX antgain - combined gain
					Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, sector_idx, aa, zz);
					random_RSRP_antgain = Get_RSRP(channel_of_interest, sector_idx, zz, aa, 1);
					/////////////////////////////////////////////////////////////////
				}

				signal = bs_maxpower + random_RSRP_antgain - pathloss;
				//
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////

				if ((bs_idx == self_bs_idx) && (sector_idx == self_sector_idx))
				{
					//
				}
				else
				{
					if (Configuration_Type == 0)  // config A
					{
						
					}
					else
					{
						interference += dBm2linear(signal);//pow(10., (signal) / 10.);
					}
				}

				if (self_ms_idx == 0)
				{
					intf_w_rnd_RSRP[bs_idx*3 + sector_idx] = signal;
				}
				else
				{
					intf_w_rnd_RSRP[bs_idx*3 + sector_idx] = signal;
				}
			}			
		}
	}
	//----------------------------------------------------------------------------------------------------------------------------//
	else if (TYPE == 12 || TYPE == 13)
	{
		int ms_to_bs_wrap_idx = (self_ms_idx / (3*num_MS_persector)); // -> ms_idx�� bs_idx��	
		for (int wrap_bs_idx = 0; wrap_bs_idx < simple_num_BS; wrap_bs_idx++)
		{
			int bs_idx                          = simple_wrap_mat[ms_to_bs_wrap_idx][wrap_bs_idx];
			channel_of_interest                 = &channel[bs_idx][self_ms_idx];
			pathloss                            = channel_of_interest->pathloss; //channel���� ����� pathloss�����´� // self_ms_idx�� ���� ����� 19��BS�� ���� ��
			RMS_delay                           = channel_of_interest->RMS_delay_spread;
			AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
			distance                            = channel_of_interest->distance;
			pathloss                            = pathloss + incar_loss + Otoi_loss + ms[self_ms_idx].LSPs[bs_idx](0);   ////// use stored LSP SF
			channel_of_interest->pathloss_final = pathloss;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
				{
					for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
					{
						antgain = Get_antgain(channel_of_interest, bs_idx, sector_idx, a, z);   ///// MAX antgain - combined gain
						Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, sector_idx, a, z);

						RSRP_antgain = Get_RSRP(channel_of_interest, sector_idx, z, a, 0);

						if (a == 0 && z == 0)  // first loop
						{
							MAX_RSRP_antgain = RSRP_antgain;

							azimuth_angle_idx_selected = selected_a;
							zenith_angle_idx_selected  = selected_z;
							panel_idx_selected         = selected_p;

							sector_a = a;
							sector_z = z;
						}
						else if (MAX_RSRP_antgain < RSRP_antgain)   // find max 
						{
							MAX_RSRP_antgain = RSRP_antgain;

							azimuth_angle_idx_selected = selected_a;
							zenith_angle_idx_selected  = selected_z;
							panel_idx_selected         = selected_p;

							sector_a = a;
							sector_z = z;
						}
						else
						{
							MAX_RSRP_antgain = MAX_RSRP_antgain;
						}
					}
				}

				signal = bs_maxpower + MAX_RSRP_antgain - pathloss;

				// 220808 jhnoh
				static_gain  [bs_idx * 3 + (sector_idx)].first  = signal;
				static_gain  [bs_idx * 3 + (sector_idx)].second = bs_idx * 3 + (sector_idx);

				analog_beam_selection[bs_idx * 3 + (sector_idx)].a = azimuth_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].z = zenith_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].p = panel_idx_selected;

				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_a = sector_a;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_z = sector_z;				

				only_pathloss_signal = bs_maxpower - pathloss;

				//if (bs_idx == 0 && sector_idx == 1)
				if (wrap_bs_idx == 0 && sector_idx == 0)
				{
					strongest_signal               = signal;
					self_bs_idx                    = bs_idx;
					self_sector_idx                = sector_idx;
					link_pathloss                  = pathloss;
					link_RSRP                      = MAX_RSRP_antgain;
					link_RMS_delay_spread          = RMS_delay;
					link_AOA_spread                = AOA_spread;
					link_AOD_spread                = AOD_spread;
					link_distance                  = distance; //channel_of_interest->distance;

					link_los                       = channel_of_interest->LOS;
					link_indoor                    = channel_of_interest->Indoor;

					azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
					zenith_angle_idx_selected_for_interference  = zenith_angle_idx_selected;
					panel_idx_selected_for_interference         = panel_idx_selected;

					sector_azimuth_angle_idx       = sector_a;
					sector_zenith_angle_idx        = sector_z;

					_sector_in_control             = self_bs_idx * 3 + (self_sector_idx);
				}
				else if (signal > strongest_signal)
				{
					// In single_cell_mode, only center BS (bs_idx 0) sectors can be serving sectors
					if (single_cell_mode == 1 && bs_idx != 0)
					{
						// Skip: this BS cannot be a serving BS in single_cell_mode
					}
					else
					{
						strongest_signal               = signal;
						self_bs_idx                    = bs_idx; // serving BS idx
						self_sector_idx                = sector_idx;
						link_pathloss                  = pathloss;
						link_RSRP                      = MAX_RSRP_antgain;
						link_RMS_delay_spread          = RMS_delay;
						link_AOA_spread                = AOA_spread;
						link_AOD_spread                = AOD_spread;
						link_distance                  = distance;

						link_los                       = channel_of_interest->LOS;
						link_indoor                    = channel_of_interest->Indoor;

						azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
						zenith_angle_idx_selected_for_interference  = zenith_angle_idx_selected;
						panel_idx_selected_for_interference         = panel_idx_selected;

						sector_azimuth_angle_idx = sector_a;
						sector_zenith_angle_idx  = sector_z;

						_sector_in_control = self_bs_idx * 3 + (self_sector_idx);
					}
				}
				else
				{
					// interference += dBm2linear(signal - MAX_RSRP_antgain + random_RSRP_antgain);
					// pow(10., (signal - MAX_RSRP_antgain + random_RSRP_antgain) / 10.);
					// jhnoh 220826
					// intf_w_rnd_RSRP[bs_idx*3 + sector_idx] = signal - MAX_RSRP_antgain + random_RSRP_antgain;
				}
			}			
		}

		// FOR mTRP mode 
		#if 0
		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2) 
		{
			for (int mTRP_idx = 0; mTRP_idx < num_mTRP; mTRP_idx++)
			{
				channel_of_interest                 = &channel[mTRP_idx+num_BS][self_ms_idx];
				pathloss                            = channel_of_interest->pathloss; //channel���� ����� pathloss�����´� // self_ms_idx�� ���� ����� 19��BS�� ���� ��
				RMS_delay                           = channel_of_interest->RMS_delay_spread;
				AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
				AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
				distance                            = channel_of_interest->distance;
				pathloss                            = pathloss + incar_loss + Otoi_loss + ms[self_ms_idx].LSPs[mTRP_idx+num_BS](0);   ////// use stored LSP SF
				channel_of_interest->pathloss_final = pathloss;

				for (int sectorized_panel_idx = 0; sectorized_panel_idx < 3; sectorized_panel_idx++)
				{
					for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
					{
						for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
						{
							antgain = Get_antgain(channel_of_interest, mTRP_idx+num_BS, sectorized_panel_idx, a, z);   ///// MAX antgain - combined gain
							Get_TX_SmallScale_antgain(channel_of_interest, mTRP_idx+num_BS, sectorized_panel_idx, a, z);

							RSRP_antgain = Get_RSRP(channel_of_interest, sectorized_panel_idx, z, a, 0);

							if (a == 0 && z == 0)  // first loop
							{
								MAX_RSRP_antgain = RSRP_antgain;

								azimuth_angle_idx_selected = selected_a;
								zenith_angle_idx_selected  = selected_z;
								panel_idx_selected         = selected_p;

								sector_a = a;
								sector_z = z;
							}
							else if (MAX_RSRP_antgain < RSRP_antgain)   // find max 
							{
								MAX_RSRP_antgain = RSRP_antgain;

								azimuth_angle_idx_selected = selected_a;
								zenith_angle_idx_selected  = selected_z;
								panel_idx_selected         = selected_p;

								sector_a = a;
								sector_z = z;
							}
							else
							{
								MAX_RSRP_antgain = MAX_RSRP_antgain;
							}
						}
					}

					signal = micro_bs_power + MAX_RSRP_antgain - pathloss;

					// 220808 jhnoh
					static_gain  [(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].first  = signal;
					static_gain  [(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].second = (mTRP_idx+num_BS) * 3 + (sectorized_panel_idx);

					analog_beam_selection[(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].a = azimuth_angle_idx_selected;
					analog_beam_selection[(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].z = zenith_angle_idx_selected;
					analog_beam_selection[(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].p = panel_idx_selected;

					analog_beam_selection[(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].sector_a = sector_a;
					analog_beam_selection[(mTRP_idx+num_BS) * 3 + (sectorized_panel_idx)].sector_z = sector_z;				

					only_pathloss_signal = micro_bs_power - pathloss;

					if (signal > strongest_signal)
					{
						strongest_signal               = signal;
						self_bs_idx                    = (mTRP_idx+num_BS); // serving BS idx
						self_sector_idx                = sectorized_panel_idx;
						link_pathloss                  = pathloss;
						link_RSRP                      = MAX_RSRP_antgain;
						link_RMS_delay_spread          = RMS_delay;
						link_AOA_spread                = AOA_spread;
						link_AOD_spread                = AOD_spread;
						link_distance                  = distance;

						link_los                       = channel_of_interest->LOS;
						link_indoor                    = channel_of_interest->Indoor;

						azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
						zenith_angle_idx_selected_for_interference  = zenith_angle_idx_selected;
						panel_idx_selected_for_interference         = panel_idx_selected;

						sector_azimuth_angle_idx = sector_a;
						sector_zenith_angle_idx  = sector_z;

						_sector_in_control = self_bs_idx * 3 + (self_sector_idx);
					}
					else
					{
						// interference += dBm2linear(signal - MAX_RSRP_antgain + random_RSRP_antgain);
						// pow(10., (signal - MAX_RSRP_antgain + random_RSRP_antgain) / 10.);
						// jhnoh 220826
						// intf_w_rnd_RSRP[bs_idx*3 + sector_idx] = signal - MAX_RSRP_antgain + random_RSRP_antgain;
					}
				}
			}
		}
		#endif		

		#if 0
		// Select Comp Sector
		int  comp_Rx_panel_idx            = 0;
		Real temp_signal_strength         = 0;

		comp_sector_idx                     = -1;
		comp_signal_strength[1]             = -INFINITY;
		comp_interf_strength[1]             = 0;

        comp_max_RSRP                       = 0;
        comp_inf_RSRP                       = 0;
        comp_azimuth_angle_idx_selected     = 0;
        comp_zenith_angle_idx_selected      = 0;
        comp_panel_idx_selected             = 0;
        comp_sector_a                       = 0;
        comp_sector_z                       = 0;

		if ( panel_idx_selected_for_interference == 0 )
		{
			comp_Rx_panel_idx = 1;
		}
		else
		{
			comp_Rx_panel_idx = 0;
		}

		for (int wrap_bs_idx = 0; wrap_bs_idx < simple_num_BS; wrap_bs_idx++)
		{
			int bs_idx = simple_wrap_mat[ms_to_bs_wrap_idx][wrap_bs_idx];
			// 220811 jhnoh
			channel_of_interest                 = &channel[bs_idx][self_ms_idx];
			pathloss                            = channel_of_interest->pathloss;
			RMS_delay                           = channel_of_interest->RMS_delay_spread;
			AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
			distance                            = channel_of_interest->distance;
			pathloss                            = channel_of_interest->pathloss_final;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				if ((bs_idx != self_bs_idx) || ( sector_idx != self_sector_idx))
				{			
					// Analog Beamforming 
					Real Tmp_RSRP_gain = 0;
					Real Max_RSRP_gain = 0;
					Real Inf_RSRP_gain = 0;

					int ue_a_idx_selected  = -1;
					int ue_z_idx_selected  = -1;
					int ue_p_idx_selected  = -1;
					int sec_a_idx_selected = -1;
					int sec_z_idx_selected = -1;

					// Sector beam index loop
					for (int sec_a_idx = 0; sec_a_idx < tilt_azimuth_angle_LCS_size; sec_a_idx++)
					{
						for (int sec_z_idx = 0; sec_z_idx < tilt_zenith_angle_LCS_size; sec_z_idx++)
						{
							// UE beam index loop
							if (ue_antenna_element_gain == 0)
							{
								Tmp_RSRP_gain = 10*log10(channel_of_interest->signal_RSRP_gain
												[sector_idx][sec_z_idx][sec_a_idx]
												[0][0][0]);
								if ( Tmp_RSRP_gain > Max_RSRP_gain  || 
								( sec_a_idx == 0 && sec_z_idx == 0 ))
								{
									Max_RSRP_gain       = Tmp_RSRP_gain;
									ue_a_idx_selected   = 0;
									ue_z_idx_selected   = 0;
									ue_p_idx_selected   = 0;
									sec_a_idx_selected  = sec_a_idx;
									sec_z_idx_selected  = sec_z_idx;
								}
							}
							else
							{
								for (int ue_z_idx = 0; ue_z_idx < ue_tilt_zenith_angle_LCS_size; ue_z_idx++)
								{
									for (int ue_a_idx = 0; ue_a_idx < ue_tilt_azimuth_angle_LCS_size; ue_a_idx++)
									{						
										Tmp_RSRP_gain = 10*log10(channel_of_interest->signal_RSRP_gain
																[sector_idx][sec_z_idx][sec_a_idx]
																[ue_z_idx][ue_a_idx][comp_Rx_panel_idx]);


										if ( Tmp_RSRP_gain > Max_RSRP_gain  || 
										( sec_a_idx == 0 && sec_z_idx == 0 && ue_z_idx == 0 && ue_a_idx == 0 ))
										{
											Max_RSRP_gain       = Tmp_RSRP_gain;
											ue_a_idx_selected   = ue_a_idx;
											ue_z_idx_selected   = ue_z_idx;
											ue_p_idx_selected   = comp_Rx_panel_idx;
											sec_a_idx_selected  = sec_a_idx;
											sec_z_idx_selected  = sec_z_idx;
										}
									}
								} // UE beam index loop end
							}
						}
					} // Sector beam index loop end

					temp_signal_strength = bs_maxpower + Max_RSRP_gain - pathloss;
					if ( comp_sector_idx == -1 )
					{
						comp_signal_strength[1] = temp_signal_strength;

						// The interference induced from the comp sector to the signal of the main sector  
						comp_interf_strength[1] = bs_maxpower + 10*log10(channel_of_interest->signal_RSRP_gain
						[sector_idx][sec_z_idx_selected][sec_a_idx_selected]
						[zenith_angle_idx_selected_for_interference][azimuth_angle_idx_selected_for_interference][panel_idx_selected_for_interference])
						- pathloss;

						comp_azimuth_angle_idx_selected   = ue_a_idx_selected;
						comp_zenith_angle_idx_selected    = ue_z_idx_selected;
						comp_panel_idx_selected           = ue_p_idx_selected;
						comp_sector_a                     = sec_a_idx_selected;
						comp_sector_z                     = sec_z_idx_selected;
						comp_sector_idx                   = bs_idx * 3 + (sector_idx);						
					}
					if ( temp_signal_strength > comp_signal_strength[1] )
					{						
						Real temp_interf_strength = bs_maxpower + 10*log10(channel_of_interest->signal_RSRP_gain
						[sector_idx][sec_z_idx_selected][sec_a_idx_selected]
						[ue_z_idx_selected][ue_a_idx_selected][panel_idx_selected_for_interference])
						- pathloss;

						if (temp_signal_strength > temp_interf_strength)
						{
							comp_signal_strength[1] = temp_signal_strength;
							comp_interf_strength[1] = temp_interf_strength;
							comp_azimuth_angle_idx_selected   = ue_a_idx_selected;
							comp_zenith_angle_idx_selected    = ue_z_idx_selected;
							comp_panel_idx_selected           = ue_p_idx_selected;
							comp_sector_a                     = sec_a_idx_selected;
							comp_sector_z                     = sec_z_idx_selected;
							comp_sector_idx                   = bs_idx * 3 + (sector_idx);
						}
					}
				} // Exclude Serving BS & Sector 
			} // Sector index loop end
		} // BS index loop end
		

		// The interference induced from the main sector to the signal of the comp sector
		comp_interf_strength[0] = micro_bs_power + 10*log10(
							channel[self_bs_idx][self_ms_idx].signal_RSRP_gain[self_sector_idx]
							[sector_zenith_angle_idx][sector_azimuth_angle_idx]
							[comp_azimuth_angle_idx_selected][comp_zenith_angle_idx_selected][comp_panel_idx_selected])
						  - channel[self_bs_idx][self_ms_idx].pathloss_final;
		#endif

		//-------------------------------------------------------------------------
		//---------------------------  Get interference ---------------------------
		//-------------------------------------------------------------------------
		for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		{
			pathloss = channel[bs_idx][self_ms_idx].pathloss_final;
			channel_of_interest = &channel[bs_idx][self_ms_idx];

			intf_w_rnd_RSRP[bs_idx*3  ] = linear2dBm(0);
			intf_w_rnd_RSRP[bs_idx*3+1] = linear2dBm(0);
			intf_w_rnd_RSRP[bs_idx*3+2] = linear2dBm(0);

			rand_sec_a[3*bs_idx  ] = -1;
			rand_sec_a[3*bs_idx+1] = -1;
			rand_sec_a[3*bs_idx+2] = -1;
			
			rand_sec_z[3*bs_idx  ] = -1;
			rand_sec_z[3*bs_idx+1] = -1;
			rand_sec_z[3*bs_idx+2] = -1;

			//int sector_idx = (int)(3*randnum.u());
			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				if ((bs_idx == self_bs_idx) && (sector_idx == self_sector_idx))
				{
				}
				else
				{
					//--------------------------- random beam interference ---------------------------
					int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
					int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

					rand_sec_a[3*bs_idx+sector_idx] = aa;
					rand_sec_z[3*bs_idx+sector_idx] = zz;

					//--------------------------------------------------------------------------------				
					random_RSRP_antgain = linear2dB(channel[bs_idx][self_ms_idx].signal_RSRP_gain[sector_idx][zz][aa]
									[zenith_angle_idx_selected_for_interference]
									[azimuth_angle_idx_selected_for_interference]
									[panel_idx_selected_for_interference]);
					
					signal = bs_maxpower + random_RSRP_antgain - pathloss;

					interference += dBm2linear(signal);
					intf_w_rnd_RSRP[bs_idx*3 + sector_idx] = signal;
				}				
			}
		}

	}
	#if 0
	else if (TYPE == 13)
	{
		int ms_to_bs_wrap_idx = (self_ms_idx / 30); // -> ms_idx�� bs_idx��


		//for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		for (int wrap_bs_idx = 0; wrap_bs_idx < simple_num_BS; wrap_bs_idx++)
		{
			int bs_idx = simple_wrap_mat[ms_to_bs_wrap_idx][wrap_bs_idx];
			channel_of_interest = &channel[bs_idx][self_ms_idx];

			pathloss   = channel_of_interest->pathloss; //channel���� ����� pathloss�����´� // self_ms_idx�� ���� ����� 19��BS�� ���� ��
			RMS_delay  = channel_of_interest->RMS_delay_spread;
			AOA_spread = channel_of_interest->circular_angle_spread_AOA;
			AOD_spread = channel_of_interest->circular_angle_spread_AOD;
			distance   = channel_of_interest->distance;
			//cout << distance << endl;

			pathloss = pathloss + incar_loss + Otoi_loss + ms[self_ms_idx].LSPs[bs_idx](0);   ////// use stored LSP SF

			channel_of_interest->pathloss_final = pathloss;

			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{

				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////   Antenna gain - NEW version //////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////				
				for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
				{
					for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
					{
						antgain = Get_antgain(channel_of_interest, bs_idx, sector_idx, a, z);   ///// MAX antgain - combined gain
						Get_TX_SmallScale_antgain(channel_of_interest, bs_idx, sector_idx, a, z);

						RSRP_antgain = Get_RSRP(channel_of_interest, (sector_idx), z, a, 0);

						MAX_RSRP_antgain = RSRP_antgain;



						azimuth_angle_idx_selected = selected_a;
						zenith_angle_idx_selected = selected_z;
						panel_idx_selected = selected_p;

						sector_a = a;
						sector_z = z;


					}
				}
				
				// 220808 jhnoh
				static_gain  [bs_idx * 3 + (sector_idx)].first  = signal;
				static_gain  [bs_idx * 3 + (sector_idx)].second = bs_idx * 3 + (sector_idx);

				analog_beam_selection[bs_idx * 3 + (sector_idx)].a = azimuth_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].z = zenith_angle_idx_selected;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].p = panel_idx_selected;

				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_a = sector_a;
				analog_beam_selection[bs_idx * 3 + (sector_idx)].sector_z = sector_z;	

				//////////////////////////////////////////////////////////////////////////////////////////////////////////

				//if (bs_idx == 0 && sector_idx == 1)
				if (wrap_bs_idx == 0 && sector_idx == 0)
				{

					strongest_signal = signal;
					link_prev_MAX_RSRP = MAX_RSRP_antgain;

					self_bs_idx = bs_idx;
					self_sector_idx = sector_idx;
					link_pathloss = pathloss;
					link_coupling_loss = pathloss - antgain;
					//link_pathloss = pathloss - ms[self_ms_idx].LSPs[bs_idx](0);
					link_antgain = antgain;
					serving_bs = channel_of_interest->self_bs;
					link_RMS_delay_spread = RMS_delay;
					link_AOA_spread = AOA_spread;
					link_AOD_spread = AOD_spread;
					link_distance = channel_of_interest->distance;

					azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
					zenith_angle_idx_selected_for_interference = zenith_angle_idx_selected;
					panel_idx_selected_for_interference = panel_idx_selected;

					sector_azimuth_angle_idx = sector_a;
					sector_zenith_angle_idx = sector_z;

					_sector_in_control = self_bs_idx * 3 + (self_sector_idx);
				}
				else if (signal > strongest_signal)
				{
					interference += dBm2linear(strongest_signal); //pow(10., (strongest_signal) / 10.);
					// jhnoh 220826
					intf_w_rnd_RSRP[wrap_bs_idx*3 + sector_idx] = strongest_signal;

					strongest_signal = signal;
					link_prev_MAX_RSRP = MAX_RSRP_antgain;

					self_bs_idx = bs_idx; // serving BS idx
					self_sector_idx = sector_idx;
					link_pathloss = pathloss;
					link_coupling_loss = pathloss - antgain;
					//link_pathloss = pathloss - ms[self_ms_idx].LSPs[bs_idx](0);
					link_antgain = antgain;
					serving_bs = channel_of_interest->self_bs;
					link_RMS_delay_spread = RMS_delay;
					link_AOA_spread = AOA_spread;
					link_AOD_spread = AOD_spread;
					link_distance = channel_of_interest->distance;

					azimuth_angle_idx_selected_for_interference = azimuth_angle_idx_selected;
					zenith_angle_idx_selected_for_interference = zenith_angle_idx_selected;
					panel_idx_selected_for_interference = panel_idx_selected;

					sector_azimuth_angle_idx = sector_a;
					sector_zenith_angle_idx = sector_z;

					_sector_in_control = self_bs_idx * 3 + (self_sector_idx);
				}
				else
				{
					interference += dBm2linear(signal); //pow(10., (signal) / 10.);
					intf_w_rnd_RSRP[wrap_bs_idx*3 + sector_idx] = signal;
				}
			}
		}
	}
	#endif
	str_signal   = strongest_signal;
	signal       = strongest_signal;
	interference = linear2dBm(interference); // 10. * log10(interference);
	Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
	geometry = dBm2linear(signal)/(dBm2linear(interference) + noise);
}
#endif  // deprecated Get_signal_interference_old

void LINK::Get_adj_SECTORS()
{
	Real pathlosss;
	Real antgainn;

	sort(static_gain,static_gain + num_SECTORS, comparator);


	for (int idx = 0; idx < num_SECTORS; idx++)
		adj_sector[idx] = static_gain[idx].second;

	// jhnoh 230109 
	// Consider that ue select the other Rx antenna pannel for avoiding inter-layer interference
	// adj_sector[1] = comp_sector_idx;
}


void LINK::UE_Initial_Setting(void)
{
  // Deterministic per-UE RNG (reproducibility fix). UE_Initial_Setting runs inside
  // the OpenMP "Get RSRP" loop (Initiallization.cpp), where the O2I / in-car loss
  // draws below would otherwise data-race the global `randnum` across threads AND
  // consume it by a run-dependent amount (work distribution varies per run) — the
  // dominant source of non-reproducible coupling loss / SINR. Shadowing `randnum`
  // with a generator seeded from this UE's identity (base_seed, drop, ue) makes the
  // O2I realization a pure function of the work item: bit-reproducible across runs
  // and thread counts. Draws below are sequential, so their order is fixed.
  Rand randnum = make_workitem_rng((unsigned long long)_seed,
                                   (unsigned long long)drop_idx,
                                   (unsigned long long)self_ms_idx,
                                   0x4F32491ULL /* salt: O2I stream */);

  /*------------------------ Car penetration loss ----------------------------*/
  if (TYPE == 12) // Dense Urban 
  {
	  if (ms[self_ms_idx].Indoor == true)   ///// Indoor
	  {
		  incar_loss = 0.;
	  }
	  else if (ms[self_ms_idx].Indoor == false)                              ////// outdoor
	  {
		  // 3GPP TR 38.901 Table A-2: UMa scenario assumes outdoor UEs are pedestrians
		  // Pedestrian UEs have no in-car penetration loss
		  incar_loss = 0.;
	  }
	  else
	  {
		  cout << "Channel Model Type ERROR!! - in_car_loss 1111 !!" << endl;
		  cout << "Indoor type for ms[" << self_ms_idx <<"] is "<<  ms[self_ms_idx].Indoor << endl;
		  getchar();
	  }
  }
  else if (TYPE == 13) // Rural
  {
	  if (ms[self_ms_idx].Indoor == true)   ///// Indoor
	  {
		  incar_loss = 0.;
	  }
	  else                              ////// outdoor
	  {
		  if (Configuration_Type == 2)     /////// Rural configuration C = 40% outdoor pedestrian, 20% outdoor incar
		  {
			  if (randnum.u() < (0.4/0.6))  //// Pedestrian
			  {
				  incar_loss = 0.;
			  }
			  else                   //// In Car
			  {
				  incar_loss = 9. + 5. * randnum.n();
			  }
		  }
		  else    /// 50% All outdoor
		  {
			  incar_loss = 9. + 5. * randnum.n();
		  }
	  }
  }
  else //// Indoor case
  {
    incar_loss = 0.;
  }


  ///////////////////////////// Building penetration loss /////////////////////////////////
  Real glass_L    = 2 + 0.2 * carrier_freq / (1000000000.);      ///// Standard multi-pane glass , carrier_freq = GHz
  Real IRRglass_L = 23 + 0.3 * carrier_freq / (1000000000.);     ///// Infrared Reflective(IRR) glass , carrier_freq = GHz
  Real concrete_L = 5 + 4 * carrier_freq / (1000000000.);        ///// Concrete glass , carrier_freq = GHz
  Real wood_L     = 4.85 + 0.12 * carrier_freq / (1000000000.);  ///// Wood glass , carrier_freq = GHz

  Real indoor_L;

  if (TYPE == 12) // Dense Urban UMa
  {
	  if (ms[self_ms_idx].Indoor == false) // Outdoor
	  {
		  Otoi_loss = 0;
	  }
	  else // Indoor
	  {
		  if (Channel_Model_Type == 1) // Model B
		  {
			  if (randnum.u() > 0.5) // 80% Low Loss  O2I
			  {
				  // Low-loss model
				  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
				  indoor_L  = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));
				  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();
			  }
			  else // 20& High Loss O2I
			  {
				  // High-loss model
				  Otoi_loss = 5 - 10 * log10(0.7 * pow(10, -1 * (IRRglass_L / 10)) + 0.3 * pow(10, -1 * (concrete_L / 10)));
				  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));

				  Otoi_loss = Otoi_loss + indoor_L + 6.5 * randnum.n();
			  }
		  }
		  else if (Channel_Model_Type == 0 || Channel_Model_Type == 2)  //// Model A 
		  {
			  if (carrier_freq >= 6000000000.)  //// above 6GHz
			  {
				  if (randnum.u() > 0.2) 
				  // 3GPP TR 38.901 50% Low Loss 
				  // ITU-R M.2412 80% Low Loss  O2I
				  {
					  // Low-loss model
					  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
					  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));
					  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();
					  high_loss_flag = 0;
				  }
				  else
				  {
					  // High-loss model
					  Otoi_loss = 5 - 10 * log10(0.7 * pow(10, -1 * (IRRglass_L / 10)) + 0.3 * pow(10, -1 * (concrete_L / 10)));
					  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));
					  Otoi_loss = Otoi_loss + indoor_L + 6.5 * randnum.n();
					  high_loss_flag = 1;
				  }
			  }
			  else // below 6GHz
			  {
				  Otoi_loss = 20.;
				  indoor_L = 0.5 * ms[self_ms_idx].ms_d_in;
				  Otoi_loss = Otoi_loss + indoor_L;
			  }
		  }
		  else
		  {
			  cout << "Channel Model Type ERROR!! - Otoi Loss" << endl;
		  }
	  }
  }
  else if (TYPE == 13) // Rural
  {
	  if (ms[self_ms_idx].Indoor == false)   ///// Outdoor
	  {
		  Otoi_loss = 0;
	  }
	  else                               /////// Indoor
	  {
		  if (Channel_Model_Type == 1)   ///// Model B
		  {
			  /////// 100% Low-loss model
			  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
			  indoor_L = 0.5 * MIN((10 * randnum.u()), (10 * randnum.u()));

			  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();

		  }
		  else if (Channel_Model_Type == 0)  //// Model A 
		  {
			  if (carrier_freq > 6000000000)  //// above 6GHz
			  {
				  /////// 100% Low-loss model
				  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
				  indoor_L = 0.5 * MIN((10 * randnum.u()), (10 * randnum.u()));

				  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();
			  }
			  else                         //// below 6GHz
			  {
				  Otoi_loss = 10;
				  indoor_L = 0.5 * (10 * randnum.u());

				  Otoi_loss = Otoi_loss + indoor_L;
			  }
		  }
		  else
		  {
			  cout << "Channel Model Type ERROR!! - Otoi Loss" << endl;
		  }
	  }
  }
  else   ///// Indoor case
  {
	  Otoi_loss = 0.;
  }


  // initial setting
  _throughput = 0 ;
  CQI_offset = 0;
  num_added = 0;
  num_re_tx =    0;
  ACK = true;

  //H_m = NULL;
  //HI_m = NULL;
}


Real LINK::Get_angle(LOCATION interferer, LOCATION victim, LOCATION interfere_rx  )
{
	Real angle = 0;
	Real vector1_x, vector1_y, vector2_x, vector2_y;
	Real acos_input = 0;

	vector1_x = victim.x - interferer.x;  /// BS->MS vector
	vector1_y = victim.y - interferer.y;
	vector2_x = interfere_rx.x - interferer.x;  /// Boresight vector
	vector2_y = interfere_rx.y - interferer.y;
	acos_input = (vector1_x*vector2_x + vector1_y*vector2_y)
		/ (sqrt(vector1_x*vector1_x + vector1_y*vector1_y)*sqrt(vector2_x*vector2_x + vector2_y*vector2_y));

	if (acos_input >= 1)
	{
		angle = 0;
	}
	else
	{
		angle = acos(acos_input);

		if ((victim.y-interfere_rx.y)<0.)
		{
			angle = -angle;
		}

	}

	return angle * 180. / pi;

}

Real LINK::Get_antgain(CHANNEL * channel_of_interest, int _bs_idx, int sector_idx, int tilt_a_idx, int tilt_z_idx)  ////// LOS gain BS
{
	Real h_angle_pi     = channel_of_interest->LOS_AOD_GCS * (pi / 180.);  //// pi in GCS
	Real v_angle_theta  = channel_of_interest->LOS_ZOD_GCS * (pi / 180.);  ///// theta in GCS

	Real F_theta_GCS_P1 = 0;
	Real F_pi_GCS_P1    = 0;
	Real F_theta_GCS_P2 = 0;
	Real F_pi_GCS_P2    = 0;

	int _sector_idx       = sector_idx;

	ComplexReal  F_tx_theta  = {0, 0};
	ComplexReal  F_tx_pi     = {0, 0};

	Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, _sector_idx, 
	F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

	// SNS: pass LOS VR info to field pattern for per-element attenuation
	const ClusterVR* los_vr_ptr = (channel_of_interest->sns_any_limited) ? &channel_of_interest->sns_vr_los : nullptr;
	Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, _sector_idx, tilt_z_idx, tilt_a_idx,
	v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi, los_vr_ptr, g_sns_rolloff_C);

	TransmitterAntennaGainXLOS_theta = F_tx_theta;
	TransmitterAntennaGainXLOS_pi    = F_tx_pi;

	return combined_antenna_gain;
}


void LINK::Get_TX_SmallScale_antgain(CHANNEL * channel_of_interest, int _bs_idx, int sector_idx, int tilt_a_idx, int tilt_z_idx)
{
	if (TYPE == 11)   /// Indoor
	{
		SmallScale_TX_AntennaGainXLOS_theta.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_TX_AntennaGainXLOS_pi.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);

		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < MAX_NUM_RAYS; j++)
			{
				SmallScale_TX_AntennaGainXLOS_theta(i, j) = 0.;
				SmallScale_TX_AntennaGainXLOS_pi(i, j) = 0.;
			}
		}


		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
			{

				Real h_angle_pi    = channel_of_interest->ray_AOD[i][j][0] * (pi / 180.);
				Real v_angle_theta = channel_of_interest->ray_AOD[i][j][1] * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_tx_theta = {0, 0};
				ComplexReal  F_tx_pi = {0, 0};

				Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

				// SNS: pass per-cluster VR info to field pattern for per-element attenuation
				{
					const ClusterVR* vr_ptr = (channel_of_interest->sns_any_limited) ? &channel_of_interest->sns_vr[i] : nullptr;
					Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi, vr_ptr, g_sns_rolloff_C);
				}


				SmallScale_TX_AntennaGainXLOS_theta(i, j) = F_tx_theta;
				SmallScale_TX_AntennaGainXLOS_pi(i, j) = F_tx_pi;
			}
		}
	}
	else if (TYPE == 12 || TYPE == 13)   ////// Dense Urban & Rural
	{
		SmallScale_TX_AntennaGainXLOS_theta.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_TX_AntennaGainXLOS_pi.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);

		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < MAX_NUM_RAYS; j++)
			{
				SmallScale_TX_AntennaGainXLOS_theta(i, j) = 0.;
				SmallScale_TX_AntennaGainXLOS_pi(i, j) = 0.;
			}
		}


		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
			{
				Real h_angle_pi    = channel_of_interest->ray_AOD[i][j][0] * (pi / 180.);
				Real v_angle_theta = channel_of_interest->ray_AOD[i][j][1] * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_tx_theta = {0, 0};
				ComplexReal  F_tx_pi = {0, 0};

				Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

				// SNS: pass per-cluster VR info to field pattern for per-element attenuation
				{
					const ClusterVR* vr_ptr = (channel_of_interest->sns_any_limited) ? &channel_of_interest->sns_vr[i] : nullptr;
					Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi, vr_ptr, g_sns_rolloff_C);
				}

				SmallScale_TX_AntennaGainXLOS_theta(i,j) = F_tx_theta;
				SmallScale_TX_AntennaGainXLOS_pi(i,j) = F_tx_pi;

			}
		}
	}

	else
	{
		cout << "TX small scale antgain ERROR!" << endl;
		getchar();
	}

}

void LINK::Get_RX_SmallScale_antgain(int M, int N, int P, CHANNEL * channel_of_interest, int sector_idx, int tilt_z_idx, int tilt_a_idx, int ue_interference_index)
{
	int ms_idx = self_ms_idx;

	if (TYPE == 11)   /// Indoor
	{
		SmallScale_RX_AntennaGainXLOS_theta.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_pi.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_theta_panel_2.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_pi_panel_2.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);

		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < MAX_NUM_RAYS; j++)
			{
				SmallScale_RX_AntennaGainXLOS_theta(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_pi(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j) = 0.;
			}
		}


		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
			{
				Real h_angle_pi = channel_of_interest->ray_AOA[i][j][0] * (pi / 180.);
				Real v_angle_theta = channel_of_interest->ray_AOA[i][j][1] * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_rx_theta   = {0, 0};
				ComplexReal  F_rx_pi      = {0, 0};
				ComplexReal  F_rx_theta_2 = {0, 0};
				ComplexReal  F_rx_pi_2    = {0, 0};

				Real combined_antenna_gain;


				if (ue_antenna_element_gain == 0 && !handheld_mode) // UE omni-antenna case
				{
					combined_antenna_gain = 0;

					if (P == 0)
					{
						SmallScale_RX_AntennaGainXLOS_theta(i, j) = 1.;
						SmallScale_RX_AntennaGainXLOS_pi(i, j) = 0.;
					}
					else
					{
						SmallScale_RX_AntennaGainXLOS_theta(i, j) = 0.;
						SmallScale_RX_AntennaGainXLOS_pi(i, j) = 1.;
					}
				}
				else
				{

					combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, handheld_mode ? M : -1);

					Get_UE_antenna_field_pattern(channel_of_interest->self_ms, M, N, P, ms_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_rx_theta, F_rx_pi, F_rx_theta_2, F_rx_pi_2);

					SmallScale_RX_AntennaGainXLOS_theta(i, j) = F_rx_theta;
					SmallScale_RX_AntennaGainXLOS_pi(i, j) = F_rx_pi;

					SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) = F_rx_theta_2;
					SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j) = F_rx_pi_2;
				}


			}
		}
	}
	else if (TYPE == 12 || TYPE == 13)   ////// Dense Urban & Rural
	{

		SmallScale_RX_AntennaGainXLOS_theta.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_pi.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_theta_panel_2.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);
		SmallScale_RX_AntennaGainXLOS_pi_panel_2.resize(channel_of_interest->NUM_PATH_for_channelcoeff, MAX_NUM_RAYS);

		// initialize
		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < MAX_NUM_RAYS; j++)
			{
				SmallScale_RX_AntennaGainXLOS_theta(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_pi(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) = 0.;
				SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j) = 0.;
			}
		}

		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
			{
				Real h_angle_pi = channel_of_interest->ray_AOA[i][j][0] * (pi / 180.);
				Real v_angle_theta = channel_of_interest->ray_AOA[i][j][1] * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_rx_theta   = {0, 0};
				ComplexReal  F_rx_pi      = {0, 0};
				ComplexReal  F_rx_theta_2 = {0, 0};
				ComplexReal  F_rx_pi_2    = {0, 0};

				Real combined_antenna_gain;

				if (ue_antenna_element_gain == 0 && !handheld_mode) // UE omni-antenna case
				{
					combined_antenna_gain = 0;

					/// UE antenna = omni case
					if (P == 0)
					{
						SmallScale_RX_AntennaGainXLOS_theta(i, j) = 1.;
						SmallScale_RX_AntennaGainXLOS_pi(i, j) = 0.;
					}
					else
					{
						SmallScale_RX_AntennaGainXLOS_theta(i, j) = 0.;
						SmallScale_RX_AntennaGainXLOS_pi(i, j) = 1.;
					}


				}
				else  // UE directional antenna case (or handheld)
				{

					combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, handheld_mode ? M : -1);

					Get_UE_antenna_field_pattern(channel_of_interest->self_ms, M, N, P, ms_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_rx_theta, F_rx_pi, F_rx_theta_2, F_rx_pi_2);

					SmallScale_RX_AntennaGainXLOS_theta(i, j) = F_rx_theta;
					SmallScale_RX_AntennaGainXLOS_pi(i, j) = F_rx_pi;

					SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) = F_rx_theta_2;
					SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j) = F_rx_pi_2;
				}
			}
		}
	}

	else
	{
		cout << "RX small scale antgain ERROR!" << endl;
		getchar();
	}
}


Real LINK::Get_RSRP(CHANNEL * channel_of_interest, int sec_number, int sec_z_idx, int sec_a_idx, int ue_interference)
{
	// ue_interference = 0 -> find max ue beam & RSRP
	// ue_interference = 1 -> for interference calculation, ue beam already selected
	complex<Real> Big_PI_LOS             (0. , channel_of_interest->random_phase_vv_LOS * (pi / 180.0));   // use stored initial phase

	// 220815 jhnoh
    //CHANNEL * channel_of_interest = &channel[_bs_idx][_ms_idx];
	Real K_linear = channel_of_interest->K_linear;

	Real XPR = pow(10.0, (channel_of_interest->XPR) / 10.0);

	Real alpha = 0;
	Real alpha_panel_2 = 0;
	Real alpha_zero = 0.;
	Real alpha_zero_panel_2 = 0.;
	ComplexReal  alpha_zero_temp = 0.;
	ComplexReal  alpha_zero_temp_panel_2 = 0.;
	Real alpha_nmup = 0.;
	Real alpha_nmup_panel_2 = 0.;
	ComplexReal  alpha_nmup_temp = 0.;
	ComplexReal  alpha_nmup_temp_panel_2 = 0.;


	Real UE_antenna_gain;

	Real kappa;
	Real _1_over_sqrt_K;

	if (ue_antenna_element_gain == 0)  //// 1 to 1 mapping, NOT USE ue beam, ue PORT
	{
		int M = MS_M;
		int N = MS_N;
		int P = MS_P;

		for (int m = 0; m < M; m++)
		{
			for (int n = 0; n < N; n++)
			{
				for (int p = 0; p < P; p++)
				{
					alpha_zero = 0;
					alpha_nmup = 0;

					UE_antenna_gain = Get_MS_antgain(m, n, p, channel_of_interest, 0, 0, 0, 0);   ///// MAX antgain - combined gain;
					Get_RX_SmallScale_antgain(m, n, p, channel_of_interest, 0, 0, 0, 0);

					//if (p == 0)  // polar 1
					//{
					if (channel_of_interest->Propagation == LOS_propagation)
					{
						alpha_zero_temp = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);

						alpha_zero = abs(alpha_zero_temp)*abs(alpha_zero_temp);
						// SNS: per-element attenuation now applied inside Get_BS_antenna_field_pattern
					}


					for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
					{
						Real cluster_power = 0.;
						for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
						{
							kappa = channel_of_interest->kappa[i][j];
							_1_over_sqrt_K = 1.0 / sqrt(kappa);

							complex<Real> Big_pi_NLOS_thetatheta(0, channel_of_interest->random_phase_vv[i][j] * (pi / 180.0));   // use stored initial phase
							complex<Real> Big_pi_NLOS_thetapi(0, channel_of_interest->random_phase_vh[i][j] * (pi / 180.0));
							complex<Real> Big_pi_NLOS_pitheta(0, channel_of_interest->random_phase_hv[i][j] * (pi / 180.0));
							complex<Real> Big_pi_NLOS_pipi(0, channel_of_interest->random_phase_hh[i][j] * (pi / 180.0));

							alpha_nmup_temp = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
								(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
								(SmallScale_RX_AntennaGainXLOS_theta(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
								);

							cluster_power += (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));
						}
						// SNS: per-element attenuation now applied inside Get_BS_antenna_field_pattern
						alpha_nmup += cluster_power;
					}
					alpha += alpha_zero + alpha_nmup;
				}
			}
		}
		channel_of_interest->signal_RSRP_gain[sec_number][sec_z_idx][sec_a_idx][0][0][0] = alpha / (Real)(M*N*P); // *P

		selected_a = 0;
		selected_z = 0;
		selected_p = 0;

		Real RSRP_antenna_gain = (alpha) / (M*N*P); // *P

		return 10 * log10(RSRP_antenna_gain);
			
	}

	///// UE beamforming 
	else  ///// UE beamforming 
	{
		if (ue_interference == 1)  /// for interference calculation, ue beam already selected
		{
			int M = MS_Mp;
			int N = MS_Np;
			int P = MS_P;
			Real MAX_alpha;

			alpha = 0;
			alpha_panel_2 = 0;


			for (int m = 0; m < M; m++)
			{
				for (int n = 0; n < N; n++)
				{
					for (int polar = 0; polar < P; polar++)
					{
						alpha_zero = 0;
						alpha_nmup = 0;
						alpha_zero_panel_2 = 0;
						alpha_nmup_panel_2 = 0;

						UE_antenna_gain = Get_MS_antgain(m, n, polar, channel_of_interest, 0, zenith_angle_idx_selected_for_interference, azimuth_angle_idx_selected_for_interference, ue_interference);   ///// MAX antgain - combined gain;
						Get_RX_SmallScale_antgain(m, n, polar, channel_of_interest, 0, zenith_angle_idx_selected_for_interference, azimuth_angle_idx_selected_for_interference, ue_interference);

						if (channel_of_interest->Propagation == LOS_propagation)
						{
							// panel 1
							alpha_zero_temp = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
							alpha_zero = abs(alpha_zero_temp)*abs(alpha_zero_temp);
							// SNS: per-element attenuation now in Get_BS_antenna_field_pattern

							// panel 2
							alpha_zero_temp_panel_2 = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
							alpha_zero_panel_2 = abs(alpha_zero_temp_panel_2)*abs(alpha_zero_temp_panel_2);
							// SNS: per-element attenuation now in Get_BS_antenna_field_pattern
						}

						for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
						{
							Real cluster_power_p1 = 0.;
							Real cluster_power_p2 = 0.;
							for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
							{
								kappa = channel_of_interest->kappa[i][j];
								_1_over_sqrt_K = 1.0 / sqrt(kappa);

								complex<Real> Big_pi_NLOS_thetatheta(0, channel_of_interest->random_phase_vv[i][j] * (pi / 180.0));   // use stored initial phase
								complex<Real> Big_pi_NLOS_thetapi(0, channel_of_interest->random_phase_vh[i][j] * (pi / 180.0));
								complex<Real> Big_pi_NLOS_pitheta(0, channel_of_interest->random_phase_hv[i][j] * (pi / 180.0));
								complex<Real> Big_pi_NLOS_pipi(0, channel_of_interest->random_phase_hh[i][j] * (pi / 180.0));

								// panel 1
								alpha_nmup_temp = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
									(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
									(SmallScale_RX_AntennaGainXLOS_theta(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
									);

								cluster_power_p1 += (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));

								/// panel 2
								alpha_nmup_temp_panel_2 = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
									(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
									(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
									);

								cluster_power_p2 += (abs(alpha_nmup_temp_panel_2) * abs(alpha_nmup_temp_panel_2));

							}
							alpha_nmup += cluster_power_p1;  // SNS: per-element in Get_BS_antenna_field_pattern
							alpha_nmup_panel_2 += cluster_power_p2;  // SNS: per-element in Get_BS_antenna_field_pattern
						}
						alpha += alpha_zero + alpha_nmup;
						alpha_panel_2 += alpha_zero_panel_2 + alpha_nmup_panel_2;
					}
				}
			}
			Real alpha_interference;

			if (panel_idx_selected_for_interference == 0)
			{
				alpha_interference = alpha;
			}
			else if (panel_idx_selected_for_interference == 1)
			{
				alpha_interference = alpha_panel_2;
			}

			Real RSRP_antenna_gain = (alpha_interference) / (M*N*P);

			return 10 * log10(RSRP_antenna_gain);
		}
		else  /// find max ue beam selection
		{
			int M = MS_Mp;
			int N = MS_Np;
			int P = MS_P;
			Real MAX_alpha;

			for (int z = 0; z < ue_tilt_zenith_angle_LCS_size; z++)
			{
				for (int a = 0; a < ue_tilt_azimuth_angle_LCS_size; a++)
				{
					alpha = 0;
					alpha_panel_2 = 0;

					for (int m = 0; m < M; m++)
					{
						for (int n = 0; n < N; n++)
						{
							for (int polar = 0; polar < P; polar++)
							{
								alpha_zero = 0;
								alpha_nmup = 0;
								alpha_zero_panel_2 = 0;
								alpha_nmup_panel_2 = 0;

								UE_antenna_gain = Get_MS_antgain(m, n, polar, channel_of_interest, 0, z, a, ue_interference);   ///// MAX antgain - combined gain;
								Get_RX_SmallScale_antgain(m, n, polar, channel_of_interest, 0, z, a, ue_interference);

								if (channel_of_interest->Propagation == LOS_propagation)
								{
									// panel 1
									alpha_zero_temp = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
									alpha_zero = abs(alpha_zero_temp)*abs(alpha_zero_temp);
									// SNS: per-element attenuation now in Get_BS_antenna_field_pattern

									// panel 2
									alpha_zero_temp_panel_2 = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
									alpha_zero_panel_2 = abs(alpha_zero_temp_panel_2)*abs(alpha_zero_temp_panel_2);
									// SNS: per-element attenuation now in Get_BS_antenna_field_pattern
								}

								for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
								{
									Real cluster_power_p1 = 0.;
									Real cluster_power_p2 = 0.;
									for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
									{
										kappa = channel_of_interest->kappa[i][j];
										_1_over_sqrt_K = 1.0 / (Real)sqrt(kappa);

										complex<Real> Big_pi_NLOS_thetatheta(0, channel_of_interest->random_phase_vv[i][j] * (pi / 180.0));   // use stored initial phase
										complex<Real> Big_pi_NLOS_thetapi(0, channel_of_interest->random_phase_vh[i][j] * (pi / 180.0));
										complex<Real> Big_pi_NLOS_pitheta(0, channel_of_interest->random_phase_hv[i][j] * (pi / 180.0));
										complex<Real> Big_pi_NLOS_pipi(0, channel_of_interest->random_phase_hh[i][j] * (pi / 180.0));

										// panel 1
										alpha_nmup_temp = (Real)sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
											(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
											(SmallScale_RX_AntennaGainXLOS_theta(i, j) * (Real)_1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
											);

										cluster_power_p1 += (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));

										/// panel 2
										alpha_nmup_temp_panel_2 = (Real)sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
											(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)* _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
											(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * (Real)_1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
											);

										cluster_power_p2 += (abs(alpha_nmup_temp_panel_2) * abs(alpha_nmup_temp_panel_2));

									}
									alpha_nmup += cluster_power_p1;  // SNS: per-element in Get_BS_antenna_field_pattern
									alpha_nmup_panel_2 += cluster_power_p2;  // SNS: per-element in Get_BS_antenna_field_pattern
								}
								alpha += alpha_zero + alpha_nmup;
								alpha_panel_2 += alpha_zero_panel_2 + alpha_nmup_panel_2;
							}
						}
					}

					channel_of_interest->signal_RSRP_gain[sec_number][sec_z_idx][sec_a_idx][z][a][0] = alpha / (Real)(M*N*P);          //[sec_idx] [sec_zenith][sec_azimuth][ue_zenith_idx][ue_azimuth_idx][pannel];
					channel_of_interest->signal_RSRP_gain[sec_number][sec_z_idx][sec_a_idx][z][a][1] = alpha_panel_2 / (Real)(M*N*P);  //[sec_idx] [sec_zenith][sec_azimuth][ue_zenith_idx][ue_azimuth_idx][pannel];

					if (a == 0 && z == 0)  // first loop
					{
						if (alpha > alpha_panel_2)
						{
							MAX_alpha = alpha;

							selected_a  = a;
							selected_z  = z;
							selected_p  = 0;
						}
						else
						{
							MAX_alpha = alpha_panel_2;

							selected_a  = a;
							selected_z  = z;
							selected_p  = 1;
						}
					}
					else if ((MAX_alpha < alpha) || (MAX_alpha < alpha_panel_2))
					{
						if (alpha > alpha_panel_2)
						{
							MAX_alpha = alpha;

							selected_a = a;
							selected_z = z;
							selected_p = 0;
						}
						else
						{
							MAX_alpha = alpha_panel_2;

							selected_a = a;
							selected_z = z;
							selected_p = 1;
						}
					}
					else
					{
						// 220810 jhnoh
						//MAX_alpha = MAX_alpha;

						//selected_a = selected_a;
						//selected_z = selected_z;
						//selected_p = selected_p;
					}
				}
			}

			Real RSRP_antenna_gain = (MAX_alpha) / (M*N*P);
			return 10 * log10(RSRP_antenna_gain);
		}
	}
}


Real Get_LCS_theta(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{
	Real x = cos(beta)*cos(gamma)*cos(GCS_theta) + (sin(beta)*cos(gamma)*cos(GCS_pi-alpha) - sin(gamma)*sin(GCS_pi-alpha))*sin(GCS_theta);
	// 부동소수점 오차로 |x| > 1 → acos(NaN) 방지
	if (x > 1.0) x = 1.0;
	if (x < -1.0) x = -1.0;
	return acosl(x);
}

Real Get_LCS_pi(Real alpha, Real beta, Real gamma, Real GCS_theta, Real GCS_pi)
{

	// arg(x+iy) = atan2(y,x)


	Real real;
	Real imag;

	real = (cos(beta)*sin(GCS_theta)*cos(GCS_pi - alpha) - sin(beta)*cos(GCS_theta));
	imag = (cos(beta)*sin(gamma)*cos(GCS_theta) + (sin(beta)*sin(gamma)*cos(GCS_pi - alpha) + cos(gamma)*sin(GCS_pi - alpha))*sin(GCS_theta));

	return atan2l(imag,real);

}


MatrixXReal Get_distance_angular(Real a, Real b, Real c, Real x, Real y, Real z)
{
	///// TR 36.873 page 9

	MatrixXReal Rd(3, 1);

	Rd(0, 0) = cos(a)*cos(b) * x + cos(a)*sin(b)*sin(c) - sin(a)*cos(c) * y + cos(a)*sin(b)*cos(c) + sin(a)*sin(c) * z;
	Rd(1, 0) = sin(a)*cos(b) * x + sin(a)*sin(b)*sin(c) + cos(a)*cos(c) * y + sin(a)*sin(b)*cos(c) - cos(a)*sin(c) * z;
	Rd(2, 0) = -sin(b)		 * x + cos(b)*sin(c)					    * y + cos(b)*cos(c)						   * z;

	return Rd;
}


Real LINK::Get_MS_antgain(int M, int N, int P, 
	                      CHANNEL * channel_of_interest , int sector_idx, 
						  int tilt_z_idx, int tilt_a_idx, int ue_interference_index)
{
	int ms_idx = self_ms_idx;

	Real h_angle_pi = channel_of_interest->LOS_AOA_GCS * (pi / 180.);  //// pi in GCS, rad
	Real v_angle_theta = channel_of_interest->LOS_ZOA_GCS * (pi / 180.);  ///// theta in GCS, rad

	Real F_theta_GCS_P1 = 0;
	Real F_pi_GCS_P1 = 0;
	Real F_theta_GCS_P2 = 0;
	Real F_pi_GCS_P2 = 0;

	ComplexReal  F_rx_theta   = {0, 0};
	ComplexReal  F_rx_pi      = {0, 0};
	ComplexReal  F_rx_theta_2 = {0, 0};
	ComplexReal  F_rx_pi_2    = {0, 0};

	Real combined_antenna_gain;

	if (ue_antenna_element_gain == 0 && !handheld_mode) // UE omni-antenna case
	{
		combined_antenna_gain = 0;
		if (P == 0)
		{
			ReceiverAntennaGainXLOS_theta = 1.;
			ReceiverAntennaGainXLOS_pi = 0.;
		}
		else
		{
			ReceiverAntennaGainXLOS_theta = 0.;
			ReceiverAntennaGainXLOS_pi = 1.;
		}
	}
	else
	{
		combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, handheld_mode ? M : -1);
		Get_UE_antenna_field_pattern(channel_of_interest->self_ms, M, N, P, ms_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_rx_theta, F_rx_pi, F_rx_theta_2, F_rx_pi_2);

		ReceiverAntennaGainXLOS_theta = F_rx_theta;
		ReceiverAntennaGainXLOS_pi = F_rx_pi;

		ReceiverAntennaGainXLOS_theta_panel_2 = F_rx_theta_2;
		ReceiverAntennaGainXLOS_pi_panel_2 = F_rx_pi_2;
	}


	return combined_antenna_gain;
}

Real Get_BS_antenna_pattern(Real theta_GCS, Real pi_GCS, int bs_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2)
{
	if (TYPE == 11)   ///InH
	{
		Real v_angle_theta_LCS_temp = 0.;
		Real h_angle_pi_LCS_temp = 0.;

		/////  GCS -> LCS with mechanical down tile angle 
		/////  TR 36.814 / A.2.1.6.2 / equation (1),(2)
		/////GCS -> LCS , use alpha, beta, gamma
		///// TR 36.873

		Real alpha = bs[bs_idx].ant[sector_index][0];
		Real beta  = bs[bs_idx].ant[sector_index][1];
		Real gamma = bs[bs_idx].ant[sector_index][2];

		// Polarization slant angle: cross-pol(P=2) ±45°, co-pol(P=1) 0°
		Real zeta_1, zeta_2;
		if (BS_P == 2) { zeta_1 = pi / 4; zeta_2 = -pi / 4; }
		else           { zeta_1 = 0;       zeta_2 = 90;        }

		Real F_theta_LCS_P1;
		Real F_pi_LCS_P1;
		Real F_theta_LCS_P2;
		Real F_pi_LCS_P2;

		Real GCS_field_theta_P1;
		Real GCS_field_pi_P1;
		Real GCS_field_theta_P2;
		Real GCS_field_pi_P2;

		v_angle_theta_LCS_temp = Get_LCS_theta(alpha, beta, gamma, theta_GCS, pi_GCS);
		h_angle_pi_LCS_temp = Get_LCS_pi(alpha, beta, gamma, theta_GCS, pi_GCS);

		v_angle_theta_LCS_temp = v_angle_theta_LCS_temp * (180. / pi);
		h_angle_pi_LCS_temp = h_angle_pi_LCS_temp * (180. / pi);

		Real h_antgain = -MIN(12. * (h_angle_pi_LCS_temp / 90.)*(h_angle_pi_LCS_temp / 90.), 25.);
		Real v_antgain = -MIN(12. * ((v_angle_theta_LCS_temp - 90.) / 90.) * ((v_angle_theta_LCS_temp - 90.) / 90.), 25.);

		Real combined_antenna_gain = max_antgain - MIN(-(h_antgain + v_antgain), 25.);

		/// Polarized antennas Model-2
		/// 36.873 -> 7.1.1
		Real A = sqrt(pow(10.0, (combined_antenna_gain / 10.0)));

		F_theta_LCS_P1 = A * cos(zeta_1);
		F_pi_LCS_P1    = A * sin(zeta_1);
		F_theta_LCS_P2 = A * cos(zeta_2);
		F_pi_LCS_P2    = A * sin(zeta_2);
		//// 

		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1, F_pi_LCS_P1, GCS_field_theta_P1, GCS_field_pi_P1);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2, F_pi_LCS_P2, GCS_field_theta_P2, GCS_field_pi_P2);

		F_theta_GCS_P1 = GCS_field_theta_P1;
		F_pi_GCS_P1    = GCS_field_pi_P1;
		F_theta_GCS_P2 = GCS_field_theta_P2;
		F_pi_GCS_P2    = GCS_field_pi_P2;

		return combined_antenna_gain;
	}
	else if (TYPE == 12) /// Dense Urban
	{
		int _sector_idx = sector_index;

		Real v_angle_theta_LCS_temp = 0.;
		Real h_angle_pi_LCS_temp = 0.;

		Real alpha = bs[bs_idx].ant[_sector_idx][0];
		Real beta = bs[bs_idx].ant[_sector_idx][1];
		Real gamma = bs[bs_idx].ant[_sector_idx][2];

		// Polarization slant angle: cross-pol(P=2) ±45°, co-pol(P=1) 0°
		Real zeta_1, zeta_2;
		if (BS_P == 2) { zeta_1 = pi / 4; zeta_2 = -pi / 4; }
		else           { zeta_1 = 0;       zeta_2 = 0;        }

		Real F_theta_LCS_P1;
		Real F_pi_LCS_P1;
		Real F_theta_LCS_P2;
		Real F_pi_LCS_P2;

		Real GCS_field_theta_P1;
		Real GCS_field_pi_P1;
		Real GCS_field_theta_P2;
		Real GCS_field_pi_P2;

		v_angle_theta_LCS_temp = Get_LCS_theta(alpha, beta, gamma, theta_GCS, pi_GCS);    /// radian
		h_angle_pi_LCS_temp = Get_LCS_pi(alpha, beta, gamma, theta_GCS, pi_GCS);

		v_angle_theta_LCS_temp = v_angle_theta_LCS_temp * (180. / pi);                    /// degree
		h_angle_pi_LCS_temp = h_angle_pi_LCS_temp * (180. / pi);

		Real h_antgain = -MIN(12. * (h_angle_pi_LCS_temp / 65.)*(h_angle_pi_LCS_temp / 65.), 30.);
		Real v_antgain = -MIN(12. * ((v_angle_theta_LCS_temp - 90.) / 65.) * ((v_angle_theta_LCS_temp - 90.) / 65.), 30.);

		Real combined_antenna_gain = max_antgain - MIN(-(h_antgain + v_antgain), 30.);

		/// Polarized antennas Model-2
		/// 36.873 -> 7.1.1
		Real A = sqrt(pow(10.0, (combined_antenna_gain / 10.0)));

		F_theta_LCS_P1 = A * cos(zeta_1);
		F_pi_LCS_P1    = A * sin(zeta_1);
		F_theta_LCS_P2 = A * cos(zeta_2);
		F_pi_LCS_P2    = A * sin(zeta_2);
		//// 

		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1, F_pi_LCS_P1, GCS_field_theta_P1, GCS_field_pi_P1);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2, F_pi_LCS_P2, GCS_field_theta_P2, GCS_field_pi_P2);

		F_theta_GCS_P1 = GCS_field_theta_P1;
		F_pi_GCS_P1    = GCS_field_pi_P1;
		F_theta_GCS_P2 = GCS_field_theta_P2;
		F_pi_GCS_P2    = GCS_field_pi_P2;

		return combined_antenna_gain;

	}
	else if (TYPE == 13) /// Rural
	{
		int _sector_idx = sector_index;

		Real v_angle_theta_LCS_temp = 0.;
		Real h_angle_pi_LCS_temp = 0.;

		Real alpha = bs[bs_idx].ant[_sector_idx][0];
		Real beta = bs[bs_idx].ant[_sector_idx][1];
		Real gamma = bs[bs_idx].ant[_sector_idx][2];

		// Polarization slant angle: cross-pol(P=2) ±45°, co-pol(P=1) 0°
		Real zeta_1, zeta_2;
		if (BS_P == 2) { zeta_1 = pi / 4; zeta_2 = -pi / 4; }
		else           { zeta_1 = 0;       zeta_2 = 0;        }

		Real F_theta_LCS_P1;
		Real F_pi_LCS_P1;
		Real F_theta_LCS_P2;
		Real F_pi_LCS_P2;

		Real GCS_field_theta_P1;
		Real GCS_field_pi_P1;
		Real GCS_field_theta_P2;
		Real GCS_field_pi_P2;

		v_angle_theta_LCS_temp = Get_LCS_theta(alpha, beta, gamma, theta_GCS, pi_GCS);    /// radian
		h_angle_pi_LCS_temp = Get_LCS_pi(alpha, beta, gamma, theta_GCS, pi_GCS);

		v_angle_theta_LCS_temp = v_angle_theta_LCS_temp * (180. / pi);                    /// degree
		h_angle_pi_LCS_temp = h_angle_pi_LCS_temp * (180. / pi);

		Real h_antgain = -MIN(12. * (h_angle_pi_LCS_temp / 65.)*(h_angle_pi_LCS_temp / 65.), 30.);
		Real v_antgain = -MIN(12. * ((v_angle_theta_LCS_temp - 90.) / 65.) * ((v_angle_theta_LCS_temp - 90.) / 65.), 30.);

		Real combined_antenna_gain = max_antgain - MIN(-(h_antgain + v_antgain), 30.);


		/// Polarized antennas Model-2
		/// 36.873 -> 7.1.1
		Real A = sqrt(pow(10.0, (combined_antenna_gain / 10.0)));

		F_theta_LCS_P1 = A * cos(zeta_1);
		F_pi_LCS_P1    = A * sin(zeta_1);
		F_theta_LCS_P2 = A * cos(zeta_2);
		F_pi_LCS_P2    = A * sin(zeta_2);
		//// 

		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1, F_pi_LCS_P1, GCS_field_theta_P1, GCS_field_pi_P1);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2, F_pi_LCS_P2, GCS_field_theta_P2, GCS_field_pi_P2);

		F_theta_GCS_P1 = GCS_field_theta_P1;
		F_pi_GCS_P1    = GCS_field_pi_P1;
		F_theta_GCS_P2 = GCS_field_theta_P2;
		F_pi_GCS_P2    = GCS_field_pi_P2;

		return combined_antenna_gain;
	}
	else {
		cout << "Something wrong in Get_BS_antenna_pattern" << endl;
		getchar();
		return -1;
	}
}


ComplexReal  Get_BS_antenna_field_pattern(LOCATION interferer,
	int _bs_idx,
	int _sec_idx,
	int tilt_z_idx,
	int tilt_a_idx,
	Real v_angle_theta,
	Real h_angle_pi,
	Real F_theta_GCS_P1,
	Real F_pi_GCS_P1,
	Real F_theta_GCS_P2,
	Real F_pi_GCS_P2,
	ComplexReal  &F_tx_theta,
	ComplexReal  &F_tx_pi,
	const ClusterVR* sns_vr_ptr,
	Real sns_rolloff_C)
{
	Real dH;
	Real dV;

	dH = BS_dH;
	dV = BS_dV;

	int K = BS_M / BS_Mp; //   4/4 = 1
	int L = BS_N / BS_Np; //   4/4 = 1

	ComplexReal jay(0.0, 1.0);
	ComplexReal  F_theta_temp = {0,0};
	ComplexReal  F_pi_temp    = {0,0};
	ComplexReal  weight = 0;
	ComplexReal  w;
	LOCATION3D d_tx;
	LOCATION3D r_tx;

	Real Theta = v_angle_theta;
	Real Phi = h_angle_pi;

	r_tx.x = sin(Theta) * cos(Phi);
	r_tx.y = sin(Theta) * sin(Phi);
	r_tx.z = cos(Theta);

	weight = {0, 0};
	ComplexReal weight_optimized = {0, 0};  // For comparison

	// Precompute constant factor outside loop
	const Real phase_factor = REAL(2.0) * pi / Wavelength;

	// Set to 1 to enable comparison (WARNING: slower due to running both versions!)
	// Set to 0 to use default method (see USE_OPTIMIZED_PHASE below)
	#define COMPARE_PHASE_CALCULATION 0

	// Set to 1 to use optimized cos/sin method, 0 to use original exp() method
	// Only applies when COMPARE_PHASE_CALCULATION is 0
	#define USE_OPTIMIZED_PHASE 1

	#if COMPARE_PHASE_CALCULATION
	auto start_original = std::chrono::high_resolution_clock::now();
	#endif

	for (int k = 0; k < K; k++) // vertical element
	{
		for (int l = 0; l < L; l++) // horizental element
		{
			w = virtualization_weight_wv[tilt_z_idx][tilt_a_idx][k][l];
			d_tx.x = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].x;// +interferer.x; //// sector(rx) [M][N][P][Mg][Ng] <port0>   channel[bs_idx][ms_idx].self_bs.x
			d_tx.y = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].y;// +interferer.y; //// sector(rx) [M][N][P][Mg][Ng] <port0>   channel[bs_idx][ms_idx].self_bs.y
			d_tx.z = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].z;// +bs_height;    //// sector(rx) [M][N][P][Mg][Ng] <port0>   bs_height

			// SNS: per-element per-cluster attenuation (when VR info provided)
			Real sns_w = REAL(1.0);
			if (sns_vr_ptr != nullptr) {
				Real pos_h = l * dH;
				Real pos_v = k * dV;
				sns_w = compute_sns_attenuation(pos_h, pos_v, *sns_vr_ptr, sns_rolloff_C);
			}

			#if COMPARE_PHASE_CALCULATION
			// Original version: using exp()
			weight += w * sns_w * exp(jay * phase_factor * REAL(dot(r_tx, d_tx)));
			#elif USE_OPTIMIZED_PHASE
			// Optimized version: Use Euler's formula directly
			// exp(i*theta) = cos(theta) + i*sin(theta)
			Real theta = phase_factor * REAL(dot(r_tx, d_tx));
			Real cos_theta = cos(theta);
			Real sin_theta = sin(theta);
			weight += w * sns_w * ComplexReal(cos_theta, sin_theta);
			#else
			// Original version: using exp()
			weight += w * sns_w * exp(jay * phase_factor * REAL(dot(r_tx, d_tx)));
			#endif
		}
	}

	#if COMPARE_PHASE_CALCULATION
	auto end_original = std::chrono::high_resolution_clock::now();

	// Now test optimized version
	auto start_optimized = std::chrono::high_resolution_clock::now();
	for (int k = 0; k < K; k++) // vertical element
	{
		for (int l = 0; l < L; l++) // horizental element
		{
			w = virtualization_weight_wv[tilt_z_idx][tilt_a_idx][k][l];
			d_tx.x = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].x;
			d_tx.y = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].y;
			d_tx.z = bs[_bs_idx].d_tx[_sec_idx][k][l][0][0][0].z;

			// SNS: per-element per-cluster attenuation (when VR info provided)
			Real sns_w2 = REAL(1.0);
			if (sns_vr_ptr != nullptr) {
				Real pos_h = l * dH;
				Real pos_v = k * dV;
				sns_w2 = compute_sns_attenuation(pos_h, pos_v, *sns_vr_ptr, sns_rolloff_C);
			}

			// Optimized version: Use Euler's formula directly
			// exp(i*theta) = cos(theta) + i*sin(theta)
			Real theta = phase_factor * REAL(dot(r_tx, d_tx));
			Real cos_theta = cos(theta);
			Real sin_theta = sin(theta);
			weight_optimized += w * sns_w2 * ComplexReal(cos_theta, sin_theta);
		}
	}
	auto end_optimized = std::chrono::high_resolution_clock::now();

	// Calculate elapsed time in microseconds
	auto duration_original = std::chrono::duration_cast<std::chrono::microseconds>(end_original - start_original).count();
	auto duration_optimized = std::chrono::duration_cast<std::chrono::microseconds>(end_optimized - start_optimized).count();

	// Compare results (only for first few calls to avoid spam)
	static int call_count = 0;
	if (call_count < 3 && (_bs_idx == 0 && _sec_idx == 0))
	{
		cout << "========== Phase Calculation Comparison (Call #" << call_count + 1 << ") ==========" << endl;
		cout << "Original  weight: " << weight.real() << " + " << weight.imag() << "i" << endl;
		cout << "Optimized weight: " << weight_optimized.real() << " + " << weight_optimized.imag() << "i" << endl;
		cout << "Difference (real): " << abs(weight.real() - weight_optimized.real()) << endl;
		cout << "Difference (imag): " << abs(weight.imag() - weight_optimized.imag()) << endl;
		cout << "Relative error: " << abs(weight - weight_optimized) / abs(weight) * 100.0 << "%" << endl;
		cout << endl;
		cout << "Execution time comparison:" << endl;
		cout << "  Original (exp):  " << duration_original << " μs" << endl;
		cout << "  Optimized (cos/sin): " << duration_optimized << " μs" << endl;
		cout << "  Speedup: " << (Real)duration_original / duration_optimized << "x" << endl;
		cout << "================================================================" << endl;
		if (call_count == 0) {
			cout << "Press Enter to continue..." << endl;
			getchar();
		}
		call_count++;
	}
	#endif

	F_theta_temp = (F_theta_GCS_P1  * weight);
	F_pi_temp    = (F_pi_GCS_P1     * weight);


	/// F, antenna gain for RSRP
	F_tx_theta = F_theta_temp;
	F_tx_pi    = F_pi_temp;

	return 0;
}


Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2, int port_idx)
{
	// ================================================================
	// Handheld UT antenna model (TR 38.901 Section 7.8)
	// 2-stage rotation: Reference pattern → UT LCS → GCS
	// ================================================================
	if (handheld_mode && port_idx >= 0)
	{
		// Per-antenna Euler angles (α_u): boresight direction from center to antenna position
		// For 8 candidate positions on 15cm×7cm device
		// α_u = atan2(pos.x, -pos.y) maps center→antenna direction to azimuth
		static const Real HANDHELD_ALPHA_U[8] = {
			atan2(-0.075, 0.035),   // pos 1: (-L/2, -W/2) → ≈ -65°
			atan2( 0.0,   0.035),   // pos 2: (0, -W/2)    → 0°
			atan2( 0.075, 0.035),   // pos 3: (+L/2, -W/2) → ≈ +65°
			atan2( 0.075, 0.0  ),   // pos 4: (+L/2, 0)    → 90°
			atan2( 0.075,-0.035),   // pos 5: (+L/2, +W/2) → ≈ +115°
			atan2( 0.0,  -0.035),   // pos 6: (0, +W/2)    → 180°
			atan2(-0.075,-0.035),   // pos 7: (-L/2, +W/2) → ≈ -115°
			atan2(-0.075, 0.0  ),   // pos 8: (-L/2, 0)    → -90°
		};
		static const Real HANDHELD_BETA_U  = pi / 2.0;   // 90° for all antennas
		static const Real HANDHELD_GAMMA_U = 0.0;         // 0° for all antennas

		int ant_idx = handheld_port_indices[port_idx] - 1;  // 1-based → 0-based
		Real alpha_u = HANDHELD_ALPHA_U[ant_idx];
		Real beta_u  = HANDHELD_BETA_U;
		Real gamma_u = HANDHELD_GAMMA_U;

		// UE device orientation (Ω_α, Ω_β, Ω_γ)
		Real omega_alpha = ms[ms_idx].alpha;
		Real omega_beta  = ms[ms_idx].beta;
		Real omega_gamma = ms[ms_idx].gamma;

		// Stage 2 inverse: GCS → UT LCS
		Real theta_p = Get_LCS_theta(omega_alpha, omega_beta, omega_gamma, theta_GCS, pi_GCS);
		Real phi_p   = Get_LCS_pi(omega_alpha, omega_beta, omega_gamma, theta_GCS, pi_GCS);

		// Stage 1 inverse: UT LCS → Reference pattern frame
		Real theta_pp = Get_LCS_theta(alpha_u, beta_u, gamma_u, theta_p, phi_p);
		Real phi_pp   = Get_LCS_pi(alpha_u, beta_u, gamma_u, theta_p, phi_p);

		// Element pattern: isotropic (0 dBi) or Table 7.3-2 directive (5.3 dBi max)
		Real A;
		Real combined_gain;
		if (ue_antenna_element_gain == 0) {
			// Isotropic element: 0 dBi amplitude
			// Per-antenna polarization direction still applied via (α_u, β_u, γ_u) rotation
			A = 1.0;
			combined_gain = 0.0;
		} else {
			// Table 7.3-2 element pattern (θ_3dB=φ_3dB=125°, A_max=22.5 dB, G_max=5.3 dBi)
			Real theta_pp_deg = theta_pp * (180.0 / pi);
			Real phi_pp_deg   = phi_pp * (180.0 / pi);
			Real v_gain = -MIN(12.0 * ((theta_pp_deg - 90.0) / 125.0) * ((theta_pp_deg - 90.0) / 125.0), 22.5);
			Real h_gain = -MIN(12.0 * (phi_pp_deg / 125.0) * (phi_pp_deg / 125.0), 22.5);
			combined_gain = 5.3 - MIN(-(h_gain + v_gain), 22.5);  // dBi
			A = sqrt(pow(10.0, combined_gain / 10.0));
		}

		// Single-pol: ζ=0 → F'_θ = A, F'_φ = 0
		// Each antenna's unique polarization direction comes from α_u rotation (Section 4.5)
		Real F_theta_ref = A;
		Real F_phi_ref   = 0.0;

		// Stage 1 field rotation: reference frame → UT LCS
		Real F_theta_ut, F_phi_ut;
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha_u, beta_u, gamma_u,
			theta_p, phi_p, F_theta_ref, F_phi_ref, F_theta_ut, F_phi_ut);

		// Stage 2 field rotation: UT LCS → GCS
		LCS_Antenna_field_to_GCS_antenna_pattern(omega_alpha, omega_beta, omega_gamma,
			theta_GCS, pi_GCS, F_theta_ut, F_phi_ut, F_theta_GCS_P1, F_pi_GCS_P1);

		// Handheld = single polarization per port, no panel 2
		F_theta_GCS_P2 = 0.0;
		F_pi_GCS_P2    = 0.0;

		return combined_gain;
	}

	if (TYPE == 11 || TYPE == 12)   ///InH
	{
		Real v_angle_theta_LCS_temp = 0.;
		Real h_angle_pi_LCS_temp = 0.;
		Real v_angle_theta_LCS_temp_Panel2 = 0.;
		Real h_angle_pi_LCS_temp_Panel2 = 0.;

		/////  GCS -> LCS with mechanical down tile angle 
		/////  TR 36.814 / A.2.1.6.2 / equation (1),(2)
		Real alpha = ms[ms_idx].alpha;
		Real beta = ms[ms_idx].beta;
		Real gamma = ms[ms_idx].gamma;
		Real alpha_plus_pi = ms[ms_idx].alpha + pi;

		// Polarization slant angle: 0°/90° (radian)
		Real zeta_1, zeta_2;
		if (MS_P == 2) { zeta_1 = 0; zeta_2 = pi / 2; }
		else           { zeta_1 = 0; zeta_2 = pi / 2;  }

		Real F_theta_LCS_P1;
		Real F_pi_LCS_P1;
		Real F_theta_LCS_P2;
		Real F_pi_LCS_P2;

		Real F_theta_LCS_P1_Panel2;
		Real F_pi_LCS_P1_Panel2;
		Real F_theta_LCS_P2_Panel2;
		Real F_pi_LCS_P2_Panel2;

		Real GCS_field_theta_P1;
		Real GCS_field_pi_P1;
		Real GCS_field_theta_P2;
		Real GCS_field_pi_P2;

		Real GCS_field_theta_P1_Panel2;
		Real GCS_field_pi_P1_Panel2;
		Real GCS_field_theta_P2_Panel2;
		Real GCS_field_pi_P2_Panel2;

		v_angle_theta_LCS_temp = Get_LCS_theta(alpha, beta, gamma, theta_GCS, pi_GCS);
		h_angle_pi_LCS_temp    = Get_LCS_pi(alpha, beta, gamma, theta_GCS, pi_GCS);
		v_angle_theta_LCS_temp_Panel2 = Get_LCS_theta(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS);
		h_angle_pi_LCS_temp_Panel2    = Get_LCS_pi(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS);

		v_angle_theta_LCS_temp = v_angle_theta_LCS_temp * (180. / pi);   // degree
		h_angle_pi_LCS_temp    = h_angle_pi_LCS_temp    * (180. / pi);
		v_angle_theta_LCS_temp_Panel2 = v_angle_theta_LCS_temp_Panel2 * (180. / pi);   // degree
		h_angle_pi_LCS_temp_Panel2    = h_angle_pi_LCS_temp_Panel2    * (180. / pi);

		Real h_antgain        = -MIN(12. * (h_angle_pi_LCS_temp / 90.)*(h_angle_pi_LCS_temp / 90.), 25.);
		Real v_antgain        = -MIN(12. * ((v_angle_theta_LCS_temp - 90.) / 90.) * ((v_angle_theta_LCS_temp - 90.) / 90.), 25.);
		Real h_antgain_Panel2 = -MIN(12. * (h_angle_pi_LCS_temp_Panel2 / 90.)*(h_angle_pi_LCS_temp_Panel2 / 90.), 25.);
		Real v_antgain_Panel2 = -MIN(12. * ((v_angle_theta_LCS_temp_Panel2 - 90.) / 90.) * ((v_angle_theta_LCS_temp_Panel2 - 90.) / 90.), 25.);

		Real combined_antenna_gain        = ue_antenna_element_gain - MIN(-(h_antgain + v_antgain), 25.);
		Real combined_antenna_gain_Panel2 = ue_antenna_element_gain - MIN(-(h_antgain_Panel2 + v_antgain_Panel2), 25.);

		if (carrier_freq <= 6000000000)  //below 6Ghz
		{
			combined_antenna_gain = 0.;
			combined_antenna_gain_Panel2 = 0.;
		}

		Real A  = sqrt(pow(10.0, (combined_antenna_gain / 10.0)));
		Real A2 = sqrt(pow(10.0, (combined_antenna_gain_Panel2 / 10.0)));

		F_theta_LCS_P1 = A * cos(zeta_1);
		F_pi_LCS_P1    = A * sin(zeta_1);
		F_theta_LCS_P2 = A * cos(zeta_2);
		F_pi_LCS_P2    = A * sin(zeta_2);

		F_theta_LCS_P1_Panel2 = A2 * cos(zeta_1);
		F_pi_LCS_P1_Panel2    = A2 * sin(zeta_1);
		F_theta_LCS_P2_Panel2 = A2 * cos(zeta_2);
		F_pi_LCS_P2_Panel2    = A2 * sin(zeta_2);
		////

		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1, F_pi_LCS_P1, GCS_field_theta_P1, GCS_field_pi_P1);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2, F_pi_LCS_P2, GCS_field_theta_P2, GCS_field_pi_P2);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1_Panel2, F_pi_LCS_P1_Panel2, GCS_field_theta_P1_Panel2, GCS_field_pi_P1_Panel2);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2_Panel2, F_pi_LCS_P2_Panel2, GCS_field_theta_P2_Panel2, GCS_field_pi_P2_Panel2);
		
		if (ue_antenna_element_gain == 0)
		{
			F_theta_GCS_P1 = GCS_field_theta_P1;
			F_pi_GCS_P1    = GCS_field_pi_P1;
			F_theta_GCS_P2 = GCS_field_theta_P2;  
			F_pi_GCS_P2    = GCS_field_pi_P2;
		}
		else
		{
			if (P == 0)
			{
				F_theta_GCS_P1 = GCS_field_theta_P1;
				F_pi_GCS_P1    = GCS_field_pi_P1;
				F_theta_GCS_P2 = GCS_field_theta_P1_Panel2;  //// P2�� Panel2�� �޾ƿ����� �ߵ�. ������ Polarization������. �ϳ��� Polarization ���� ��� �����ϱ⋚����...
				F_pi_GCS_P2    = GCS_field_pi_P1_Panel2;
			}
			else
			{
				F_theta_GCS_P1 = GCS_field_theta_P2;
				F_pi_GCS_P1    = GCS_field_pi_P2;
				F_theta_GCS_P2 = GCS_field_theta_P2_Panel2;  //// P2�� Panel2�� �޾ƿ����� �ߵ�. ������ Polarization������. �ϳ��� Polarization ���� ��� �����ϱ⋚����...
				F_pi_GCS_P2    = GCS_field_pi_P2_Panel2;
			}

		}
		return combined_antenna_gain;
	}
	else if (TYPE == 13) /// Rural
	{
		Real v_angle_theta_LCS_temp = 0.;
		Real h_angle_pi_LCS_temp = 0.;
		Real v_angle_theta_LCS_temp_Panel2 = 0.;
		Real h_angle_pi_LCS_temp_Panel2 = 0.;

		/////  GCS -> LCS with mechanical down tile angle 
		/////  TR 36.814 / A.2.1.6.2 / equation (1),(2)
		Real alpha = ms[ms_idx].alpha;
		Real beta = ms[ms_idx].beta;
		Real gamma = ms[ms_idx].gamma;
		Real alpha_plus_pi = ms[ms_idx].alpha + pi;

		// Polarization slant angle: 0°/90°
		Real zeta_1, zeta_2;
		if (MS_P == 2) { zeta_1 = 0; zeta_2 = pi / 2; }
		else           { zeta_1 = 0; zeta_2 = 0;        }

		Real F_theta_LCS_P1;
		Real F_pi_LCS_P1;
		Real F_theta_LCS_P2;
		Real F_pi_LCS_P2;

		Real F_theta_LCS_P1_Panel2;
		Real F_pi_LCS_P1_Panel2;
		Real F_theta_LCS_P2_Panel2;
		Real F_pi_LCS_P2_Panel2;

		Real GCS_field_theta_P1;
		Real GCS_field_pi_P1;
		Real GCS_field_theta_P2;
		Real GCS_field_pi_P2;

		Real GCS_field_theta_P1_Panel2;
		Real GCS_field_pi_P1_Panel2;
		Real GCS_field_theta_P2_Panel2;
		Real GCS_field_pi_P2_Panel2;

		v_angle_theta_LCS_temp = Get_LCS_theta(alpha, beta, gamma, theta_GCS, pi_GCS);
		h_angle_pi_LCS_temp = Get_LCS_pi(alpha, beta, gamma, theta_GCS, pi_GCS);
		v_angle_theta_LCS_temp_Panel2 = Get_LCS_theta(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS);
		h_angle_pi_LCS_temp_Panel2 = Get_LCS_pi(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS);

		v_angle_theta_LCS_temp = v_angle_theta_LCS_temp * (180. / pi);   // degree
		h_angle_pi_LCS_temp = h_angle_pi_LCS_temp    * (180. / pi);
		v_angle_theta_LCS_temp_Panel2 = v_angle_theta_LCS_temp_Panel2 * (180. / pi);   // degree
		h_angle_pi_LCS_temp_Panel2 = h_angle_pi_LCS_temp_Panel2    * (180. / pi);

		Real h_antgain = -MIN(12. * (h_angle_pi_LCS_temp / 90.)*(h_angle_pi_LCS_temp / 90.), 25.);
		Real v_antgain = -MIN(12. * ((v_angle_theta_LCS_temp - 90.) / 90.) * ((v_angle_theta_LCS_temp - 90.) / 90.), 25.);
		Real h_antgain_Panel2 = -MIN(12. * (h_angle_pi_LCS_temp_Panel2 / 90.)*(h_angle_pi_LCS_temp_Panel2 / 90.), 25.);
		Real v_antgain_Panel2 = -MIN(12. * ((v_angle_theta_LCS_temp_Panel2 - 90.) / 90.) * ((v_angle_theta_LCS_temp_Panel2 - 90.) / 90.), 25.);

		Real combined_antenna_gain = ue_antenna_element_gain - MIN(-(h_antgain + v_antgain), 25.);
		Real combined_antenna_gain_Panel2 = ue_antenna_element_gain - MIN(-(h_antgain_Panel2 + v_antgain_Panel2), 25.);

		if (carrier_freq < 6000000000)  //below 6Ghz
		{
			combined_antenna_gain = 0.;
			combined_antenna_gain_Panel2 = 0.;
		}

		Real A = sqrt(pow(10.0, (combined_antenna_gain / 10.0)));
		Real A2 = sqrt(pow(10.0, (combined_antenna_gain_Panel2 / 10.0)));

		F_theta_LCS_P1 = A * cos(zeta_1);
		F_pi_LCS_P1 = A * sin(zeta_1);
		F_theta_LCS_P2 = A * cos(zeta_2);
		F_pi_LCS_P2 = A * sin(zeta_2);

		F_theta_LCS_P1_Panel2 = A2 * cos(zeta_1);
		F_pi_LCS_P1_Panel2 = A2 * sin(zeta_1);
		F_theta_LCS_P2_Panel2 = A2 * cos(zeta_2);
		F_pi_LCS_P2_Panel2 = A2 * sin(zeta_2);
		//// 


		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1, F_pi_LCS_P1, GCS_field_theta_P1, GCS_field_pi_P1);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2, F_pi_LCS_P2, GCS_field_theta_P2, GCS_field_pi_P2);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P1_Panel2, F_pi_LCS_P1_Panel2, GCS_field_theta_P1_Panel2, GCS_field_pi_P1_Panel2);
		LCS_Antenna_field_to_GCS_antenna_pattern(alpha_plus_pi, beta, gamma, theta_GCS, pi_GCS, F_theta_LCS_P2_Panel2, F_pi_LCS_P2_Panel2, GCS_field_theta_P2_Panel2, GCS_field_pi_P2_Panel2);

		if (ue_antenna_element_gain == 0)
		{
			F_theta_GCS_P1 = GCS_field_theta_P1;
			F_pi_GCS_P1 = GCS_field_pi_P1;
			F_theta_GCS_P2 = GCS_field_theta_P2;
			F_pi_GCS_P2 = GCS_field_pi_P2;
		}
		else
		{
			if (P == 0)
			{
				F_theta_GCS_P1 = GCS_field_theta_P1;
				F_pi_GCS_P1 = GCS_field_pi_P1;
				F_theta_GCS_P2 = GCS_field_theta_P1_Panel2;  //// P2 Panel2
				F_pi_GCS_P2 = GCS_field_pi_P1_Panel2;
			}
			else
			{
				F_theta_GCS_P1 = GCS_field_theta_P2;
				F_pi_GCS_P1 = GCS_field_pi_P2;
				F_theta_GCS_P2 = GCS_field_theta_P2_Panel2;  //// P2�� Panel2
				F_pi_GCS_P2 = GCS_field_pi_P2_Panel2;
			}
		}
		return combined_antenna_gain;
	}
	else {
		cout << " Something wrong with Get_UE_antenna_pattern" << endl;
		getchar();
		return -1;
	}
}

ComplexReal  Get_UE_antenna_field_pattern(LOCATION UE, int M, int N, int P, int ms_idx, int sector_idx, int tilt_z_idx, int tilt_a_idx, Real v_angle_theta, Real h_angle_pi, Real F_theta_GCS_P1, Real F_pi_GCS_P1, Real F_theta_GCS_P2, Real F_pi_GCS_P2, ComplexReal  &F_rx_theta_panel_1, ComplexReal  &F_rx_pi_panel_1, ComplexReal  &F_rx_theta_panel_2, ComplexReal  &F_rx_pi_panel_2)
{
	Real dH;
	Real dV;

	dH = MS_dH;
	dV = MS_dV;

	int K = MS_M / MS_Mp;    /////   4 / 4 = 1
	int L = MS_N / MS_Np;   //////   4/4 = 1

	if (ue_antenna_element_gain == 0)
	{
		K = MS_M;
		L = MS_N;
	}

	ComplexReal jay(0, 1);

	ComplexReal  F_theta_temp   = {0, 0};
	ComplexReal  F_pi_temp      = {0, 0};
	ComplexReal  F_theta_temp_2 = {0, 0};
	ComplexReal  F_pi_temp_2    = {0, 0};


	ComplexReal  weight[2];
	weight[0] = {0, 0};
	weight[1] = {0, 0};

	ComplexReal  w;
	LOCATION3D d_rx[2];
	LOCATION3D r_rx;

	Real Theta = v_angle_theta;
	Real Phi = h_angle_pi;

	r_rx.x = sin(Theta) * cos(Phi);
	r_rx.y = sin(Theta) * sin(Phi);
	r_rx.z = cos(Theta);

	if (ue_antenna_element_gain == 0)
	{

		w = 1.;
		d_rx[0].x = ms[ms_idx].d_rx[M][N][0][0][0].x;// +UE.x; //// [M][N][P][Mg][Ng] <panel0>	
		d_rx[0].y = ms[ms_idx].d_rx[M][N][0][0][0].y;// +UE.y; //// [M][N][P][Mg][Ng] <panel0>	
		d_rx[0].z = ms[ms_idx].d_rx[M][N][0][0][0].z;// +ms[ms_idx].MS_HEIGHT_FINAL - bs_height; //// [M][N][P][Mg][Ng] <panel0>	




		weight[0] = FAST_EXP(jay*REAL(2.0)*pi / Wavelength * REAL(dot(r_rx, d_rx[0])));


		F_theta_temp = (F_theta_GCS_P1  * weight[0]);
		F_pi_temp    = (F_pi_GCS_P1     * weight[0]);
		F_theta_temp_2 = (F_theta_GCS_P2  * weight[0]);   // P2�� ���⼭�� Polarization
		F_pi_temp_2    = (F_pi_GCS_P2     * weight[0]);


		/// F, antenna gain for RSRP
		F_rx_theta_panel_1 = F_theta_temp;
		F_rx_pi_panel_1 = F_pi_temp;

		F_rx_theta_panel_2 = F_theta_temp_2;
		F_rx_pi_panel_2 = F_pi_temp_2;
	}
	else
	{
		for (int k = 0; k < K; k++) // vertical element
		{
			for (int l = 0; l < L; l++) // horizental element
			{
				w = ue_virtualization_weight_wv[tilt_z_idx][tilt_a_idx][k][l];
				d_rx[0].x = ms[ms_idx].d_rx[k][l][0][0][0].x; // +UE.x; //// [M][N][P][Mg][Ng] <panel0>
				d_rx[0].y = ms[ms_idx].d_rx[k][l][0][0][0].y; // +UE.y; //// [M][N][P][Mg][Ng] <panel0>
				d_rx[0].z = ms[ms_idx].d_rx[k][l][0][0][0].z; // +ms[ms_idx].MS_HEIGHT_FINAL - bs_height; //// [M][N][P][Mg][Ng] <panel0>

				d_rx[1].x = ms[ms_idx].d_rx[k][l][0][0][1].x; // +UE.x; //// [M][N][P][Mg][Ng] <panel1>
				d_rx[1].y = ms[ms_idx].d_rx[k][l][0][0][1].y; // +UE.y; //// [M][N][P][Mg][Ng] <panel1>
				d_rx[1].z = ms[ms_idx].d_rx[k][l][0][0][1].z; // +ms[ms_idx].MS_HEIGHT_FINAL - bs_height; //// [M][N][P][Mg][Ng] <panel1>


				weight[0] += w * FAST_EXP(jay*REAL(2.0)*pi / Wavelength * REAL(dot(r_rx, d_rx[0])));
				weight[1] += w * FAST_EXP(jay*REAL(2.0)*pi / Wavelength * REAL(dot(r_rx, d_rx[1])));
			}
		}
		F_theta_temp = (F_theta_GCS_P1  * weight[0]);
		F_pi_temp    = (F_pi_GCS_P1     * weight[0]);

		F_theta_temp_2 = (F_theta_GCS_P2  * weight[1]);   // P2�� Panel2
		F_pi_temp_2    = (F_pi_GCS_P2     * weight[1]);


		/// F, antenna gain for RSRP
		F_rx_theta_panel_1 = F_theta_temp;
		F_rx_pi_panel_1 = F_pi_temp;

		F_rx_theta_panel_2 = F_theta_temp_2;
		F_rx_pi_panel_2 = F_pi_temp_2;
	}


	return 0;
}



void LCS_Antenna_field_to_GCS_antenna_pattern(Real alpha, Real beta, Real gamma, Real GCS_angle_theta, Real GCS_angle_pi, Real LCS_field_theta, Real LCS_field_pi, Real &GCS_field_theta, Real &GCS_field_pi) //GCS_angle degree
{
	Real cos_psi;
	Real sin_psi;
	Real theta = GCS_angle_theta; // rad
	Real phi = GCS_angle_pi; // rad

	if ((pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2) == 1) || (pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2) == 1))
	{
		theta += 0.0001;
	}

	Real x = cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta);
	Real denom_sq = 1.0 - x * x;
	if (denom_sq < 1e-15) denom_sq = 1e-15;  // sqrt(음수) → NaN 방지
	Real denom = sqrt(denom_sq);

	cos_psi = (cos(beta)*cos(gamma)*sin(theta) - (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*cos(theta)) / denom;
	sin_psi = (sin(beta)*cos(gamma)*sin(phi - alpha) + sin(gamma)*cos(phi - alpha)) / denom;

	GCS_field_theta = cos_psi * LCS_field_theta - sin_psi * LCS_field_pi;
	GCS_field_pi = sin_psi * LCS_field_theta + cos_psi * LCS_field_pi;

	if (std::isnan(GCS_field_theta) || std::isnan(GCS_field_pi) ||
	    std::isinf(GCS_field_theta) || std::isinf(GCS_field_pi)) {
		static int nan_count_lcs = 0;
		if (nan_count_lcs < 10) {
			nan_count_lcs++;
			cout << "\n*** NaN/Inf in LCS_to_GCS ***" << endl;
			cout << "  theta=" << theta << " phi=" << phi
			     << " alpha=" << alpha << " beta=" << beta << " gamma=" << gamma << endl;
			cout << "  x=" << x << " x^2=" << x*x << " denom_sq=" << denom_sq << " denom=" << denom << endl;
			cout << "  cos_psi=" << cos_psi << " sin_psi=" << sin_psi << endl;
			cout << "  LCS_field: theta=" << LCS_field_theta << " pi=" << LCS_field_pi << endl;
			cout << "  GCS_field: theta=" << GCS_field_theta << " pi=" << GCS_field_pi << endl;
		}
	}
}


#if 0  // unused: PIot_sector_antenna_gain, PIot_ue_antenna_gain — never called
void PIot_sector_antenna_gain()
{
	ComplexReal  w;
	LOCATION3D r_rx;
	LOCATION3D d_rx[3];
	ComplexReal  weight[3];
	ComplexReal  j(0.0, 1.0);
	SPHERICAL_ANGLE A;
	SPHERICAL_ANGLE A_LCS[3];
	Real antenna_gain[3];
	ANTENNA_FIELD F[3];
	ANTENNA_FIELD F_LCS[3];

	ofstream x;
	ofstream y;

	int K = BS_M / BS_Mp;    /////   4 / 4 = 1
	int L = BS_N / BS_Np;   //////   4/4 = 1

	MatrixXcReal wei(K, L);
	Real dH;
	Real dV;

	dH = BS_dH;
	dV = BS_dV;

	x.open("output/x.dat");
	y.open("output/y.dat");
	for (int theta = 0; theta < 181; theta++)
	{
		for (int phi = -180; phi < 181; phi++)
		{

			r_rx = Transform_angle_to_spheical_vector(Real(phi), Real(theta));
			x << r_rx.x << '\t';
			y << r_rx.y << '\t';
		}
		x << endl;
		y << endl;
	}
	x.close();
	y.close();

	ofstream out[3];
	string filename[3];



	for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
	{
		cout << z << endl;
		for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
		{
			stringstream fn[3];
			fn[0] << "output/sector_gain/" << "sector_0_zenith_" << z << "_azimuth_" << a << ".dat";
			fn[1] << "output/sector_gain/" << "sector_1_zenith_" << z << "_azimuth_" << a << ".dat";
			fn[2] << "output/sector_gain/" << "sector_2_zenith_" << z << "_azimuth_" << a << ".dat";

			filename[0] = fn[0].str();
			filename[1] = fn[1].str();
			filename[2] = fn[2].str();
			fn[0].clear();
			fn[1].clear();
			fn[2].clear();

			out[0].open(filename[0]);
			out[1].open(filename[1]);
			out[2].open(filename[2]);

			for (int theta = 0; theta<181; theta++)
			{
				for (int phi = -180; phi<181; phi++)
				{
					A.phi = Real(phi) * (pi / 180.);;
					A.theta = Real(theta) * (pi / 180.);;

					for (int sector_idx = 0; sector_idx < 3; sector_idx++)
					{
						Real	F_theta_GCS_P1;
						Real F_pi_GCS_P1; 
						Real F_theta_GCS_P2;
						Real F_pi_GCS_P2;

						ComplexReal  F_tx_theta;
						ComplexReal  F_tx_pi;

						LOCATION AA;
						AA.x = 0;
						AA.y = 0;

						Real combined_antenna_gain = Get_BS_antenna_pattern(A.theta, A.phi, 0, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

						F[sector_idx].theta = F_theta_GCS_P1;
						F[sector_idx].phi = F_pi_GCS_P1;



						Get_BS_antenna_field_pattern(AA, 0, sector_idx, z, a, A.theta, A.phi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi);

						out[sector_idx] << pow(abs(F_tx_theta), 2) + pow(abs(F_tx_pi), 2) << '\t';

					}
				}
				out[0] << endl;
				out[1] << endl;
				out[2] << endl;
			}
			out[0].close();
			out[1].close();
			out[2].close();
		}
	}
	cout << "sector gain--" << endl;
}


void PIot_ue_antenna_gain()
{

	ComplexReal  w;
	LOCATION3D r_rx;
	LOCATION3D d_rx[3];
	ComplexReal  weight[3];
	ComplexReal  j(0.0, 1.0);
	SPHERICAL_ANGLE A;
	SPHERICAL_ANGLE A_LCS[3];
	Real antenna_gain[3];
	ANTENNA_FIELD F[3];
	ANTENNA_FIELD F_LCS[3];

	ofstream out0;
	ofstream out1;

	ofstream x;
	ofstream y;



	int K = MS_M / MS_Mp;    /////   4 / 4 = 1
	int L = MS_N / MS_Np;   //////   4/4 = 1

	MatrixXcReal wei(K, L);
	Real dH;
	Real dV;

	dH = MS_dH;
	dV = MS_dV;

	x.open("output/x.dat");
	y.open("output/y.dat");
	for (int theta = 0; theta < 181; theta++)
	{
		for (int phi = -180; phi < 181; phi++)
		{

			r_rx = Transform_angle_to_spheical_vector(Real(phi), Real(theta));
			x << r_rx.x << '\t';
			y << r_rx.y << '\t';
		}
		x << endl;
		y << endl;
	}
	x.close();
	y.close();


	string filename[2];



	for (int z = 0; z < ue_tilt_zenith_angle_LCS_size; z++)
	{
		cout << z << endl;
		for (int a = 0; a < ue_tilt_azimuth_angle_LCS_size; a++)
		{
			stringstream fn;
			stringstream fn1;

			fn << "output/" << "ue_pannel_0_zenith_" << z << "_azimuth_" << a << ".dat";
			fn1 << "output/" << "ue_pannel_1_zenith_" << z << "_azimuth_" << a << ".dat";
			filename[0] = fn.str();
			filename[1] = fn1.str();

			fn.clear();
			fn1.clear();

			out0.open(filename[0]);
			out1.open(filename[1]);

			for (int theta = 0; theta<181; theta++)
			{
				for (int phi = -180; phi<181; phi++)
				{
					A.phi = Real(phi) * (pi / 180.);;
					A.theta = Real(theta) * (pi / 180.);;


						Real	F_theta_GCS_P1;
						Real F_pi_GCS_P1;
						Real F_theta_GCS_P2;
						Real F_pi_GCS_P2;

						ComplexReal  F_tx_theta;
						ComplexReal  F_tx_pi;
						ComplexReal  F_tx_theta_Panel2;
						ComplexReal  F_tx_pi_Panel2;

						Real combined_antenna_gain = Get_UE_antenna_pattern(0, A.theta, A.phi, 0, 0, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

						F[0].theta = F_theta_GCS_P1;
						F[0].phi = F_pi_GCS_P1;

						LOCATION AA;
						AA.x = 0;
						AA.y = 0;

						Get_UE_antenna_field_pattern(AA, 0, 0, 0, 0, 0, z, a, A.theta, A.phi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi, F_tx_theta_Panel2, F_tx_pi_Panel2);

						out0 << pow(abs(F_tx_theta), 2) + pow(abs(F_tx_pi), 2) << '\t';
						out1 << pow(abs(F_tx_theta_Panel2), 2) + pow(abs(F_tx_pi_Panel2), 2) << '\t';


				}
				out0 << endl;
				out1 << endl;

			}
			out0.close();
			out1.close();
		}
	}

	cout << "ue gain--" << endl;
}
#endif  // unused: PIot_sector_antenna_gain, PIot_ue_antenna_gain

void LINK::Delete_link_memory()
{
	//delete[] LS_gain;
	//delete[] static_gain;

	SmallScale_TX_AntennaGainXLOS_theta.resize(0,0);
	SmallScale_TX_AntennaGainXLOS_pi.resize(0, 0);
	SmallScale_RX_AntennaGainXLOS_theta.resize(0, 0);
	SmallScale_RX_AntennaGainXLOS_pi.resize(0, 0);
	SmallScale_RX_AntennaGainXLOS_theta_panel_2.resize(0, 0);
	SmallScale_RX_AntennaGainXLOS_pi_panel_2.resize(0, 0);
}




#if 0  // unused: Antenna_field_in_LCS — never called
ANTENNA_FIELD Antenna_field_in_LCS(Real LCS_antenna_gain_in_dB, Real zeta)
{
	ANTENNA_FIELD field_pattern;
	Real LCS_antenna_gain = dB2linear(LCS_antenna_gain_in_dB);

	field_pattern.theta = pow(LCS_antenna_gain, 0.5)*cos(zeta);
	field_pattern.phi = pow(LCS_antenna_gain, 0.5)*sin(zeta);

	return field_pattern;
}
#endif  // unused: Antenna_field_in_LCS

Real dot(LOCATION3D a, LOCATION3D b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}



LOCATION3D Transform_angle_to_spheical_vector(Real phi, Real theta)
{
	LOCATION3D r;

	Real Theta = theta / 180 * pi;
	Real Phi = phi / 180 * pi;

	r.x = sin(Theta) * cos(Phi);
	r.y = sin(Theta) * sin(Phi);
	r.z = cos(Theta);

	return r;
}

#if 0  // unused: GCS_Angle_to_LCS_Angle, Sector_Antenna_gain_in_dB, UE_Antenna_gain_in_dB, LCS_Antenna_field_to_GCS_antenna_pattern_rev, Generate_vertical/horizontal_virtualization_weight_link
SPHERICAL_ANGLE GCS_Angle_to_LCS_Angle(Real alpha, Real beta, Real gamma, SPHERICAL_ANGLE GCS_angle) //GCS_angle : degree  ||  alpha, beta gamma : radian
{
	SPHERICAL_ANGLE LCS_angle;
	Real theta = GCS_angle.theta * pi / 180;
	Real phi = GCS_angle.phi * pi / 180;

	LCS_angle.theta = acosl(cos(beta) * cos(gamma) * cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta));

	LCS_angle.phi = atan2l((cos(beta)*sin(gamma)*cos(theta) + (sin(beta)*sin(gamma)*cos(phi - alpha) + cos(gamma)*sin(phi - alpha))*sin(theta)), (cos(beta)*sin(theta)*cos(phi - alpha) - sin(beta)*cos(theta)));

	LCS_angle.theta = LCS_angle.theta * 180 / pi;

	LCS_angle.phi = LCS_angle.phi * 180 / pi;
	return LCS_angle;
}


Real Sector_Antenna_gain_in_dB(SPHERICAL_ANGLE LCS_angle) //LCS_angle : degree
{
	Real _3dB_angle;
	Real SLA;
	Real Am;
	Real antenna_gain_Vertical;
	Real antenna_gain_Horizontal;
	Real antenna_gain;


	_3dB_angle = 65;
	SLA = 30;
	Am = 30;;



	antenna_gain_Vertical = - std::min(REAL(12) * pow(((LCS_angle.theta - REAL(90.0)) / _3dB_angle), REAL(2)), SLA);
	antenna_gain_Horizontal = - std::min(REAL(12) * pow((LCS_angle.phi / _3dB_angle), REAL(2)), Am);
	antenna_gain = -(std::min(-(antenna_gain_Vertical + antenna_gain_Horizontal), Am));

	antenna_gain = antenna_gain + max_antgain;

	return antenna_gain;

}

Real UE_Antenna_gain_in_dB(SPHERICAL_ANGLE LCS_angle) //LCS_angle : degree
{
	Real _3dB_angle;
	Real SLA;
	Real Am;
	Real antenna_gain_Vertical;
	Real antenna_gain_Horizontal;
	Real antenna_gain;

	_3dB_angle = 90;
	SLA = 25;
	Am = 25;;

	antenna_gain_Vertical = -min(REAL(12 * pow(((LCS_angle.theta - 90.0) / _3dB_angle), 2)), SLA);
	antenna_gain_Horizontal = -min(REAL(12 * pow((LCS_angle.phi / _3dB_angle), 2)), Am);
	antenna_gain = -(min(-(antenna_gain_Vertical + antenna_gain_Horizontal), Am));

	antenna_gain = antenna_gain + ue_antenna_element_gain;

	return antenna_gain;

}

ANTENNA_FIELD LCS_Antenna_field_to_GCS_antenna_pattern_rev(Real alpha, Real beta, Real gamma, SPHERICAL_ANGLE GCS_angle, ANTENNA_FIELD LCS_field) //GCS_angle degree
{
	ANTENNA_FIELD GCS_field;
	Real cos_psi;
	Real sin_psi;
	Real theta = GCS_angle.theta * pi / 180;
	Real phi = GCS_angle.phi * pi / 180;

	if ((pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2) == 1) || (pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2) == 1))
	{
		theta += 0.0001;
	}

	cos_psi = (cos(beta)*cos(gamma)*sin(theta) - (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*cos(theta)) / sqrt(1.0 - pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2));
	sin_psi = (sin(beta)*cos(gamma)*sin(phi - alpha) + sin(gamma)*cos(phi - alpha)) / sqrt(1.0 - pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2));;

	GCS_field.theta = cos_psi * LCS_field.theta - sin_psi * LCS_field.phi;
	GCS_field.phi = sin_psi * LCS_field.theta + cos_psi * LCS_field.phi;

	return GCS_field;

}



ComplexReal  Generate_vertical_virtualization_weight_link(Real theta, Real dv, Real k, Real K)
{
	ComplexReal  j(0.0, 1.0);
	ComplexReal  w;

	w = REAL(1.0) / REAL(sqrt(K))*FAST_EXP(-j * REAL(2.0) * pi / Wavelength * k * dv*REAL(cos(theta)));

	return w;
}


ComplexReal  Generate_horizontal_virtualization_weight_link(Real theta, Real phi, Real dH, Real l, Real L)
{
	ComplexReal  j(0.0, 1.0);
	ComplexReal  v;

	v = (Real)1.0 / sqrt(Real(L))*FAST_EXP(-j * ((Real)2.0 * pi / Wavelength * l * dH * (Real)sin(theta) * (Real)sin(phi)));

	return v;
}
#endif  // unused: GCS_Angle_to_LCS_Angle ... Generate_horizontal_virtualization_weight_link

Real Transform_angle_minus_180_to_plus_180(Real x)
{

	x = fmod(x, 360);

	if (x > 180)
	{
		x = x - 360;
	}
	else if (x < -180)
	{
		x = x + 360;
	}

	return x;
}

Real Transform_angle_0_to_plus_180(Real x)
{
	x = fmod(x, 360);
	if (x > 180)
	{
		x = 360 - x;
	}
	else if ((x < 0) && (x >= -180))
	{
		x = -x;
	}
	else if (x < -180)
	{
		x = x + 360;
	}
	return x;

}
// ============================================================================
// Profiling and Debugging Functions
// ============================================================================

#if PROFILE_FAST_EXP
void Print_FAST_EXP_Profile() {
    #if USE_FAST_COMPLEX_EXP
    std::cout << "\n========== FAST_EXP Profiling Results ==========" << std::endl;
    std::cout << "Optimized Version (cos/sin)" << std::endl;
    std::cout << "  Total calls: " << fast_exp_call_count << std::endl;
    std::cout << "  Total time: " << fast_exp_total_time << " μs" << std::endl;
    if (fast_exp_call_count > 0) {
        std::cout << "  Average time per call: " << fast_exp_total_time / fast_exp_call_count << " μs" << std::endl;
    }
    #else
    std::cout << "\n========== Standard exp() Profiling Results ==========" << std::endl;
    std::cout << "Original Version (exp())" << std::endl;
    std::cout << "  Total calls: " << exp_call_count << std::endl;
    std::cout << "  Total time: " << exp_total_time << " μs" << std::endl;
    if (exp_call_count > 0) {
        std::cout << "  Average time per call: " << exp_total_time / exp_call_count << " μs" << std::endl;
    }
    #endif
    std::cout << "================================================\n" << std::endl;
}
#endif


// =============================================================================
// Refactored signal/interference computation
// Decomposed from the original monolithic Get_signal_interference()
// into 11 smaller functions with clear responsibilities.
// =============================================================================

// Step 1: Build scenario configuration from global variables
ScenarioConfig LINK::build_scenario_config() const
{
	ScenarioConfig cfg;

	if (TYPE == 11 && num_Indoor_TRxP == 1)       // InH 1TRxP (12 TRxP)
	{
		cfg.num_candidate_bs                = num_BS + num_mTRP;
		cfg.sectors_per_bs                  = 1;
		cfg.use_wraparound                  = false;
		cfg.apply_penetration_loss          = false;
		cfg.use_single_cell_mode            = false;
		cfg.use_rsrp_table_for_interference = false;
		cfg.store_los_indoor_info           = false;
		cfg.configuration_type              = Configuration_Type;
		cfg.ms_to_bs_wrap_idx               = 0;
	}
	else if (TYPE == 11 && num_Indoor_TRxP == 3)   // InH 3TRxP (36 TRxP)
	{
		cfg.num_candidate_bs                = num_BS + num_mTRP;
		cfg.sectors_per_bs                  = 3;
		cfg.use_wraparound                  = false;
		cfg.apply_penetration_loss          = false;
		cfg.use_single_cell_mode            = false;
		cfg.use_rsrp_table_for_interference = false;
		cfg.store_los_indoor_info           = false;
		cfg.configuration_type              = Configuration_Type;
		cfg.ms_to_bs_wrap_idx               = 0;
	}
	else // Dense Urban (TYPE==12) or Rural (TYPE==13)
	{
		cfg.num_candidate_bs                = simple_num_BS;
		cfg.sectors_per_bs                  = 3;
		cfg.use_wraparound                  = true;
		cfg.apply_penetration_loss          = true;
		cfg.use_single_cell_mode            = (single_cell_mode == 1);
		cfg.use_rsrp_table_for_interference = true;
		cfg.store_los_indoor_info           = true;
		cfg.configuration_type              = Configuration_Type;
		cfg.ms_to_bs_wrap_idx               = self_ms_idx / (3 * num_MS_persector);
	}

	return cfg;
}


// Step 2: Compute final pathloss for a given BS
Real LINK::compute_pathloss_final(const ScenarioConfig& cfg, int bs_idx)
{
	CHANNEL* ch = &channel[bs_idx][self_ms_idx];
	Real pl = ch->pathloss;

	if (cfg.apply_penetration_loss)
		pl += incar_loss + Otoi_loss;

	pl += ms[self_ms_idx].LSPs[bs_idx](0);   // Add Shadow Fading
	ch->pathloss_final = pl;

	return pl;
}


// Step 3: Compute TX antenna gains (wrapper around Get_antgain)
Real LINK::compute_tx_antenna_gains(CHANNEL* ch, int bs_idx, int sector_idx, int tilt_a, int tilt_z)
{
	return Get_antgain(ch, bs_idx, sector_idx, tilt_a, tilt_z);
	// Side effect: sets TransmitterAntennaGainXLOS_theta/pi
}


// Step 4: Compute TX small-scale gains (unified, no TYPE branch needed)
void LINK::compute_tx_smallscale_gains(CHANNEL* ch, int bs_idx, int sector_idx, int tilt_a, int tilt_z)
{
	// The original Get_TX_SmallScale_antgain has TYPE==11 and TYPE==12/13 branches
	// but they are identical in logic. We simply delegate to the existing function.
	Get_TX_SmallScale_antgain(ch, bs_idx, sector_idx, tilt_a, tilt_z);
}


// Step 5: Compute RSRP (delegates to existing Get_RSRP)
Real LINK::compute_rsrp(CHANNEL* ch, int sec_number, int sec_z_idx, int sec_a_idx, int mode)
{
	// mode: 0 = find best UE beam & compute RSRP, 1 = fixed UE beam (for interference)
	return Get_RSRP(ch, sec_number, sec_z_idx, sec_a_idx, mode);
	// Side effect: sets selected_a, selected_z, selected_p (when mode==0)
	// Side effect: fills ch->signal_RSRP_gain table (when mode==0)
}


// Step 6: Find best TX beam for a given (BS, sector) pair
//
// Optimized RSRP computation following TR 36.873 Eq. (8.1-1)~(8.1-5).
//
// Key insight: The RSRP formula computes |α_{n,m,u,p}|² per ray (INCOHERENT sum),
// and there is NO RX spatial phase in the formula. Therefore:
//
//   |α_{n,m,u,p}|² = (P_n / (M_n · (K_R+1))) · |raysPreComp[u_p][n][m]|² · |AF(n,m,beam)|²
//
// where:
//   raysPreComp = F_rx,u × PolarizationMatrix × F_tx,element   (beam-INDEPENDENT)
//   AF(n,m,beam) = Σ_{k,l} w_{k,l}(beam) · exp(j·2π/λ · r̂_{n,m} · d_tx(k,l,0))  (beam-DEPENDENT)
//
// IMPORTANT: This INCOHERENT per-ray summation (Σ_m |h_m|²) differs from the
// full channel coefficient formula (|Σ_m h_m|² — coherent). Using H_usn (which
// stores coherent ray sums) would produce different RSRP values.
//
// This allows us to:
//   Phase 1: Pre-compute antenna patterns and |raysPreComp|² for all rays ONCE
//   Phase 2: Sweep beams — only compute |AF|² per ray (cheap weight × exp sum)
//
BeamSearchResult LINK::find_best_tx_beam(CHANNEL* ch, int bs_idx, int sector_idx)
{
	const int N = ch->NUM_PATH_for_channelcoeff;
	if (N <= 0) {
		BeamSearchResult empty;
		return empty;
	}

	const int K = BS_M / BS_Mp;   // vertical elements per port 0
	const int L = BS_N / BS_Np;   // horizontal elements per port 0
	const int totalRx = MS_M * MS_N * MS_P / (MS_Mp*MS_Np);
	const Real deg2rad = pi / REAL(180.0);
	const Real k_2pi = REAL(2.0) * pi / Wavelength;

	// ================================================================
	// Phase 1: Pre-compute beam-independent ray data (ONCE)
	// Follows TR 36.873 Eq. (8.1-2), (8.1-3) structure
	// ================================================================

	// 1a. Cache port-0 element positions: d_tx(k,l,p=0,Mg=0,Ng=0)
	LOCATION3D port0_pos[8][8];   // [K_MAX][L_MAX]
	for (int k = 0; k < K; k++)
		for (int l = 0; l < L; l++)
			port0_pos[k][l] = bs[bs_idx].d_tx[sector_idx][k][l][0][0][0];

	// 1b. Cache TX departure direction unit vectors per ray
	//     (for computing array factor AF in Phase 2)
	Real sinZoD_cosAoD[MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	Real sinZoD_sinAoD[MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	Real cosZoD       [MAX_NUM_CLUSTERS][MAX_NUM_RAYS];

	// 1b'. Cache RX arrival direction unit vectors per ray (for UE AF)
	Real sinZoA_cosAoA[MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	Real sinZoA_sinAoA[MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	Real cosZoA       [MAX_NUM_CLUSTERS][MAX_NUM_RAYS];

	// 1c. Compute |raysPreComp[polRx][n][m]|² — Eq. (8.1-2) without port virtualization
	//     raysPreComp = [F_rx,u]^T × [PolMatrix] × [F_tx,element(polTx=0)]
	//     This is the beam-INDEPENDENT part of |α_{n,m,u,p}|²
	Real ray_power[2][MAX_NUM_CLUSTERS][MAX_NUM_RAYS];   // [polRx][cluster][ray]
	memset(ray_power, 0, sizeof(ray_power));

	// Handheld: per-port ray power (each port has different element pattern)
	static const int MAX_HANDHELD_PORTS = 8;
	Real hh_ray_power[MAX_HANDHELD_PORTS][MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
	Real hh_los_ray_power[MAX_HANDHELD_PORTS];
	if (handheld_mode) {
		memset(hh_ray_power, 0, sizeof(hh_ray_power));
		memset(hh_los_ray_power, 0, sizeof(hh_los_ray_power));
	}

	for (int n = 0; n < N; n++) {
		int M_n = ch->NUM_RAY_per_ClusterNUM[n];
		for (int m = 0; m < M_n; m++) {
			Real aod_rad = ch->ray_AOD[n][m][0] * deg2rad;
			Real zod_rad = ch->ray_AOD[n][m][1] * deg2rad;
			Real aoa_rad = ch->ray_AOA[n][m][0] * deg2rad;
			Real zoa_rad = ch->ray_AOA[n][m][1] * deg2rad;

			// Cache TX direction unit vector (for AF computation in Phase 2)
			sinZoD_cosAoD[n][m] = sin(zod_rad) * cos(aod_rad);
			sinZoD_sinAoD[n][m] = sin(zod_rad) * sin(aod_rad);
			cosZoD[n][m]        = cos(zod_rad);

			// Cache RX direction unit vector (for UE AF in Phase 2)
			sinZoA_cosAoA[n][m] = sin(zoa_rad) * cos(aoa_rad);
			sinZoA_sinAoA[n][m] = sin(zoa_rad) * sin(aoa_rad);
			cosZoA[n][m]        = cos(zoa_rad);

			// TX element pattern (BS) — polTx=0 for port 0
			Real tx_Ft_P1, tx_Fp_P1, tx_Ft_P2, tx_Fp_P2;
			Get_BS_antenna_pattern(zod_rad, aod_rad, bs_idx, sector_idx,
				tx_Ft_P1, tx_Fp_P1, tx_Ft_P2, tx_Fp_P2);

			// RX element pattern (UE)
			if (handheld_mode) {
				// Handheld: per-port element pattern (2-stage rotation)
				Real inv_sqrt_kappa = REAL(1.0) / sqrt(ch->kappa[n][m]);
				ComplexReal exp_vv(cos(ch->random_phase_vv[n][m] * deg2rad), sin(ch->random_phase_vv[n][m] * deg2rad));
				ComplexReal exp_vh(cos(ch->random_phase_vh[n][m] * deg2rad), sin(ch->random_phase_vh[n][m] * deg2rad));
				ComplexReal exp_hv(cos(ch->random_phase_hv[n][m] * deg2rad), sin(ch->random_phase_hv[n][m] * deg2rad));
				ComplexReal exp_hh(cos(ch->random_phase_hh[n][m] * deg2rad), sin(ch->random_phase_hh[n][m] * deg2rad));

				for (int port = 0; port < handheld_num_ports; port++) {
					Real rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2;
					Get_UE_antenna_pattern(0, zoa_rad, aoa_rad, self_ms_idx, 0,
						rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2, port);
					// Single pol (MS_P=1): use P1 only
					ComplexReal rpc =
						rx_Ft_P1 * exp_vv * tx_Ft_P1 +
						rx_Ft_P1 * inv_sqrt_kappa * exp_vh * tx_Fp_P1 +
						rx_Fp_P1 * inv_sqrt_kappa * exp_hv * tx_Ft_P1 +
						rx_Fp_P1 * exp_hh * tx_Fp_P1;
					hh_ray_power[port][n][m] = std::norm(rpc);
				}
				// Also fill ray_power[0] with port-average for compatibility
				Real sum_rp = 0;
				for (int port = 0; port < handheld_num_ports; port++)
					sum_rp += hh_ray_power[port][n][m];
				ray_power[0][n][m] = sum_rp;
			} else {
				// For omni UE (ue_antenna_element_gain==0): isotropic {1,0}/{0,1}
				// For directional UE: actual pattern via Get_UE_antenna_pattern
				Real rx_Ft[2], rx_Fp[2];
				if (ue_antenna_element_gain == 0) {
					rx_Ft[0] = 1.0; rx_Fp[0] = 0.0;  // polRx=0: θ-polarized
					rx_Ft[1] = 0.0; rx_Fp[1] = 1.0;  // polRx=1: φ-polarized
				} else {
					Real rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2;
					Get_UE_antenna_pattern(0, zoa_rad, aoa_rad, self_ms_idx, 0,
						rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2);
					rx_Ft[0] = rx_Ft_P1; rx_Fp[0] = rx_Fp_P1;
					rx_Ft[1] = rx_Ft_P2; rx_Fp[1] = rx_Fp_P2;
				}

				// Initial random phases and XPR  — Eq. (8.1-2) polarization matrix
				Real inv_sqrt_kappa = REAL(1.0) / sqrt(ch->kappa[n][m]);
				ComplexReal exp_vv(cos(ch->random_phase_vv[n][m] * deg2rad), sin(ch->random_phase_vv[n][m] * deg2rad));
				ComplexReal exp_vh(cos(ch->random_phase_vh[n][m] * deg2rad), sin(ch->random_phase_vh[n][m] * deg2rad));
				ComplexReal exp_hv(cos(ch->random_phase_hv[n][m] * deg2rad), sin(ch->random_phase_hv[n][m] * deg2rad));
				ComplexReal exp_hh(cos(ch->random_phase_hh[n][m] * deg2rad), sin(ch->random_phase_hh[n][m] * deg2rad));

				// raysPreComp[polRx] for polTx=0 (port 0 uses first polarization)
				for (int polRx = 0; polRx < MS_P; polRx++) {
					ComplexReal rpc =
						rx_Ft[polRx] * exp_vv * tx_Ft_P1 +
						rx_Ft[polRx] * inv_sqrt_kappa * exp_vh * tx_Fp_P1 +
						rx_Fp[polRx] * inv_sqrt_kappa * exp_hv * tx_Ft_P1 +
						rx_Fp[polRx] * exp_hh * tx_Fp_P1;
					ray_power[polRx][n][m] = std::norm(rpc);   // |rpc|²
				}
			}
		}
	}

	// 1d. LOS pre-computation — Eq. (8.1-3)
	Real los_ray_power[2] = {0.0, 0.0};
	Real los_sinZoD_cosAoD = 0, los_sinZoD_sinAoD = 0, los_cosZoD = 0;
	Real los_sinZoA_cosAoA = 0, los_sinZoA_sinAoA = 0, los_cosZoA = 0;
	Real K_linear = ch->K_linear;
	bool is_los = (ch->Propagation == LOS_propagation);

	if (is_los) {
		Real los_aod_rad = ch->LOS_AOD_GCS * deg2rad;
		Real los_zod_rad = ch->LOS_ZOD_GCS * deg2rad;
		Real los_aoa_rad = ch->LOS_AOA_GCS * deg2rad;
		Real los_zoa_rad = ch->LOS_ZOA_GCS * deg2rad;

		// LOS TX direction unit vector
		los_sinZoD_cosAoD = sin(los_zod_rad) * cos(los_aod_rad);
		los_sinZoD_sinAoD = sin(los_zod_rad) * sin(los_aod_rad);
		los_cosZoD        = cos(los_zod_rad);

		// LOS RX direction unit vector (for UE AF)
		los_sinZoA_cosAoA = sin(los_zoa_rad) * cos(los_aoa_rad);
		los_sinZoA_sinAoA = sin(los_zoa_rad) * sin(los_aoa_rad);
		los_cosZoA        = cos(los_zoa_rad);

		// LOS antenna patterns
		Real tx_Ft, tx_Fp, tx_Ft2, tx_Fp2;
		Get_BS_antenna_pattern(los_zod_rad, los_aod_rad, bs_idx, sector_idx,
			tx_Ft, tx_Fp, tx_Ft2, tx_Fp2);

		if (handheld_mode) {
			// Per-port LOS ray power
			Real los_phase_rad = ch->random_phase_vv_LOS * deg2rad;
			ComplexReal los_exp(cos(los_phase_rad), sin(los_phase_rad));

			for (int port = 0; port < handheld_num_ports; port++) {
				Real rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2;
				Get_UE_antenna_pattern(0, los_zoa_rad, los_aoa_rad, self_ms_idx, 0,
					rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2, port);
				ComplexReal los_rpc = rx_Ft_P1 * los_exp * tx_Ft - rx_Fp_P1 * los_exp * tx_Fp;
				hh_los_ray_power[port] = std::norm(los_rpc);
			}
			// Also fill los_ray_power[0] with port-sum for compatibility
			Real sum_lrp = 0;
			for (int port = 0; port < handheld_num_ports; port++)
				sum_lrp += hh_los_ray_power[port];
			los_ray_power[0] = sum_lrp;
		} else {
			// RX pattern for LOS: same omni/directional handling as NLOS above
			Real los_rx_Ft[2], los_rx_Fp[2];
			if (ue_antenna_element_gain == 0) {
				los_rx_Ft[0] = 1.0; los_rx_Fp[0] = 0.0;
				los_rx_Ft[1] = 0.0; los_rx_Fp[1] = 1.0;
			} else {
				Real rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2;
				Get_UE_antenna_pattern(0, los_zoa_rad, los_aoa_rad, self_ms_idx, 0,
					rx_Ft_P1, rx_Fp_P1, rx_Ft_P2, rx_Fp_P2);
				los_rx_Ft[0] = rx_Ft_P1; los_rx_Fp[0] = rx_Fp_P1;
				los_rx_Ft[1] = rx_Ft_P2; los_rx_Fp[1] = rx_Fp_P2;
			}

			// LOS polarization matrix: [[e^jΦ, 0]; [0, -e^jΦ]]  — Eq. (8.1-3)
			Real los_phase_rad = ch->random_phase_vv_LOS * deg2rad;
			ComplexReal los_exp(cos(los_phase_rad), sin(los_phase_rad));

			for (int polRx = 0; polRx < MS_P; polRx++) {
				ComplexReal los_rpc =
					los_rx_Ft[polRx] * los_exp * tx_Ft -     // (1,1): +e^jΦ
					los_rx_Fp[polRx] * los_exp * tx_Fp;      // (2,2): -e^jΦ
				los_ray_power[polRx] = std::norm(los_rpc);
			}
		}
	}

	// ================================================================
	// Phase 2: Beam search — only virtualization weights change
	// ================================================================

	Real K_denom = is_los ? (K_linear + REAL(1.0)) : REAL(1.0);  // (K_R + 1)
	Real K_los_num = is_los ? K_linear : REAL(0.0);               // K_R

	// UE port-0 panel-0 element positions (for UE AF in BF case)
	int K_ue = MS_M / MS_Mp;
	int L_ue = MS_N / MS_Np;
	LOCATION3D ue_pos[8][8];
	if (ue_antenna_element_gain != 0 && !handheld_mode) {
		// Runtime guard: ue_pos and the ue_w/ue_v/ue_virtualization_weight_wv tables
		// are 8-deep per dimension, but MAX_MS_M(16) admits larger subarrays. The
		// scenario may resolve ue_antenna_element_gain after cfg validation, so the
		// parse-time check cannot cover every path — enforce here where it matters.
		if (K_ue > 8 || L_ue > 8) {
			cout << "FATAL: directional UE subarray " << K_ue << "x" << L_ue
			     << " exceeds the 8x8 beam-search buffers" << endl;
			exit(1);
		}
		for (int k = 0; k < K_ue; k++)
			for (int l = 0; l < L_ue; l++)
				ue_pos[k][l] = ms[self_ms_idx].d_rx[k][l][0][0][0];  // panel 0
	}

	BeamSearchResult best;
	bool first = true;

	if (handheld_mode)
	{
		// --- Handheld UE: per-port element pattern, no UE beam search ---
		// Each port has a unique element pattern from 2-stage rotation.
		// RSRP = (1/num_ports) × Σ_port |h_port|²
		for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
		{
			for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
			{
				Real alpha_total = 0.0;

				// NLOS: sum over ports
				for (int n = 0; n < N; n++) {
					int M_n = ch->NUM_RAY_per_ClusterNUM[n];
					Real P_n = ch->power[n];
					for (int m = 0; m < M_n; m++) {
						// TX array factor (same for all ports)
						ComplexReal AF(0.0, 0.0);
						for (int k = 0; k < K; k++) {
							for (int l = 0; l < L; l++) {
								Real phase = k_2pi * (
									sinZoD_cosAoD[n][m] * port0_pos[k][l].x +
									sinZoD_sinAoD[n][m] * port0_pos[k][l].y +
									cosZoD[n][m]        * port0_pos[k][l].z);
								// SNS: per-element per-cluster attenuation inside AF sum
								Real sns_w = REAL(1.0);
								if (ch->sns_any_limited) {
									Real pos_h = l * BS_dH;
									Real pos_v = k * BS_dV;
									sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr[n], g_sns_rolloff_C);
								}
								AF += virtualization_weight_wv[z][a][k][l] * sns_w *
									ComplexReal(cos(phase), sin(phase));
							}
						}
						Real AF_sq = std::norm(AF);

						// Sum per-port ray power
						Real sum_port_power = 0.0;
						for (int port = 0; port < handheld_num_ports; port++)
							sum_port_power += hh_ray_power[port][n][m];

						alpha_total += (P_n / (Real)M_n) * AF_sq * sum_port_power;
					}
				}
				alpha_total /= K_denom;

				// LOS component
				if (is_los) {
					ComplexReal AF_LOS(0.0, 0.0);
					for (int k = 0; k < K; k++) {
						for (int l = 0; l < L; l++) {
							Real phase = k_2pi * (
								los_sinZoD_cosAoD * port0_pos[k][l].x +
								los_sinZoD_sinAoD * port0_pos[k][l].y +
								los_cosZoD        * port0_pos[k][l].z);
							// SNS: per-element LOS attenuation inside AF sum
							Real sns_w = REAL(1.0);
							if (ch->sns_any_limited) {
								Real pos_h = l * BS_dH;
								Real pos_v = k * BS_dV;
								sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr_los, g_sns_rolloff_C);
							}
							AF_LOS += virtualization_weight_wv[z][a][k][l] * sns_w *
								ComplexReal(cos(phase), sin(phase));
						}
					}
					Real AF_LOS_sq = std::norm(AF_LOS);

					Real sum_port_los = 0.0;
					for (int port = 0; port < handheld_num_ports; port++)
						sum_port_los += hh_los_ray_power[port];

					alpha_total += (K_los_num / K_denom) * AF_LOS_sq * sum_port_los;
				}

				// Average over ports
				Real rsrp_gain = alpha_total / (Real)handheld_num_ports;
				Real rsrp_dB = 10.0 * log10(rsrp_gain);

				ch->signal_RSRP_gain[sector_idx][z][a][0][0][0] = rsrp_gain;

				if (first || rsrp_dB > best.max_rsrp_dB)
				{
					best.max_rsrp_dB     = rsrp_dB;
					best.sec_azimuth_idx = a;
					best.sec_zenith_idx  = z;
					best.ue_azimuth_idx  = 0;
					best.ue_zenith_idx   = 0;
					best.ue_panel_idx    = 0;
					first = false;
				}
			}
		}
	}
	else if (ue_antenna_element_gain == 0)
	{
		// --- Omni UE: no UE beam search, no UE AF ---
		for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
		{
			for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
			{
				Real alpha = 0.0;

				for (int n = 0; n < N; n++) {
					int M_n = ch->NUM_RAY_per_ClusterNUM[n];
					Real P_n = ch->power[n];
					for (int m = 0; m < M_n; m++) {
						ComplexReal AF(0.0, 0.0);
						for (int k = 0; k < K; k++) {
							for (int l = 0; l < L; l++) {
								Real phase = k_2pi * (
									sinZoD_cosAoD[n][m] * port0_pos[k][l].x +
									sinZoD_sinAoD[n][m] * port0_pos[k][l].y +
									cosZoD[n][m]        * port0_pos[k][l].z);
								// SNS: per-element per-cluster attenuation inside AF sum
								Real sns_w = REAL(1.0);
								if (ch->sns_any_limited) {
									Real pos_h = l * BS_dH;
									Real pos_v = k * BS_dV;
									sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr[n], g_sns_rolloff_C);
								}
								AF += virtualization_weight_wv[z][a][k][l] * sns_w *
									ComplexReal(cos(phase), sin(phase));
							}
						}
						Real AF_sq = std::norm(AF);

						Real sum_u_ray_power = 0.0;
						for (int p = 0; p < MS_P; p++)
							sum_u_ray_power += ray_power[p][n][m];

						sum_u_ray_power *= (Real)(MS_M * MS_N);

						alpha += (P_n / (Real)M_n) * AF_sq * sum_u_ray_power;
					}
				}
				alpha /= K_denom;

				if (is_los) {
					ComplexReal AF_LOS(0.0, 0.0);
					for (int k = 0; k < K; k++) {
						for (int l = 0; l < L; l++) {
							Real phase = k_2pi * (
								los_sinZoD_cosAoD * port0_pos[k][l].x +
								los_sinZoD_sinAoD * port0_pos[k][l].y +
								los_cosZoD        * port0_pos[k][l].z);
							// SNS: per-element LOS attenuation inside AF sum
							Real sns_w = REAL(1.0);
							if (ch->sns_any_limited) {
								Real pos_h = l * BS_dH;
								Real pos_v = k * BS_dV;
								sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr_los, g_sns_rolloff_C);
							}
							AF_LOS += virtualization_weight_wv[z][a][k][l] * sns_w *
								ComplexReal(cos(phase), sin(phase));
						}
					}
					Real AF_LOS_sq = std::norm(AF_LOS);

					Real sum_u_los_power = 0.0;
					for (int p = 0; p < MS_P; p++)
						sum_u_los_power += los_ray_power[p];
					sum_u_los_power *= (Real)(MS_M * MS_N);

					alpha += (K_los_num / K_denom) * AF_LOS_sq * sum_u_los_power;
				}

				Real rsrp_gain = alpha / (Real)totalRx;
				Real rsrp_dB = 10.0 * log10(rsrp_gain);

				ch->signal_RSRP_gain[sector_idx][z][a][0][0][0] = rsrp_gain;

				if (first || rsrp_dB > best.max_rsrp_dB)
				{
					best.max_rsrp_dB     = rsrp_dB;
					best.sec_azimuth_idx = a;
					best.sec_zenith_idx  = z;
					best.ue_azimuth_idx  = 0;
					best.ue_zenith_idx   = 0;
					best.ue_panel_idx    = 0;
					first = false;
				}
			}
		}
	}
	else
	{
		// --- Directional UE: joint TX beam × UE beam search with UE AF ---
		// Panel 0 only. UE AF = Σ_{k,l} ue_wv[ue_z][ue_a][k][l] * exp(j·k_2pi·r̂_rx·d_rx[k][l])

		for (int ue_a = 0; ue_a < ue_tilt_azimuth_angle_LCS_size; ue_a++)
		{
			for (int ue_z = 0; ue_z < ue_tilt_zenith_angle_LCS_size; ue_z++)
			{
				// Pre-compute |AF_rx|² per NLOS ray for this UE beam
				Real ue_af_sq[MAX_NUM_CLUSTERS][MAX_NUM_RAYS];
				for (int n = 0; n < N; n++) {
					int M_n = ch->NUM_RAY_per_ClusterNUM[n];
					for (int m = 0; m < M_n; m++) {
						ComplexReal AF_rx(0.0, 0.0);
						for (int k = 0; k < K_ue; k++) {
							for (int l = 0; l < L_ue; l++) {
								Real phase = k_2pi * (
									sinZoA_cosAoA[n][m] * ue_pos[k][l].x +
									sinZoA_sinAoA[n][m] * ue_pos[k][l].y +
									cosZoA[n][m]        * ue_pos[k][l].z);
								AF_rx += ue_virtualization_weight_wv[ue_z][ue_a][k][l] *
									ComplexReal(cos(phase), sin(phase));
							}
						}
						ue_af_sq[n][m] = std::norm(AF_rx);
					}
				}

				// |AF_rx|² for LOS
				Real ue_los_af_sq = 1.0;
				if (is_los) {
					ComplexReal AF_rx_los(0.0, 0.0);
					for (int k = 0; k < K_ue; k++) {
						for (int l = 0; l < L_ue; l++) {
							Real phase = k_2pi * (
								los_sinZoA_cosAoA * ue_pos[k][l].x +
								los_sinZoA_sinAoA * ue_pos[k][l].y +
								los_cosZoA        * ue_pos[k][l].z);
							AF_rx_los += ue_virtualization_weight_wv[ue_z][ue_a][k][l] *
								ComplexReal(cos(phase), sin(phase));
						}
					}
					ue_los_af_sq = std::norm(AF_rx_los);
				}

				// TX beam search
				for (int a = 0; a < tilt_azimuth_angle_LCS_size; a++)
				{
					for (int z = 0; z < tilt_zenith_angle_LCS_size; z++)
					{
						Real alpha = 0.0;

						// NLOS clusters
						for (int n = 0; n < N; n++) {
							int M_n = ch->NUM_RAY_per_ClusterNUM[n];
							Real P_n = ch->power[n];
							for (int m = 0; m < M_n; m++) {
								// TX array factor
								ComplexReal AF(0.0, 0.0);
								for (int k = 0; k < K; k++) {
									for (int l = 0; l < L; l++) {
										Real phase = k_2pi * (
											sinZoD_cosAoD[n][m] * port0_pos[k][l].x +
											sinZoD_sinAoD[n][m] * port0_pos[k][l].y +
											cosZoD[n][m]        * port0_pos[k][l].z);
										// SNS: per-element per-cluster attenuation inside AF sum
										Real sns_w = REAL(1.0);
										if (ch->sns_any_limited) {
											Real pos_h = l * BS_dH;
											Real pos_v = k * BS_dV;
											sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr[n], g_sns_rolloff_C);
										}
										AF += virtualization_weight_wv[z][a][k][l] * sns_w *
											ComplexReal(cos(phase), sin(phase));
									}
								}
								Real AF_sq = std::norm(AF);

								Real sum_pol_power = 0.0;
								for (int p = 0; p < MS_P; p++)
									sum_pol_power += ray_power[p][n][m];

								alpha += (P_n / (Real)M_n) * AF_sq * ue_af_sq[n][m] * sum_pol_power;
							}
						}
						alpha /= K_denom;

						// LOS component
						if (is_los) {
							ComplexReal AF_LOS(0.0, 0.0);
							for (int k = 0; k < K; k++) {
								for (int l = 0; l < L; l++) {
									Real phase = k_2pi * (
										los_sinZoD_cosAoD * port0_pos[k][l].x +
										los_sinZoD_sinAoD * port0_pos[k][l].y +
										los_cosZoD        * port0_pos[k][l].z);
									// SNS: per-element LOS attenuation inside AF sum
									Real sns_w = REAL(1.0);
									if (ch->sns_any_limited) {
										Real pos_h = l * BS_dH;
										Real pos_v = k * BS_dV;
										sns_w = compute_sns_attenuation(pos_h, pos_v, ch->sns_vr_los, g_sns_rolloff_C);
									}
									AF_LOS += virtualization_weight_wv[z][a][k][l] * sns_w *
										ComplexReal(cos(phase), sin(phase));
								}
							}
							Real AF_LOS_sq = std::norm(AF_LOS);

							Real sum_pol_los_power = 0.0;
							for (int p = 0; p < MS_P; p++)
								sum_pol_los_power += los_ray_power[p];

							alpha += (K_los_num / K_denom) * AF_LOS_sq * ue_los_af_sq * sum_pol_los_power;
						}

						Real rsrp_gain = alpha / (Real)MS_P;
						Real rsrp_dB = 10.0 * log10(rsrp_gain);

						ch->signal_RSRP_gain[sector_idx][z][a][ue_z][ue_a][0] = rsrp_gain;

						if (first || rsrp_dB > best.max_rsrp_dB)
						{
							best.max_rsrp_dB     = rsrp_dB;
							best.sec_azimuth_idx = a;
							best.sec_zenith_idx  = z;
							best.ue_azimuth_idx  = ue_a;
							best.ue_zenith_idx   = ue_z;
							best.ue_panel_idx    = 0;
							first = false;
						}
					}
				}
			}
		}
	}

	return best;
}


// Step 7: Evaluate one candidate cell (one BS + one sector)
CandidateCell LINK::evaluate_candidate_cell(const ScenarioConfig& cfg, int bs_idx, int sector_idx)
{
	CandidateCell cell;
	CHANNEL* ch = &channel[bs_idx][self_ms_idx];

	cell.bs_idx     = bs_idx;
	cell.sector_idx = sector_idx;
	cell.flat_idx   = (cfg.sectors_per_bs == 1) ? bs_idx : bs_idx * 3 + sector_idx;

	// Pathloss is computed per-BS (shared across sectors), so only compute once per BS
	// The caller ensures pathloss_final is already computed before calling this for sector > 0
	cell.pathloss_final    = ch->pathloss_final;
	cell.distance          = ch->distance;
	cell.rms_delay_spread  = ch->RMS_delay_spread;
	cell.aoa_spread        = ch->circular_angle_spread_AOA;
	cell.aod_spread        = ch->circular_angle_spread_AOD;
	cell.is_los            = ch->LOS;
	cell.is_indoor         = ch->Indoor;

	// Beam search
	cell.beam        = find_best_tx_beam(ch, bs_idx, sector_idx);
	cell.max_rsrp_dB = cell.beam.max_rsrp_dB;

	// Signal power
	cell.signal_dBm = bs_maxpower + cell.max_rsrp_dB - cell.pathloss_final;

	return cell;
}


// Step 8: Compute all candidate cells
std::vector<CandidateCell> LINK::compute_all_candidate_cells(const ScenarioConfig& cfg)
{
	std::vector<CandidateCell> candidates;
	candidates.reserve(cfg.num_candidate_bs * cfg.sectors_per_bs);

	for (int wrap_idx = 0; wrap_idx < cfg.num_candidate_bs; wrap_idx++)
	{
		int bs_idx = cfg.use_wraparound
		           ? simple_wrap_mat[cfg.ms_to_bs_wrap_idx][wrap_idx]
		           : wrap_idx;

		// Compute pathloss once per BS (shared across sectors of same BS)
		compute_pathloss_final(cfg, bs_idx);

		for (int sec = 0; sec < cfg.sectors_per_bs; sec++)
		{
			CandidateCell cell = evaluate_candidate_cell(cfg, bs_idx, sec);
			int flat = cell.flat_idx;

			// Store in legacy arrays for backward compatibility
			static_gain[flat].first  = cell.signal_dBm;
			static_gain[flat].second = flat;

			analog_beam_selection[flat].a        = cell.beam.ue_azimuth_idx;
			analog_beam_selection[flat].z        = cell.beam.ue_zenith_idx;
			analog_beam_selection[flat].p        = cell.beam.ue_panel_idx;
			analog_beam_selection[flat].sector_a = cell.beam.sec_azimuth_idx;
			analog_beam_selection[flat].sector_z = cell.beam.sec_zenith_idx;

			candidates.push_back(cell);
		}
	}

	return candidates;
}


// Step 9: Select serving cell from candidates
void LINK::select_serving_cell(const ScenarioConfig& cfg, const std::vector<CandidateCell>& candidates)
{
	int best_idx = -1;
	Real best_signal = -1e30;

	for (int i = 0; i < (int)candidates.size(); i++)
	{
		const CandidateCell& c = candidates[i];

		// single_cell_mode: only bs_idx 0 can serve
		if (cfg.use_single_cell_mode && c.bs_idx != 0)
			continue;

		if (best_idx < 0 || c.signal_dBm > best_signal)
		{
			best_signal = c.signal_dBm;
			best_idx = i;
		}
	}

	const CandidateCell& srv = candidates[best_idx];

	self_bs_idx     = srv.bs_idx;
	self_sector_idx = srv.sector_idx;
	link_pathloss   = srv.pathloss_final;
	link_RSRP       = srv.max_rsrp_dB;
	link_antgain    = srv.max_rsrp_dB;
	link_RMS_delay_spread = srv.rms_delay_spread;
	link_AOA_spread = srv.aoa_spread;
	link_AOD_spread = srv.aod_spread;
	link_distance   = srv.distance;

	azimuth_angle_idx_selected_for_interference = srv.beam.ue_azimuth_idx;
	zenith_angle_idx_selected_for_interference  = srv.beam.ue_zenith_idx;
	panel_idx_selected_for_interference         = srv.beam.ue_panel_idx;

	sector_azimuth_angle_idx = srv.beam.sec_azimuth_idx;
	sector_zenith_angle_idx  = srv.beam.sec_zenith_idx;

	if (cfg.sectors_per_bs == 1)
		_sector_in_control = self_bs_idx;
	else
		_sector_in_control = self_bs_idx * 3 + self_sector_idx;

	if (cfg.store_los_indoor_info)
	{
		link_los    = srv.is_los;
		link_indoor = srv.is_indoor;
	}
}


// Step 10: Compute interference from non-serving cells
void LINK::compute_interference(const ScenarioConfig& cfg, const std::vector<CandidateCell>& candidates)
{
	// Deterministic per-UE RNG (reproducibility fix). Called from Get_signal_interference
	// -> Configuration inside the OpenMP "Get RSRP" loop (Initiallization.cpp). The
	// "random beam interference" tilt draws below (Config B) would otherwise data-race
	// the global `randnum` across threads and consume it by a run-dependent amount —
	// leaving interference (hence Wideband SIR / precoding SINR) non-reproducible even
	// after the O2I / HARQ fixes. Seeding from this UE's identity (base_seed, drop, ue)
	// makes the interfering-beam realization a pure function of the work item. A distinct
	// salt from UE_Initial_Setting keeps the two per-UE streams independent. Draws over
	// candidate BSs are sequential, so their order is fixed.
	Rand randnum = make_workitem_rng((unsigned long long)_seed,
	                                 (unsigned long long)drop_idx,
	                                 (unsigned long long)self_ms_idx,
	                                 0x1A7E7F3ULL /* salt: interfering-beam stream */);

	interference = 0.0;

	// --- InH scenarios ---
	if (!cfg.use_rsrp_table_for_interference)
	{
		if (cfg.configuration_type == 0) // Config A: cumulative from signal phase
		{
			// In Config A, interference = sum of all candidate signals except the serving cell
			for (int i = 0; i < (int)candidates.size(); i++)
			{
				const CandidateCell& c = candidates[i];
				bool is_serving = (c.bs_idx == self_bs_idx) &&
				                  (cfg.sectors_per_bs == 1 || c.sector_idx == self_sector_idx);
				if (!is_serving)
				{
					interference += dBm2linear(c.signal_dBm);
				}
			}
		}
		else  // Config B: random beam re-computation
		{
			int total_bs = cfg.num_candidate_bs;
			for (int wrap_idx = 0; wrap_idx < total_bs; wrap_idx++)
			{
				int bs_idx = cfg.use_wraparound
				           ? simple_wrap_mat[cfg.ms_to_bs_wrap_idx][wrap_idx]
				           : wrap_idx;

				CHANNEL* ch = &channel[bs_idx][self_ms_idx];
				Real pl = ch->pathloss_final;

				for (int sec = 0; sec < cfg.sectors_per_bs; sec++)
				{
					int flat = (cfg.sectors_per_bs == 1) ? bs_idx : bs_idx * 3 + sec;

					// Random beam selection
					int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
					int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

					if (cfg.sectors_per_bs > 1)
					{
						rand_sec_a[flat] = aa;
						rand_sec_z[flat] = zz;
					}

					compute_tx_antenna_gains(ch, bs_idx, sec, aa, zz);
					compute_tx_smallscale_gains(ch, bs_idx, sec, aa, zz);
					Real random_rsrp = compute_rsrp(ch, sec, zz, aa, 1);

					Real sig = bs_maxpower + random_rsrp - pl;

					bool is_serving = (bs_idx == self_bs_idx) &&
					                  (cfg.sectors_per_bs == 1 || sec == self_sector_idx);
					if (!is_serving)
					{
						interference += dBm2linear(sig);
					}

					intf_w_rnd_RSRP[flat] = sig;
				}
			}
		}

		// comp_interf_strength for InH 3TRxP
		if (cfg.sectors_per_bs == 3)
		{
			comp_interf_strength[0] = bs_maxpower + 10 * log10(
				channel[self_bs_idx][self_ms_idx].signal_RSRP_gain[self_sector_idx]
				[sector_zenith_angle_idx][sector_azimuth_angle_idx]
				[comp_azimuth_angle_idx_selected][comp_zenith_angle_idx_selected][comp_panel_idx_selected])
				- channel[self_bs_idx][self_ms_idx].pathloss_final;
		}
	}
	// --- Dense Urban / Rural: RSRP table lookup ---
	else
	{
		// Initialize arrays
		for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		{
			intf_w_rnd_RSRP[bs_idx * 3    ] = linear2dBm(0);
			intf_w_rnd_RSRP[bs_idx * 3 + 1] = linear2dBm(0);
			intf_w_rnd_RSRP[bs_idx * 3 + 2] = linear2dBm(0);

			rand_sec_a[3 * bs_idx    ] = -1;
			rand_sec_a[3 * bs_idx + 1] = -1;
			rand_sec_a[3 * bs_idx + 2] = -1;

			rand_sec_z[3 * bs_idx    ] = -1;
			rand_sec_z[3 * bs_idx + 1] = -1;
			rand_sec_z[3 * bs_idx + 2] = -1;
		}

		// Macro BS interference
		// NOTE: for Dense Urban / Rural this result is recomputed and overwritten by
		// Get_CouplingLoss (Initiallization.cpp) from the rand_sec_a/z stored below.
		// Row-beam mode therefore only changes the lookup there; the draws here stay
		// untouched so the RNG stream is identical with the flag on or off.
		for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		{
			Real pl = channel[bs_idx][self_ms_idx].pathloss_final;

			for (int sec = 0; sec < 3; sec++)
			{
				if (bs_idx == self_bs_idx && sec == self_sector_idx)
					continue;

				int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
				int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

				rand_sec_a[3 * bs_idx + sec] = aa;
				rand_sec_z[3 * bs_idx + sec] = zz;

				Real random_rsrp = linear2dB(
					channel[bs_idx][self_ms_idx].signal_RSRP_gain[sec][zz][aa]
					[zenith_angle_idx_selected_for_interference]
					[azimuth_angle_idx_selected_for_interference]
					[panel_idx_selected_for_interference]);

				Real sig = bs_maxpower + random_rsrp - pl;

				interference += dBm2linear(sig);
				intf_w_rnd_RSRP[bs_idx * 3 + sec] = sig;
			}
		}

	}
}


// Step 11: Main entry point
void LINK::Get_signal_interference()
{
	// Phase 0: Build scenario configuration
	ScenarioConfig cfg = build_scenario_config();

	// Phase 1: Evaluate all candidate cells (beam search + signal power)
	std::vector<CandidateCell> candidates = compute_all_candidate_cells(cfg);

	// Phase 2: Select serving cell (strongest signal)
	select_serving_cell(cfg, candidates);

	// Find the strongest signal value from serving cell
	Real strongest_signal = -1e30;
	for (const auto& c : candidates)
	{
		bool is_serving = (c.bs_idx == self_bs_idx) &&
		                  (cfg.sectors_per_bs == 1 || c.sector_idx == self_sector_idx);
		if (is_serving)
		{
			strongest_signal = c.signal_dBm;
			break;
		}
	}

	// Phase 3: Compute interference
	compute_interference(cfg, candidates);

	// Phase 4: Compute geometry (SIR/SINR)
	str_signal   = strongest_signal;
	signal       = strongest_signal;
	interference = linear2dBm(interference);
	Real noise_linear = dBm2linear(thermal_noise + 10.0 * log10(bandwidth) + MS_noisefig);
	geometry = dBm2linear(signal) / (dBm2linear(interference) + noise_linear);
}
