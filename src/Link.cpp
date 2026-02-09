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
Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2);

ComplexReal  Get_BS_antenna_field_pattern(LOCATION interferer, int _bs_idx, int _sec_idx, int tilt_z_idx, int tilt_a_idx, Real v_angle_theta, Real h_angle_pi, Real F_theta_GCS_P1, Real F_pi_GCS_P1, Real F_theta_GCS_P2, Real F_pi_GCS_P2, ComplexReal  &F_tx_theta, ComplexReal  &F_tx_pi);
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



void LINK::Get_signal_interference()
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

			//pathloss = pathloss + ms[self_ms_idx].LSPs[bs_idx](0);   ////// Add Shadow Fading
			pathloss = pathloss + channel_of_interest->sigma_SF * randnum.n();   ////// Add Shadow Fading
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
			pathloss                            = pathloss + channel_of_interest->sigma_SF * randnum.n();   ////// Add Shadow Fading
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

		if ( g_comp_mode )
		{
			static_gain  [comp_sector_idx].first     = signal;
			static_gain  [comp_sector_idx].second    = comp_sector_idx;

			analog_beam_selection[comp_sector_idx].a = comp_azimuth_angle_idx_selected;
			analog_beam_selection[comp_sector_idx].z = comp_zenith_angle_idx_selected;
			analog_beam_selection[comp_sector_idx].p = comp_panel_idx_selected;

			analog_beam_selection[comp_sector_idx].sector_a = comp_sector_a;
			analog_beam_selection[comp_sector_idx].sector_z = comp_sector_z;
		}

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
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			pathloss                            = pathloss + incar_loss + Otoi_loss + channel_of_interest->sigma_SF * randnum.n();
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
				pathloss                            = pathloss + incar_loss + Otoi_loss + channel_of_interest->sigma_SF * randnum.n();
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
		

		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		{
			for (int mTRP_idx = num_BS; mTRP_idx < num_mTRP + num_BS; mTRP_idx++)
			{
				// 220811 jhnoh
				channel_of_interest                 = &channel[mTRP_idx][self_ms_idx];
				pathloss                            = channel_of_interest->pathloss;
				RMS_delay                           = channel_of_interest->RMS_delay_spread;
				AOA_spread                          = channel_of_interest->circular_angle_spread_AOA;
				AOD_spread                          = channel_of_interest->circular_angle_spread_AOD;
				distance                            = channel_of_interest->distance;
				pathloss                            = channel_of_interest->pathloss_final;

				for (int sector_idx = 0; sector_idx < 3; sector_idx++)
				{
					if ((mTRP_idx != self_bs_idx))
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

						temp_signal_strength = micro_bs_power + Max_RSRP_gain - pathloss;
						if ( temp_signal_strength > comp_signal_strength[1] )
						{						
							Real temp_interf_strength = micro_bs_power + 10*log10(channel_of_interest->signal_RSRP_gain
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
								comp_sector_idx                   = mTRP_idx * 3 + (sector_idx);
							}
						}
					} // Exclude Serving BS & Sector 
				} // Sector index loop end
			} // BS index loop end

		}

		if ( g_comp_mode )
		{
			static_gain  [comp_sector_idx].first     = signal;
			static_gain  [comp_sector_idx].second    = comp_sector_idx;

			analog_beam_selection[comp_sector_idx].a = comp_azimuth_angle_idx_selected;
			analog_beam_selection[comp_sector_idx].z = comp_zenith_angle_idx_selected;
			analog_beam_selection[comp_sector_idx].p = comp_panel_idx_selected;

			analog_beam_selection[comp_sector_idx].sector_a = comp_sector_a;
			analog_beam_selection[comp_sector_idx].sector_z = comp_sector_z;
		}
		

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

		if ( g_mTRP_mode == 1 || g_mTRP_mode == 2) 
		{
			for (int mTRP_idx = num_BS; mTRP_idx < num_BS + num_mTRP; mTRP_idx++)
			{
				pathloss = channel[mTRP_idx][self_ms_idx].pathloss_final;
				channel_of_interest = &channel[mTRP_idx][self_ms_idx];

				//--------------------------- random beam interference ---------------------------
				intf_w_rnd_RSRP[mTRP_idx*3  ] = linear2dBm(0);
				intf_w_rnd_RSRP[mTRP_idx*3+1] = linear2dBm(0);
				intf_w_rnd_RSRP[mTRP_idx*3+2] = linear2dBm(0);

				rand_sec_a[3*mTRP_idx  ] = -1;
				rand_sec_a[3*mTRP_idx+1] = -1;
				rand_sec_a[3*mTRP_idx+2] = -1;

				rand_sec_z[3*mTRP_idx  ] = -1;
				rand_sec_z[3*mTRP_idx+1] = -1;
				rand_sec_z[3*mTRP_idx+2] = -1;				

				if ( mTRP_idx == self_bs_idx )
				{
				}
				else
				{
					int sector_idx = (int) (3 * randnum.u());
					int aa = (int)(tilt_azimuth_angle_LCS_size * randnum.u());
					int zz = (int)(tilt_zenith_angle_LCS_size * randnum.u());

					rand_sec_a[3*mTRP_idx+sector_idx] = aa;
					rand_sec_z[3*mTRP_idx+sector_idx] = zz;

					//--------------------------------------------------------------------------------
					random_RSRP_antgain = 
					linear2dB(channel[mTRP_idx][self_ms_idx].signal_RSRP_gain[sector_idx][zz][aa]
								[zenith_angle_idx_selected_for_interference]
								[azimuth_angle_idx_selected_for_interference]
								[panel_idx_selected_for_interference]);

					signal = micro_bs_power + random_RSRP_antgain - pathloss;
					interference += dBm2linear(signal);
					intf_w_rnd_RSRP[mTRP_idx*3 + sector_idx] = signal;
				}
			}

			Real test_interference = 0;
			for(int sec_idx = 0; sec_idx < num_SECTORS + num_mTRP_SECTORS; sec_idx++)
			{
				if( sec_idx == _sector_in_control )
				{
				}
				else
				{
					test_interference += dBm2linear(intf_w_rnd_RSRP[sec_idx]);
				}
			}
			
			if ( interference != test_interference)
				 cout << "Something wrong with interference" << endl;
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

			//pathloss = pathloss + ms[self_ms_idx].LSPs[bs_idx](0) + incar_loss + Otoi_loss;   ////// Add Shadow Fading
			pathloss = pathloss + incar_loss + Otoi_loss + channel_of_interest->sigma_SF * randnum.n();

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

void LINK::Get_adj_SECTORS() 
{
	Real pathlosss;
	Real antgainn;

	// jhnoh 220809
	if ( g_mTRP_mode == 1 || g_mTRP_mode == 2)
		sort(static_gain,static_gain + num_SECTORS + num_mTRP_SECTORS, comparator); 
	else
		sort(static_gain,static_gain + num_SECTORS, comparator); 


	for (int idx = 0; idx < num_SECTORS; idx++)
		adj_sector[idx] = static_gain[idx].second;

	//Real static_gain_ratio = static_gain[0].first - static_gain[1].first;
	//if ( static_gain_ratio < 5 && linear2dB(geometry) < 10 && g_comp_mode )
	//	comp_mode[self_ms_idx] = 1;
	//else
	//	comp_mode[self_ms_idx] = 0;

	// jhnoh 230109 
	// Consider that ue select the other Rx antenna pannel for avoiding inter-layer interference
	// adj_sector[1] = comp_sector_idx;
}


void LINK::UE_Initial_Setting(void) 
{

  /*------------------------ Car penetration loss ----------------------------*/
  if (TYPE == 12) // Dense Urban 
  {
	  if (ms[self_ms_idx].Indoor == true)   ///// Indoor
	  {
		  incar_loss = 0.;
	  }
	  else if (ms[self_ms_idx].Indoor == false)                              ////// outdoor
	  {
		  incar_loss = 9. + 5. * randnum.n();
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
	  if (ms[self_ms_idx].Indoor == false)   ///// Outdoor
	  {
		  Otoi_loss = 0;
	  }
	  else                               /////// Indoor
	  {
		  if (Channel_Model_Type == 1)   ///// Model B
		  {
			  if (randnum.u() > 0.5) //// 80% Low Loss  O2I
			  //if (true) //// 100% Low Loss  O2I TEST
			  {
				  /////// Low-loss model

				  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
				  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));

				  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();

				  //Otoi_loss = 0.;

			  }
			  else           ////// 20& High Loss O2I
			  {
				  /////// High-loss model

				  Otoi_loss = 5 - 10 * log10(0.7 * pow(10, -1 * (IRRglass_L / 10)) + 0.3 * pow(10, -1 * (concrete_L / 10)));
				  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));

				  Otoi_loss = Otoi_loss + indoor_L + 6.5 * randnum.n();
			  }
		  }
		  else if (Channel_Model_Type == 0 || Channel_Model_Type == 2)  //// Model A 
		  {
			  if (carrier_freq > 6000000000.)  //// above 6GHz
			  {
				  if (randnum.u() > 0.2) //// 80% Low Loss  O2I
				  //if (true) //// 100% Low Loss  O2I    TEST
				  {
					  /////// Low-loss model

					  Otoi_loss = 5 - 10 * log10(0.3 * pow(10, -1 * (glass_L / 10)) + 0.7 * pow(10, -1 * (concrete_L / 10)));
					  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));
					  Otoi_loss = Otoi_loss + indoor_L + 4.4 * randnum.n();
					  //Otoi_loss = 0.;

					  high_loss_flag = 0;
				  }
				  else
				  {
					  /////// High-loss model

					  Otoi_loss = 5 - 10 * log10(0.7 * pow(10, -1 * (IRRglass_L / 10)) + 0.3 * pow(10, -1 * (concrete_L / 10)));
					  indoor_L = 0.5 * MIN((25 * randnum.u()), (25 * randnum.u()));
					  Otoi_loss = Otoi_loss + indoor_L + 6.5 * randnum.n();

					  high_loss_flag = 1;
				  }
			  }
			  else                         //// below 6GHz
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

	Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, _sector_idx, tilt_z_idx, tilt_a_idx, 
	v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi);

	TransmitterAntennaGainXLOS_theta = F_tx_theta;
	TransmitterAntennaGainXLOS_pi = F_tx_pi;

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

				Real h_angle_pi    = channel_of_interest->AOD[i] + channel_of_interest->cluster_ASD * channel_of_interest->offset_angle_rand_coupling[i][j];
				Real v_angle_theta = channel_of_interest->ZOD[i] + (3. / 8.) * pow(10, channel_of_interest->mu_ZSD) * channel_of_interest->offset_angle_rand_coupling[i][j];

				Transform_angle_minus_180_to_plus_180(h_angle_pi);
				Transform_angle_0_to_plus_180(v_angle_theta);

				h_angle_pi = h_angle_pi * (pi / 180.);
				v_angle_theta = v_angle_theta * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_tx_theta = {0, 0};
				ComplexReal  F_tx_pi = {0, 0};

				Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

				Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi);


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
				Real h_angle_pi    = channel_of_interest->AOD[i] + channel_of_interest->cluster_ASD * channel_of_interest->offset_angle_rand_coupling[i][j];
				Real v_angle_theta = channel_of_interest->ZOD[i] + (3. / 8.) * pow(10, channel_of_interest->mu_ZSD) * channel_of_interest->offset_angle_rand_coupling[i][j];

				Transform_angle_minus_180_to_plus_180(h_angle_pi);
				Transform_angle_0_to_plus_180(v_angle_theta);

				h_angle_pi = h_angle_pi * (pi / 180.);
				v_angle_theta = v_angle_theta * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_tx_theta = {0, 0};
				ComplexReal  F_tx_pi = {0, 0};

				Real combined_antenna_gain = Get_BS_antenna_pattern(v_angle_theta, h_angle_pi, _bs_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

				Get_BS_antenna_field_pattern(channel_of_interest->self_bs, _bs_idx, sector_idx, tilt_z_idx, tilt_a_idx, v_angle_theta, h_angle_pi, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2, F_tx_theta, F_tx_pi);

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
				Real h_angle_pi = channel_of_interest->AOA[i] + channel_of_interest->cluster_ASA * channel_of_interest->offset_angle[i][j];
				Real v_angle_theta = channel_of_interest->ZOA[i] + channel_of_interest->cluster_ZSA * channel_of_interest->offset_angle[i][j];

				Transform_angle_minus_180_to_plus_180(h_angle_pi);
				Transform_angle_0_to_plus_180(v_angle_theta);


				h_angle_pi = h_angle_pi * (pi / 180.);
				v_angle_theta = v_angle_theta * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_rx_theta   = {0, 0};
				ComplexReal  F_rx_pi      = {0, 0};
				ComplexReal  F_rx_theta_2 = {0, 0};
				ComplexReal  F_rx_pi_2    = {0, 0};

				Real combined_antenna_gain;


				if (ue_antenna_element_gain == 0) // UE omni-antenna case
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

					combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

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

		//
		for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
		{
			for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
			{
				Real h_angle_pi = channel_of_interest->AOA[i] + channel_of_interest->cluster_ASA * channel_of_interest->offset_angle[i][j];
				Real v_angle_theta = channel_of_interest->ZOA[i] + channel_of_interest->cluster_ZSA * channel_of_interest->offset_angle[i][j];

				Transform_angle_minus_180_to_plus_180(h_angle_pi);
				Transform_angle_0_to_plus_180(v_angle_theta);

				h_angle_pi = h_angle_pi * (pi / 180.);
				v_angle_theta = v_angle_theta * (pi / 180.);

				Real F_theta_GCS_P1 = 0;
				Real F_pi_GCS_P1 = 0;
				Real F_theta_GCS_P2 = 0;
				Real F_pi_GCS_P2 = 0;

				ComplexReal  F_rx_theta   = {0, 0};
				ComplexReal  F_rx_pi      = {0, 0};
				ComplexReal  F_rx_theta_2 = {0, 0};
				ComplexReal  F_rx_pi_2    = {0, 0};

				Real combined_antenna_gain;

				if (ue_antenna_element_gain == 0) // UE omni-antenna case
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
				else  // UE directional antenna case
				{

					combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

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
	complex<Real> Big_PI_LOS             (0. , (360. * randnum.u() - 180.0)* (pi / 180.0));

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
					}


					for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
					{
						for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
						{
							kappa = channel_of_interest->kappa[i][j];
							_1_over_sqrt_K = 1.0 / sqrt(kappa);

							complex<Real> Big_pi_NLOS_thetatheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
							complex<Real> Big_pi_NLOS_thetapi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
							complex<Real> Big_pi_NLOS_pitheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
							complex<Real> Big_pi_NLOS_pipi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));

							alpha_nmup_temp = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
								(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
								(SmallScale_RX_AntennaGainXLOS_theta(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
								);

							alpha_nmup = alpha_nmup + (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));
						}
					}
					alpha += alpha_zero + alpha_nmup;
				}
			}
		}

		channel_of_interest->signal_RSRP_gain[sec_number][sec_z_idx][sec_a_idx][0][0][0] = alpha / (Real)(M*N*P);


		selected_a = 0;
		selected_z = 0;
		selected_p = 0;

		Real RSRP_antenna_gain = (alpha) / (M*N*P);

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
							//complex<Real> Big_PI_LOS(0., (360. * randnum.u() - 180.0)* (pi / 180.0));

							// panel 1
							alpha_zero_temp = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
							alpha_zero = abs(alpha_zero_temp)*abs(alpha_zero_temp);

							// panel 2
							alpha_zero_temp_panel_2 = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
							alpha_zero_panel_2 = abs(alpha_zero_temp_panel_2)*abs(alpha_zero_temp_panel_2);
						}

						for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
						{
							for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
							{
								kappa = channel_of_interest->kappa[i][j];
								_1_over_sqrt_K = 1.0 / sqrt(kappa);

								complex<Real> Big_pi_NLOS_thetatheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
								complex<Real> Big_pi_NLOS_thetapi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
								complex<Real> Big_pi_NLOS_pitheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
								complex<Real> Big_pi_NLOS_pipi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));

								// panel 1
								alpha_nmup_temp = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
									(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
									(SmallScale_RX_AntennaGainXLOS_theta(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
									);

								alpha_nmup = alpha_nmup + (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));

								/// panel 2
								alpha_nmup_temp_panel_2 = sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
									(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)* _1_over_sqrt_K *FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
									(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
									);

								alpha_nmup_panel_2 = alpha_nmup_panel_2 + (abs(alpha_nmup_temp_panel_2) * abs(alpha_nmup_temp_panel_2));

							}
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
									//complex<Real> Big_PI_LOS(0., (360. * randnum.u() - 180.0)* (pi / 180.0));

									// panel 1
									alpha_zero_temp = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
									alpha_zero = abs(alpha_zero_temp)*abs(alpha_zero_temp);

									// panel 2
									alpha_zero_temp_panel_2 = sqrt(K_linear / (K_linear + 1)) * (ReceiverAntennaGainXLOS_theta_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_theta - ReceiverAntennaGainXLOS_pi_panel_2 * FAST_EXP(Big_PI_LOS) * TransmitterAntennaGainXLOS_pi);
									alpha_zero_panel_2 = abs(alpha_zero_temp_panel_2)*abs(alpha_zero_temp_panel_2);
								}

								for (int i = 0; i < channel_of_interest->NUM_PATH_for_channelcoeff; i++)
								{
									for (int j = 0; j < channel_of_interest->NUM_RAY_per_ClusterNUM[i]; j++)
									{
										kappa = channel_of_interest->kappa[i][j];
										_1_over_sqrt_K = 1.0 / (Real)sqrt(kappa);

										complex<Real> Big_pi_NLOS_thetatheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
										complex<Real> Big_pi_NLOS_thetapi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
										complex<Real> Big_pi_NLOS_pitheta(0, (360. * randnum.u() - 180.0)*(pi / 180.0));
										complex<Real> Big_pi_NLOS_pipi(0, (360. * randnum.u() - 180.0)*(pi / 180.0));

										// panel 1
										alpha_nmup_temp = (Real)sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
											(SmallScale_RX_AntennaGainXLOS_theta(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
											(SmallScale_RX_AntennaGainXLOS_theta(i, j) * (Real)_1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi(i, j)* FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
											);

										alpha_nmup = alpha_nmup + (abs(alpha_nmup_temp) * abs(alpha_nmup_temp));

										/// panel 2
										alpha_nmup_temp_panel_2 = (Real)sqrt(channel_of_interest->power[i] / (channel_of_interest->NUM_RAY_per_ClusterNUM[i] * (K_linear + 1))) * (
											(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * FAST_EXP(Big_pi_NLOS_thetatheta) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)* _1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_pitheta)) * SmallScale_TX_AntennaGainXLOS_theta(i, j) +
											(SmallScale_RX_AntennaGainXLOS_theta_panel_2(i, j) * (Real)_1_over_sqrt_K * FAST_EXP(Big_pi_NLOS_thetapi) + SmallScale_RX_AntennaGainXLOS_pi_panel_2(i, j)*FAST_EXP(Big_pi_NLOS_pipi)) * SmallScale_TX_AntennaGainXLOS_pi(i, j)
											);

										alpha_nmup_panel_2 = alpha_nmup_panel_2 + (abs(alpha_nmup_temp_panel_2) * abs(alpha_nmup_temp_panel_2));

									}
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
	return acosl( cos(beta)*cos(gamma)*cos(GCS_theta) + (sin(beta)*cos(gamma)*cos(GCS_pi-alpha) - sin(gamma)*sin(GCS_pi-alpha))*sin(GCS_theta));
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


Real LINK::Get_MS_antgain(int M, int N, int P, CHANNEL * channel_of_interest , int sector_idx, int tilt_z_idx, int tilt_a_idx, int ue_interference_index)
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


	if (ue_antenna_element_gain == 0) // UE omni-antenna case
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

		combined_antenna_gain = Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx, F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2);

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

		Real zeta_1 = -pi / 4;
		Real zeta_2 = pi / 4;

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

		Real zeta_1 = -pi / 4;
		Real zeta_2 = pi / 4;

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

		Real zeta_1 = -pi / 4;
		Real zeta_2 = pi / 4;

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
	ComplexReal  &F_tx_pi)
{
	Real dH;
	Real dV;

	dH = BS_dH;
	dV = BS_dV;

	int K = BS_M / BS_Mp;    /////   4 / 4 = 1
	int L = BS_N / BS_Np;   //////   4/4 = 1

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

			#if COMPARE_PHASE_CALCULATION
			// Original version: using exp()
			weight += w * exp(jay * phase_factor * REAL(dot(r_tx, d_tx)));
			#elif USE_OPTIMIZED_PHASE
			// Optimized version: Use Euler's formula directly
			// exp(i*theta) = cos(theta) + i*sin(theta)
			Real theta = phase_factor * REAL(dot(r_tx, d_tx));
			Real cos_theta = cos(theta);
			Real sin_theta = sin(theta);
			weight += w * ComplexReal(cos_theta, sin_theta);
			#else
			// Original version: using exp()
			weight += w * exp(jay * phase_factor * REAL(dot(r_tx, d_tx)));
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

			// Optimized version: Use Euler's formula directly
			// exp(i*theta) = cos(theta) + i*sin(theta)
			Real theta = phase_factor * REAL(dot(r_tx, d_tx));
			Real cos_theta = cos(theta);
			Real sin_theta = sin(theta);
			weight_optimized += w * ComplexReal(cos_theta, sin_theta);
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


Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS, int ms_idx, int sector_index, Real &F_theta_GCS_P1, Real &F_pi_GCS_P1, Real &F_theta_GCS_P2, Real &F_pi_GCS_P2)
{
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

		Real zeta_1 = 0;
		Real zeta_2 = pi / 2;

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

		if (carrier_freq < 6000000000)  //below 6Ghz
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

		Real zeta_1 = 0;
		Real zeta_2 = pi / 2;

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

	cos_psi = (cos(beta)*cos(gamma)*sin(theta) - (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*cos(theta)) / sqrt(1.0 - pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2));
	sin_psi = (sin(beta)*cos(gamma)*sin(phi - alpha) + sin(gamma)*cos(phi - alpha)) / sqrt(1.0 - pow(cos(beta)*cos(gamma)*cos(theta) + (sin(beta)*cos(gamma)*cos(phi - alpha) - sin(gamma)*sin(phi - alpha))*sin(theta), 2));;

	GCS_field_theta = cos_psi * LCS_field_theta - sin_psi * LCS_field_pi;
	GCS_field_pi = sin_psi * LCS_field_theta + cos_psi * LCS_field_pi;
}


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




ANTENNA_FIELD Antenna_field_in_LCS(Real LCS_antenna_gain_in_dB, Real zeta)
{
	ANTENNA_FIELD field_pattern;
	Real LCS_antenna_gain = dB2linear(LCS_antenna_gain_in_dB);

	field_pattern.theta = pow(LCS_antenna_gain, 0.5)*cos(zeta);
	field_pattern.phi = pow(LCS_antenna_gain, 0.5)*sin(zeta);

	return field_pattern;
}

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
