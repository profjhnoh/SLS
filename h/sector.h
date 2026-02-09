#include "const.h"
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class Sector
{


public:
	Sector();
	~Sector();

	void printCombination(VectorXd, VectorXd, Real, int, int);
	void combinationUtil(VectorXd, VectorXd, VectorXd, Real, int, int, int, int r);
	void Read_Ch_Feedback(void);
	void Read_CSI_Feedback(void);  // TDD: Read full channel matrices from UEs
	int  determine_MCS(Real sinr_estimate);
	int  determine_CQI(Real sinr_estimate);
	void Initialize_schedule_map(void);
	void Set_AVR_Cqi(void);
	void Set_AVR_Cqi_NCJT( int );
	void Set_AVR_Cqi_Precoding_Based(void);  // Precoding-based SINR/CQI/MCS calculation
	void Scheduling_RoundRobin(void);
	void Scheduling_PF_MU_MIMO(void);
	void Initial_Setting(int);
	void first_max_ue(int rb_idx);
	void SUS_procedure(int rb_idx);
	void Get_User_Set(int rb_idx);
	void Filter_User_Set(int rb_idx, int newly_selected_ue);

	// Chordal Distance based scheduling
	void ChordalDistance_Scheduling(int rb_idx);
	void ChordalDistance_Scheduling_threadsafe(int rb_idx, vector<int>& local_selected_ue,
	                                           vector<int>& local_scheduled_ue);
	Real Compute_Chordal_Distance(const MatrixXcReal& P_projection, const MatrixXcReal& H_candidate);

	// Thread-safe versions with explicit state parameters
	void first_max_ue_threadsafe(int rb_idx, vector<int>& local_selected_ue);
	void SUS_procedure_threadsafe(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_scheduled_ue, vector<int>& local_user_set);
	void Get_User_Set_threadsafe(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_user_set);
	void Filter_User_Set_threadsafe(int rb_idx, int newly_selected_ue, vector<int>& local_user_set);

	// Optimized versions with pre-computed eigenvectors (cache)
	void Get_User_Set_threadsafe_cached(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_user_set,
	                                      const vector<VectorXcReal>& eigenvector_cache);
	void Filter_User_Set_threadsafe_cached(int rb_idx, int newly_selected_ue, vector<int>& local_user_set,
	                                         const vector<VectorXcReal>& eigenvector_cache);
	void SUS_procedure_threadsafe_cached(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_scheduled_ue,
	                                       vector<int>& local_user_set, const vector<VectorXcReal>& eigenvector_cache);
	void Transmit_Precoding(int rb_idx);
	void Transmit_Precoding_TDD(int rb_idx, const vector<int>& local_scheduled_ue);  // TDD: Reciprocity-based precoding using full channel matrices
	void SU_MIMO_Precoding(int, int);
	void Set_Min_Cqi(void);
	void Set_AVR_Cqi_NCJT(SCHEDULE_DECISION);
	void Scheduling_algorithm_module_MU_MIMO(void);
	void Scheduling_algorithm_module_RoundRobin(void);
	void AVR_Cqi_update(void);
	void Scheduling_mTRP( int );
	void Scheduling_NCJT( int );
	void Reset2Default( void );

	int                  self_idx                  = 0;
	bool                 schedule_done_flag        = 0;
	int                  self_bs_idx_schedule_part = 0;
	Real                 _N_pf                     = 0;
	int                  num_ue_scheduled_mumimo   = 0;
	VectorXd             ue_scheduled_mumimo          ;

	PMI_FEEDBACK **      PMI_read                  =  NULL;
	Real **              CQI_read                  =  NULL;
	Real **              CQI_comp_read             =  NULL;
	Real **              CQI_AVR                   =  NULL;
	Real **              METRIC                    =  NULL;
	int    **            MCS_decision              =  NULL;
	int    **            CQIndex_decision          =  NULL;
	VectorXcReal **      PMI_vector_read           =  NULL;
	MatrixXcReal **      CSI_matrix_read           =  NULL;  // TDD: Full channel matrices from UEs
	MatrixXcReal *       W                         =  NULL;
	SCHEDULE_DECISION **   ppschedulewrite           =  NULL;

	vector <int> ue_in_control;
	vector <int> ue_in_comp;
	vector <int> selected_ue;
	vector <int> scheduled_ue;

	vector <int> mTRP_in_control;
	
private:
	vector <int> user_set;
	Real prev_max_sinr   = 0;
	Real prev_sum_SINR   = 0;
	Real sum_of_SINR     = 0;
};
