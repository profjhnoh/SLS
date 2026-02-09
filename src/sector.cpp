#include "common.h"



Sector::Sector()
{
}



Sector::~Sector()
{
}


void Sector::Reset2Default()
{
	schedule_done_flag        = 0;
	self_bs_idx_schedule_part = 0;
	_N_pf                     = 0;
	num_ue_scheduled_mumimo   = 0;
	ue_scheduled_mumimo.resize(0);

    W                         = NULL;
    PMI_read                  = NULL;
    PMI_vector_read           = NULL;
    CQI_read                  = NULL;
    MCS_decision              = NULL;
    CQIndex_decision          = NULL;
    CQI_AVR                   = NULL;
    METRIC                    = NULL;

	ue_in_control.clear();
	ue_in_comp.clear();
	selected_ue.clear();
	scheduled_ue.clear();
}
