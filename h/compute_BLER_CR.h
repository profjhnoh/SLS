#include "const.h"  
#include "common.h"


Real bler;
int modulation_type;
Real cr;
Real x_axis;
int cr_target;
int modulation_index;
Real cr_bler_curve[7][2];




Real linear_cr(Real cr1, Real cr1_value, Real cr2, Real cr2_value, Real cr);
Real linear_approx_yais(Real a, Real b, Real c, Real d, Real x);
void bler_computing_setting(LL_CURVE_LOOK);
void bler_computing_target_setting();
void bler_curve_generation();
void bler_computation();

