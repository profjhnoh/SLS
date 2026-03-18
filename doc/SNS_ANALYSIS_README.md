# SNS (Spatial Non-Stationarity) Attenuation Analysis - Complete Report

**Date:** 2026-03-14
**Status:** Problem Identified and Documented
**Files Created:** 3 analysis documents

---

## Quick Summary

The SNS attenuation formula in `compute_sns_attenuation()` (src/channel.cpp:5938) uses:

```
gamma = exp(-d / (D_n * C))
```

where `D_n ≈ 8.94λ` and `C = 13.0`, making the denominator ~116.

**Problem:** This is 7-8x weaker than it should be.

**Root Cause:** The D_n normalization stretches the roll-off distance by D_n factor.

**Expected Behavior:** Should be `gamma = exp(-d / C)` without D_n.

**Impact:** SNS attenuation is negligible (-0.09 dB avg) instead of meaningful (-0.7 dB avg).

---

## Documents Included

### 1. SNS_Attenuation_Analysis.md (PRIMARY)
**Comprehensive technical analysis covering:**
- Array dimensions for Config 3 UMi (15.5λ × 16.1λ)
- V_n calculations (typical: 0.56-0.60)
- VR dimensions (typical: 9.3λ × 9.7λ)
- D_n values (typical: 8.94λ)
- Attenuation formula comparison (Formula A vs B)
- Numerical results across 768 BS elements
- Detailed problem assessment
- Calibration impact analysis
- Recommended fixes

**Key Finding:**
```
Formula A (current):     gamma_avg = 0.9796  →  -0.0894 dB
Formula B (correct):     gamma_avg = 0.8445  →  -0.7341 dB
Difference:              8.1x weaker than expected
```

### 2. SNS_Formula_Technical_Summary.md (TECHNICAL DETAILS)
**Deep technical analysis covering:**
- Current implementation with code snippets
- Why D_n normalization is wrong
- Problem 1: D_n scaling effect
- Problem 2: VR-size dependency issue
- Problem 3: Maximum attenuation negligible
- Quantitative impact for Config 3
- Element-by-element attenuation table
- Coupling loss impact on calibration
- Spec interpretation issues
- Hypothesis for correct formula
- Proposed code fix with diff
- Files affected by change

**Key Data:**
```
Elements at boundary:    Formula A: -0.02 dB  |  Formula B: -0.23 dB  (10x difference)
Elements at far corner:  Formula A: -0.33 dB  |  Formula B: -2.90 dB  (9x difference)
```

### 3. SNS_Code_Flow.md (IMPLEMENTATION)
**Complete code flow and call chain:**
- High-level flow diagram
- Detailed function-by-function breakdown:
  - Generate_VisibilityRegion() [Line 3344]
  - sns_setup_corner() [Line 3312]
  - Compute_SNS_RSRP_Attenuation() [Line 3451]
  - compute_sns_attenuation() [Line 5927] ← PROBLEM HERE
- Data structure ClusterVR definition
- Global variables used
- Calling context and dependencies
- Impact on calibration workflow
- Configuration-specific effects

**Code Location:**
```
File:     src/channel.cpp
Function: compute_sns_attenuation()
Line:     5938  (problem location)
Problem:  Real denom = vr.D_n * rolloff_C;  // Should be: Real denom = rolloff_C;
```

---

## Numerical Evidence

### Maximum Attenuation at Far Corner

```
Config 3 UMi parameters:
  BS array:     24×32 elements
  D_n value:    8.94λ  (large array)
  Rolloff C:    13.0

Current formula: gamma = exp(-8.94 / (8.94*13)) = exp(-1/13) = 0.926
                 Attenuation = 10*log10(0.926) = -0.334 dB

Expected formula: gamma = exp(-8.94 / 13) = exp(-0.688) = 0.513
                 Attenuation = 10*log10(0.513) = -2.90 dB

Ratio: 2.90 / 0.334 = 8.7x difference
```

### Average Attenuation Across All 768 Elements

```
Formula A (current, with D_n normalization):
  Average: 0.9796 (power)  →  -0.0894 dB
  Interpretation: Negligible effect, barely affects coupling loss

Formula B (proposed, without D_n):
  Average: 0.8445 (power)  →  -0.7341 dB
  Interpretation: Meaningful effect, measurable in calibration

Difference: 8.1x stronger with corrected formula
```

### Element Distribution

Of 768 BS elements:
- Formula A: 429 elements (55.9%) get <0.1 dB loss
- Formula B: Only 300 elements (39.1%) get <0.1 dB loss
- Formula B provides 50% of elements with >0.3 dB loss (important effect)

---

## Validation Checklist

To confirm this is a bug:

- [ ] Check j20 Section 7.6.14.1.3 specification
  - [ ] Exact formula for attenuation
  - [ ] Role of D_n parameter
  - [ ] Meaning of C parameter
  
- [ ] Compare with reference implementations
  - [ ] ns-3 CDL model
  - [ ] 3GPP reference simulator
  - [ ] Verify if they use D_n normalization

- [ ] Empirical validation
  - [ ] Run Config 3 calibration with Formula B
  - [ ] Check if coupling loss matches targets better
  - [ ] Verify RSRP improvements in test cases

- [ ] Physical reasoning
  - [ ] Does -0.09 dB SNS effect make sense? (No)
  - [ ] Does -0.73 dB SNS effect make sense? (Yes)
  - [ ] Does D_n-dependent attenuation make sense? (No)

---

## Recommended Next Steps

1. **Verify Specification** (HIGHEST PRIORITY)
   - Obtain j20 Section 7.6.14.1.3 detailed tables
   - Confirm exact attenuation formula
   - Determine if D_n appears in denominator

2. **Check Reference Implementation**
   - Find how ns-3 or 3GPP simulator implements SNS
   - Compare formula with current code

3. **Implement Tentative Fix** (IF SPEC UNCLEAR)
   - Change line 5938: `Real denom = rolloff_C;`
   - Re-run Config 3 calibration
   - Check if coupling loss improves

4. **Document Finding**
   - Add comment explaining correct formula
   - Reference j20 section and equation number
   - Note any lessons learned

---

## Files Referenced

### Source Code
- `/src/channel.cpp` - Problem location and related functions
- `/h/channel.h` - ClusterVR structure definition
- `/h/common.h` - Global SNS parameters

### Analysis Documents (Created)
- `/doc/SNS_Attenuation_Analysis.md` - Comprehensive analysis
- `/doc/SNS_Formula_Technical_Summary.md` - Technical details
- `/doc/SNS_Code_Flow.md` - Implementation flow

---

## Contact & Notes

**Analysis performed:** 2026-03-14
**Analyst:** Claude Code
**Issue Category:** SNS Formula Interpretation
**Severity:** High (affects calibration accuracy by ~0.6-0.7 dB)
**Status:** Awaiting specification review and validation

**Key Finding:** The formula appears to have a mathematical error that makes SNS effects 7-8x weaker than intended. This is likely a misinterpretation of the j20 specification regarding the role of D_n in the attenuation calculation.

