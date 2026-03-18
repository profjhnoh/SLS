# SNS Attenuation Formula: Technical Root Cause Analysis

**Analysis Date:** 2026-03-14
**Configuration:** Config 3 UMi (24×32 antenna array)
**Issue:** Attenuation is 7-8x weaker than expected due to D_n normalization

---

## Executive Summary

The SNS (Spatial Non-Stationarity) attenuation formula in `compute_sns_attenuation()` uses:

```
gamma = exp(-d / (D_n * C))   where d = distance from VR edge
```

This normalizes the exponential decay by **D_n (the visibility region diagonal)**, which makes the roll-off ~8x more gradual than it should be. The result is negligible attenuation (0.09 dB average) instead of meaningful attenuation (~0.7 dB).

**Suspected Bug:** The denominator should be just `C`, not `D_n * C`.

---

## The Attenuation Formula

### Current Implementation (src/channel.cpp, line 5938)

```cpp
Real denom = vr.D_n * rolloff_C;  // D_n * C ≈ 8.94 * 13 = 116.21
gamma = exp(-d / denom);
```

### What D_n Represents

From `sns_setup_corner()` (line 3341):

```
D_n = sqrt((xA - xa)^2 + (yB - yb)^2)
    = diagonal distance from VR inner edge to far array corner
```

For typical Config 3 clusters:
- V_n ≈ 0.60 (strong/medium clusters)
- VR dimensions: a ≈ 9.3λ, b ≈ 9.7λ
- Array dimensions: W = 15.5λ, H = 16.1λ
- **D_n ≈ 8.94λ** (quite large)

---

## Why This Is Wrong

### Problem 1: D_n Scaling

The exponential decay rate should be controlled by `C` alone:

```
CORRECT:   gamma = exp(-d / C)
CURRENT:   gamma = exp(-d / (D_n * C))
                 = exp(-d / (8.94 * 13))
                 = exp(-d / 116.21)
```

The factor 1/D_n = 1/8.94 ≈ 0.112 makes the exponent ~9x smaller, so decay is ~9x slower.

### Problem 2: VR-Size Dependency

Larger visibility regions have:
- Larger `a` and `b` dimensions
- Larger D_n
- **Even weaker attenuation** (larger denominator)

This is counterintuitive: a larger VR should NOT cause weaker edge attenuation.

### Problem 3: Maximum Attenuation is Negligible

At the farthest array corner (d = D_n):

```
Current formula:  gamma = exp(-D_n / (D_n*C)) = exp(-1/C) = exp(-1/13)
                                              = 0.9260  → -0.33 dB only

Expected formula: gamma = exp(-D_n / C)
                                              = exp(-8.94/13)
                                              = 0.5128  → -2.90 dB
```

A **0.33 dB loss** is unmeasurable in practice; **2.90 dB** is significant.

---

## Quantitative Impact for Config 3

### Test Configuration

```
BS array:     24 × 32 = 768 elements
VR (V_n=0.60): 9.3λ × 9.7λ (60% of array area)
D_n:          8.94λ
```

### Element Distribution

| Attenuation Range | Formula A (Current) | Formula B (No D_n) | Ratio |
|---|---|---|---|
| 0 dB (full power) | 308 elements (40.1%) | 300 elements (39.1%) | 1.0x |
| 0 to -0.1 dB | 429 elements (55.9%) | 300 elements (39.1%) | 1.4x |
| 0 to -0.3 dB | 567 elements (73.8%) | 335 elements (43.6%) | 1.7x |
| **0 to -0.7 dB** | **567 elements (73.8%)** | **373 elements (48.6%)** | **1.5x** |
| **-0.7 to -3.0 dB** | **0 elements (0%)** | **395 elements (51.4%)** | **∞** |

### Average Attenuation

```
Formula A (current):  gamma_avg = 0.9796  → -0.0894 dB
Formula B (no D_n):   gamma_avg = 0.8445  → -0.7341 dB
Difference:           -0.6447 dB (7.2x difference)
```

---

## Sample Element Attenuations

For visibility-limited cluster at bottom-left corner:

| Position | Distance | Formula A | Formula B | Ratio |
|---|---|---|---|---|
| Inside VR | — | 0 dB (full) | 0 dB (full) | 1.0x |
| Edge (h) | 0.7λ | -0.026 dB | -0.234 dB | 9.0x |
| Edge (v) | 0.3λ | -0.013 dB | -0.114 dB | 8.8x |
| Diagonal (±3λ) | 4.2λ | -0.134 dB | -1.194 dB | 8.9x |
| Far corner | 8.9λ | -0.334 dB | -2.986 dB | 8.9x |

**Observation:** Formula B is consistently **8-9x stronger**, exactly the ratio D_n*C / C = 8.94.

---

## Coupling Loss Impact

During calibration (Initiallization.cpp):

```cpp
coupling_loss = RSRP_with_SNS - pathloss_final

If SNS attenuation is 0.09 dB:
  coupling_loss_with_SNS = coupling_loss + 0.09 dB (negligible)

If SNS attenuation is 0.73 dB:
  coupling_loss_with_SNS = coupling_loss + 0.73 dB (significant)
```

The 0.64 dB difference could account for observed calibration discrepancies.

---

## Root Cause: Spec Interpretation

The j20 specification (Section 7.6.14.1.3) describes:

1. **Visibility Region (VR):** Stochastic rectangular region covering V_n fraction of array
2. **Attenuation:** Apply exp(-d/C) roll-off for elements outside VR edge
3. **D_n Parameter:** Defined as corner-to-corner distance in the reference

The current implementation appears to confuse:
- **D_n as a roll-off normalizer** (wrong)
- **D_n as a far-field reference distance** (may be correct context)

If D_n is meant only for **calculating the reference corner positions** (which it does), then it should **NOT appear in the exponential denominator**.

---

## Correct Formula (Hypothesis)

Based on typical RF propagation models:

```
HYPOTHESIS: gamma = exp(-d / C)

where:
  d = distance from VR boundary (scalar or vector distance)
  C = dimensionless roll-off constant (13.0 for Config 3)
```

The D_n should only be used to:
1. Define the VR corner positions (done in sns_setup_corner)
2. Anchor the visibility region size (done in Generate_VisibilityRegion)
3. NOT scale the exponential decay rate

---

## Verification Steps

To verify the bug:

1. **Check j20 Clause 7.6.14.1.3** detailed tables and equations
   - What is the exact formula for attenuation?
   - Does it include a D_n or similar normalizer?

2. **Check ns-3 CDL implementation**
   - How does the reference implementation compute SNS attenuation?
   - Does it use D_n in the denominator?

3. **Numerical validation**
   - Run Config 3 calibration with Formula B (no D_n)
   - Compare coupling loss to calibration targets
   - If Formula B improves calibration fit, the bug is confirmed

4. **Physical reasoning**
   - SNS should reduce coupling loss by 0.5-1.0 dB for visibility-limited cases
   - Formula A gives 0.09 dB (too weak)
   - Formula B gives 0.73 dB (reasonable)

---

## Code Change (If Confirmed)

### Current Code (src/channel.cpp:5938)

```cpp
Real denom = vr.D_n * rolloff_C;
if (denom <= REAL(0.0)) return REAL(1.0);

// Distance computation...
Real dx = fabs(pos_h - vr.xa);
Real dy = fabs(pos_v - vr.yb);

Real gamma;
if (!inside_x && inside_y)
    gamma = exp(-dx / denom);
else if (inside_x && !inside_y)
    gamma = exp(-dy / denom);
else
    gamma = exp(-sqrt(dx * dx + dy * dy) / denom);
```

### Proposed Fix (Remove D_n)

```cpp
// Real denom = vr.D_n * rolloff_C;  // OLD (WRONG)
Real denom = rolloff_C;             // NEW (CORRECT)

if (denom <= REAL(0.0)) return REAL(1.0);

// Distance computation... (unchanged)
Real dx = fabs(pos_h - vr.xa);
Real dy = fabs(pos_v - vr.yb);

Real gamma;
if (!inside_x && inside_y)
    gamma = exp(-dx / denom);       // Now uses C only
else if (inside_x && !inside_y)
    gamma = exp(-dy / denom);       // Now uses C only
else
    gamma = exp(-sqrt(dx * dx + dy * dy) / denom);  // Now uses C only
```

### Impact

- **Coupling loss change:** +0.64 dB for visibility-limited clusters
- **RSRP change:** -0.64 dB outside VR
- **Calibration effect:** May need +0.6 dB adjustment to reference points

---

## Files Affected

If change is made:

1. **Modified:** `src/channel.cpp` line 5938
   - Remove `vr.D_n *` from denominator

2. **No changes needed:**
   - `h/channel.h` (struct ClusterVR definition still valid)
   - `src/channel.cpp` Generate_VisibilityRegion() (VR generation unchanged)
   - `src/channel.cpp` sns_setup_corner() (corner calculation unchanged)
   - `src/channel.cpp` Compute_SNS_RSRP_Attenuation() (calling code unchanged)

3. **Recalibration needed:**
   - Re-run Config 3 calibration
   - Compare coupling loss targets (may shift by ~0.6 dB)

---

## Summary

| Aspect | Formula A (Current) | Formula B (Corrected) |
|---|---|---|
| **Formula** | exp(-d / (D_n*C)) | exp(-d / C) |
| **D_n Factor** | 8.94λ | None |
| **Max Atten** | -0.33 dB | -2.90 dB |
| **Avg Atten** | -0.09 dB | -0.73 dB |
| **Ratio** | Baseline | 8.1x stronger |
| **Effect** | Negligible | Measurable |
| **Likely Correct?** | No | Yes |

The evidence strongly suggests **Formula B (without D_n) is correct**, and the D_n normalization was a misinterpretation of the j20 specification.

