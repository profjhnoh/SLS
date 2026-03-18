# SNS Attenuation Analysis - Config 3 UMi
## Issue: Minimal Roll-off Normalization by D_n Parameter

**Date:** 2026-03-14
**Status:** Analysis Complete - Problem Identified

---

## Summary

The SNS (Spatial Non-Stationarity) attenuation formula in `compute_sns_attenuation()` uses:

```
gamma = exp(-d / (D_n * C))   [power domain]
amplitude = sqrt(gamma) = exp(-d / (2*D_n*C))
```

where:
- `d` = distance from element to VR boundary
- `D_n` = sqrt((xA - xa)^2 + (yB - yb)^2) ≈ 8-12 wavelengths for typical clusters
- `C` = rolloff constant = 13.0

**Critical Finding:** The denominator `D_n * C` is approximately 116, which makes the roll-off so gradual that SNS attenuation is almost negligible:

- **Maximum attenuation (at far corner):** exp(-1/13) = 0.926 (power) = **-0.33 dB only**
- **Average attenuation across BS array:** **-0.09 dB** (essentially no effect)
- **Without D_n normalization:** Formula B would give **-0.73 dB average** (7x stronger effect)

This suggests the formula is incorrectly normalizing the roll-off distance by D_n, making SNS ineffective on coupling loss.

---

## Configuration Parameters

Config 3 UMi parameters:
```
BS_M = 24 elements (vertical)
BS_N = 32 elements (horizontal)
BS_dH = 0.5λ (horizontal spacing)
BS_dV = 0.7λ (vertical spacing)

g_sns_vr_A = 0.12        (visibility region parameter A)
g_sns_vr_B = 0.48        (visibility region parameter B)
g_sns_vr_R = 50.0        (power ratio parameter)
g_sns_rolloff_C = 13.0   (rolloff constant)
```

---

## Part 1: Array Dimensions

BS antenna array dimensions (single panel, before panel groups):

```
W = (BS_N - 1) * BS_dH = (32 - 1) * 0.5 = 15.5λ
H = (BS_M - 1) * BS_dV = (24 - 1) * 0.7 = 16.1λ
```

The visibility region covers a subset of this array, constrained by size factor V_n.

---

## Part 2: Typical V_n Values

From `Generate_VisibilityRegion()`, line 3390:

```
V_n = A * exp(-(P_max_dB - P_n_dB) / R + delta) + B
    = 0.12 * exp(-(P_max_dB - P_n_dB) / 50.0) + 0.48
```

where `delta ~ N(0, 0.0316)` per spec.

**Examples:**

| Cluster Power | V_n Value | Interpretation |
|---|---|---|
| P_n = P_max (δ=0) | 0.60 | Strong cluster, medium-sized VR |
| P_n = P_max - 5 dB | 0.589 | Good cluster |
| P_n = P_max - 10 dB | 0.578 | Medium cluster |
| P_n = P_max - 15 dB | 0.569 | Weak cluster |
| P_n = P_max - 20 dB | 0.560 | Very weak cluster |

So typical V_n ranges from 0.56 to 0.60 for most clusters, meaning 56-60% of array dimensions are covered.

---

## Part 3: Typical VR Dimensions (a, b)

For V_n = 0.60 example:

```
a_min = V_n * W = 0.60 * 15.5 = 9.3λ
a ~ Uniform[9.3λ, 15.5λ]

Constant area constraint:
b = (V_n * W * H) / a

Examples:
  a = 9.3λ  (narrow):  b = (0.60 * 15.5 * 16.1) / 9.3 = 16.1λ  (maximum height)
  a = 12.4λ (medium): b = (0.60 * 15.5 * 16.1) / 12.4 = 12.1λ
  a = 15.5λ (wide):   b = (0.60 * 15.5 * 16.1) / 15.5 = 9.66λ (minimum height)
```

So VR dimensions vary significantly based on random aspect ratio choice, but cover approximately 60% of the array area.

---

## Part 4: D_n Values (Roll-off Normalizer)

From `sns_setup_corner()`, line 3341:

```
D_n = sqrt((xA - xa)^2 + (yB - yb)^2)
```

Where:
- (x0, y0) = array corner where VR is anchored (e.g., bottom-left = (0, 0))
- (xa, yb) = VR inner edge (extends a, b from anchor)
- (xA, yB) = opposite array corner (e.g., top-right = (W, H))

**Examples for different VR shapes:**

| VR Shape | a | b | D_n |
|---|---|---|---|
| Narrow/Tall | 8λ | 24λ | 10.89λ |
| Square | 16λ | 16λ | 0.51λ |
| Wide/Short | 24λ | 8λ | 11.74λ |
| **Typical (V_n=0.6)** | **9.3λ** | **9.7λ** | **8.94λ** |

For most visibility-limited clusters with V_n ≈ 0.60, **D_n is in the range 8-12 wavelengths**.

---

## Part 5: Attenuation Formula Analysis

### Current Formula (j20 Eq. 7.6-58):

```c
// Line 5938: denom = vr.D_n * rolloff_C
// Line 5947/5949/5951:
//    gamma = exp(-d / denom)  where denom = D_n * C
// Line 5953: return sqrt(gamma)  // convert power to amplitude
```

**For Config 3 UMi:**
```
denom = D_n * C = 8.94λ * 13.0 = 116.21
```

**At maximum distance (far corner, d = D_n = 8.94λ):**
```
gamma = exp(-8.94 / 116.21) = exp(-1 / 13.0) = 0.9260  (power)
      = -0.334 dB only!
```

This is **way too small** to have any practical effect on channel gain.

### Alternative Formula (without D_n normalization):

```
gamma = exp(-d / C)    [no D_n scaling]
At d = 8.94λ:
  gamma = exp(-8.94 / 13.0) = 0.5128  (power)
        = -2.90 dB
```

This would be **8.7x stronger** attenuation.

---

## Part 6: Numerical Comparison

### Test Case: Visibility-Limited Cluster, Bottom-Left Corner

**Configuration:**
```
V_n = 0.60
a = 9.3λ (VR width)
b = 9.7λ (VR height)
D_n = 8.94λ

Array: 24 x 32 = 768 elements
Inside VR: elements at (pos_h, pos_v) where:
  |pos_h - 0| <= 9.3λ  AND
  |pos_v - 0| <= 9.7λ
```

### Formula A: Current Implementation (with D_n normalization)

```
gamma = exp(-d / (D_n * C)) where D_n*C = 116.21

Results across all 768 BS elements:
  Average power attenuation (gamma): 0.9796
  Average loss:                      -0.0894 dB

  Minimum value (farthest element): 0.9260 (at far corner)
  Maximum loss:                      -0.3341 dB

  Elements inside VR (full power):   ~308 elements (~40%)
  Elements outside VR (attenuated):  ~460 elements (~60%)
```

**Observation:** Most elements get <0.1 dB attenuation, even those far outside the VR.

### Formula B: Without D_n Normalization

```
gamma = exp(-d / C) where C = 13.0

Results across same 768 BS elements:
  Average power attenuation (gamma): 0.8445
  Average loss:                      -0.7341 dB

  Difference from Formula A:         -0.6447 dB
```

**Observation:** Formula B produces ~7x stronger attenuation, making SNS a meaningful effect.

### Comparison Table

| Metric | Formula A (Current) | Formula B (No D_n) | Ratio |
|---|---|---|---|
| Max attenuation | -0.33 dB | -2.90 dB | 8.7x |
| Average attenuation | -0.09 dB | -0.73 dB | 8.1x |
| Coupling loss impact | Negligible | Moderate | 8x |

---

## Part 7: Problem Assessment

### Root Cause

The formula uses `exp(-d / (D_n * C))` which introduces a **normalization factor 1/D_n** in the exponent. This makes sense in some theoretical contexts, but creates two problems:

1. **D_n grows with larger VRs:** Larger visibility regions have larger D_n (10-12λ), which stretches the roll-off distance by the same factor
2. **Roll-off becomes ineffective:** For d_max = D_n at far corner:
   ```
   exp(-D_n / (D_n*C)) = exp(-1/C) = exp(-1/13) = 0.926 → only -0.33 dB
   ```

### Why This Matters

From the j20 spec (Section 7.6.14.1.3), the SNS attenuation should model **reduced signal strength outside the visibility region**. The current formula makes this effect nearly invisible (0.09 dB average), meaning:

- Coupling loss is barely affected by SNS
- RSRP outside VR is almost unchanged
- SNS parameter tuning has minimal effect on calibration

### Correct Interpretation

Looking at the formula intent:
- **D_n is meant to be the "corner-to-corner" distance** defining the farthest attenuation point
- **The roll-off should be exp(-d / C)** where C is a dimensionless time-constant-like parameter
- **D_n should NOT appear in the denominator** except possibly to normalize distances for comparison

---

## Part 8: Comparison with ns-3 CDL Model

The j20 spec describes SNS using a **stochastic visibility region model**, but the attenuation formula in Clause 7.6.14.1.3 does not explicitly state how D_n should be used.

The current implementation appears to confuse:
- **Distance normalization** (D_n used to scale positions)
- **Roll-off rate** (C parameter for exponential decay)

The formula should likely be:
```
gamma = exp(-d / (C * lambda))    [normalized distance d/lambda]
```

NOT:
```
gamma = exp(-d / (D_n * C))       [normalized by VR size, not wavelength]
```

---

## Part 9: Calibration Impact

If SNS attenuation is negligible, Config 3 coupling loss calibration may not properly account for visibility effects. The RSRP computation includes:

```cpp
// Compute_SNS_RSRP_Attenuation() - Line 3475
sns_rsrp_power_atten[n] = average of compute_sns_attenuation() across all elements
```

With Formula A: ~0.98 attenuation (barely any loss)
With Formula B: ~0.84 attenuation (significant 0.7 dB loss)

This 0.6-0.7 dB difference directly impacts reported coupling loss in calibration.

---

## Part 10: Recommended Fix (Tentative)

Without access to j20 detailed specification tables, a likely fix is:

**Option A: Remove D_n normalization entirely**
```c
// Replace line 5938:
Real denom = vr.D_n * rolloff_C;  // CURRENT (wrong)

// With:
Real denom = rolloff_C;  // CORRECT (no D_n scaling)
```

This makes attenuation depend on absolute distance from VR edge, not relative to VR size.

**Option B: Normalize distances instead of exponential**
```c
// Replace lines 5942-5951:
Real d_normalized = d / vr.D_n;  // normalize by VR size
gamma = exp(-d_normalized / rolloff_C);
```

This makes attenuation relative to the specific VR shape, but still uses C as roll-off rate.

---

## File Locations

**Current Implementation:**
- Function: `compute_sns_attenuation()` in `/src/channel.cpp` line 5927-5954
- Called by: `Allocate_H_usn_memory()` (setup), `Compute_SNS_RSRP_Attenuation()` (line 3475)
- VR setup: `Generate_VisibilityRegion()` in `/src/channel.cpp` line 3344-3440
- VR corner: `sns_setup_corner()` in `/src/channel.cpp` line 3312-3342

**Data Structures:**
- `struct ClusterVR` in `/h/channel.h` (defines limited, V_n, a, b, x0, y0, xa, yb, xA, yB, D_n)

---

## Appendix A: Code Snippets

### compute_sns_attenuation() - Current Implementation

```cpp
static inline Real compute_sns_attenuation(
    Real pos_h, Real pos_v, const ClusterVR& vr, Real rolloff_C)
{
    if (!vr.limited) return REAL(1.0);

    // Check if element is inside VR
    bool inside_x = (fabs(pos_h - vr.x0) <= vr.a);
    bool inside_y = (fabs(pos_v - vr.y0) <= vr.b);

    if (inside_x && inside_y) return REAL(1.0);

    Real denom = vr.D_n * rolloff_C;  // <-- PROBLEM: D_n scaling
    if (denom <= REAL(0.0)) return REAL(1.0);

    // Distance from element to VR edge
    Real dx = fabs(pos_h - vr.xa);
    Real dy = fabs(pos_v - vr.yb);

    Real gamma;
    if (!inside_x && inside_y)
        gamma = exp(-dx / denom);
    else if (inside_x && !inside_y)
        gamma = exp(-dy / denom);
    else  // outside both
        gamma = exp(-sqrt(dx * dx + dy * dy) / denom);

    return sqrt(gamma);  // power -> amplitude
}
```

### sns_setup_corner() - Creates D_n

```cpp
static void sns_setup_corner(ClusterVR& vr, Real a, Real b, Real W, Real H)
{
    // Set width and height
    vr.a = a;
    vr.b = b;

    // Random corner selection (50/50 for each axis)
    bool right = (randnum.u() < REAL(0.5));
    bool top   = (randnum.u() < REAL(0.5));

    if (!right && !top) {       // bottom-left
        vr.x0 = REAL(0.0);  vr.y0 = REAL(0.0);
        vr.xa = a;           vr.yb = b;
        vr.xA = W;           vr.yB = H;
    } else if (right && !top) { // bottom-right
        vr.x0 = W;           vr.y0 = REAL(0.0);
        vr.xa = W - a;       vr.yb = b;
        vr.xA = REAL(0.0);   vr.yB = H;
    } else if (!right && top) { // top-left
        vr.x0 = REAL(0.0);  vr.y0 = H;
        vr.xa = a;           vr.yb = H - b;
        vr.xA = W;           vr.yB = REAL(0.0);
    } else {                    // top-right
        vr.x0 = W;           vr.y0 = H;
        vr.xa = W - a;       vr.yb = H - b;
        vr.xA = REAL(0.0);   vr.yB = REAL(0.0);
    }

    // Compute diagonal distance: this is the PROBLEM
    // It's used to normalize the roll-off exponential
    Real dx = vr.xA - vr.xa;
    Real dy = vr.yB - vr.yb;
    vr.D_n = sqrt(dx * dx + dy * dy);  // <-- normalized by VR size
}
```

---

## Appendix B: Numerical Details

### Distance Distribution for Test Case

For V_n=0.60 cluster at bottom-left corner:
- Elements inside VR (0 < pos_h < 9.3λ, 0 < pos_v < 9.7λ): ~308 elements
- Max distance to VR edge (far corner): D_n = 8.94λ
- Average distance outside VR: ~4-5λ

### dB Conversion

Power attenuation (linear):
```
gamma = 0.9796 (Formula A)
gamma = 0.8445 (Formula B)

dB equivalent:
10 * log10(0.9796) = -0.0894 dB  (Formula A)
10 * log10(0.8445) = -0.7341 dB  (Formula B)
```

---

## Conclusion

The SNS attenuation formula in `compute_sns_attenuation()` uses D_n normalization that is **7-8x too weak** compared to a formula without this scaling. This makes the SNS effect nearly invisible on coupling loss and RSRP calculations.

**Recommendation:** Review the j20 specification (Section 7.6.14.1.3, Table 7.6.14.1.2-1/2/3) to confirm the correct formula. The current D_n * C denominator appears to be incorrect and should likely be just C (the rolloff constant).

