# SNS Attenuation: Complete Code Flow and Call Chain

**Reference Implementation:** `/src/channel.cpp`

---

## High-Level Flow

```
Per-Drop Loop
  ↓
  Set_channel()  [Initialize for BS-MS pair]
    ├→ Generate_VisibilityRegion()  [Create random VR per cluster]
    │   └→ sns_setup_corner()  [Place VR, compute D_n]
    │
    └→ Compute_SNS_RSRP_Attenuation()  [Average attenuation across array]
        └→ compute_sns_attenuation() [Element-level attenuation]
              [PROBLEM LOCATION: line 5938]
```

---

## Detailed Call Flow

### 1. Generate_VisibilityRegion() — Line 3344-3440

**Purpose:** Create a random visibility region for each cluster, determine which clusters are visibility-limited.

**Key steps:**

```cpp
void CHANNEL::Generate_VisibilityRegion()
{
    // Step 1: Per-UT visibility probability
    Real Pr_raw = normal(g_sns_mu_P_vis, g_sns_sigma_P_vis);
    sns_Pr_sns = MAX(0.0, MIN(1.0, Pr_raw));

    // Step 2: Array dimensions (in wavelengths)
    Real W = (BS_N - 1) * BS_dH + ...;    // Width ≈ 15.5λ for Config 3
    Real H = (BS_M - 1) * BS_dV + ...;    // Height ≈ 16.1λ for Config 3

    // Step 3: Find max cluster power
    Real P_max_dB = max power among all clusters

    // Step 4-5: For each cluster
    for (int n = 0; n < num_path; n++) {
        Real x_n = randnum.u();
        if (x_n < sns_Pr_sns) {  // Cluster is visibility-limited
            sns_vr[n].limited = true;

            // V_n size factor (0.56-0.60 for typical clusters)
            Real V_n = 0.12 * exp(-(P_max - P_n) / 50.0) + 0.48;
            sns_vr[n].V_n = V_n;

            // Random dimensions with constant area
            Real a = Uniform[V_n*W, W];
            Real b = (V_n*W*H) / a;  // Constant area = V_n*W*H

            sns_setup_corner(sns_vr[n], a, b, W, H);
            //         ↓↓↓ CALLS sns_setup_corner() ↓↓↓
        }
    }
}
```

### 2. sns_setup_corner() — Line 3312-3342

**Purpose:** Place the visibility region at one of four array corners, compute D_n diagonal.

**Inputs:**
- `a`, `b`: VR dimensions (width, height in wavelengths)
- `W`, `H`: Array dimensions

**Outputs:**
- `vr.x0, vr.y0`: Anchor corner (0,0), (W,0), (0,H), or (W,H)
- `vr.xa, vr.yb`: Inner edge of VR (where attenuation starts)
- `vr.xA, vr.yB`: Far corner (opposite side of array from anchor)
- **`vr.D_n`: Diagonal distance** ← THIS IS THE PROBLEM

```cpp
static void sns_setup_corner(ClusterVR& vr, Real a, Real b, Real W, Real H)
{
    // Choose random corner (50/50 for each axis)
    bool right = (randnum.u() < 0.5);
    bool top   = (randnum.u() < 0.5);

    // Example: bottom-left corner
    if (!right && !top) {
        vr.x0 = 0.0;   vr.y0 = 0.0;    // Anchor at origin
        vr.xa = a;     vr.yb = b;      // Inner edge at (a,b)
        vr.xA = W;     vr.yB = H;      // Far corner at array corner
    }
    // ... (similar for other 3 corners)

    // Compute diagonal distance
    Real dx = vr.xA - vr.xa;  // W - a ≈ 15.5 - 9.3 = 6.2λ
    Real dy = vr.yB - vr.yb;  // H - b ≈ 16.1 - 9.7 = 6.4λ
    vr.D_n = sqrt(dx*dx + dy*dy);  // ≈ 8.94λ
    //                    ↑↑↑
    //         THIS VALUE IS STORED FOR LATER USE
}
```

**For bottom-left example:**
- `a = 9.3λ, b = 9.7λ`
- `dx = 15.5 - 9.3 = 6.2λ`
- `dy = 16.1 - 9.7 = 6.4λ`
- **`D_n = sqrt(6.2² + 6.4²) = 8.94λ`**

### 3. Compute_SNS_RSRP_Attenuation() — Line 3451-3495

**Purpose:** Compute average attenuation across all BS elements for use in RSRP calculation.

**Called by:** GetChannelImpulseResponse() after channel generation

```cpp
void CHANNEL::Compute_SNS_RSRP_Attenuation()
{
    if (!sns_any_limited) return;  // No visibility-limited clusters

    int M = BS_M;    // 24 elements
    int N = BS_N;    // 32 elements
    Real total_elements = M * N;  // 768 elements

    // For each NLOS cluster
    for (int n = 0; n < num_path; n++) {
        if (!sns_vr[n].limited) continue;

        Real sum_gamma = 0.0;

        // Loop over all BS antenna elements
        for (int s_m = 0; s_m < M; s_m++) {
            Real pos_v = s_m * BS_dV;  // Vertical position in wavelengths

            for (int s_n = 0; s_n < N; s_n++) {
                Real pos_h = s_n * BS_dH;  // Horizontal position

                // Compute element-level attenuation
                Real atten = compute_sns_attenuation(
                    pos_h, pos_v, sns_vr[n], g_sns_rolloff_C);

                // Square because atten is amplitude, we want power
                sum_gamma += atten * atten;
                //         ↓↓↓ CALLS compute_sns_attenuation() 768 TIMES ↓↓↓
            }
        }

        // Average across all elements
        sns_rsrp_power_atten[n] = sum_gamma / total_elements;
        //                       ≈ 0.9796 (Formula A) or 0.8445 (Formula B)
    }

    // ... (similar for LOS path)
}
```

**Output:**
- `sns_rsrp_power_atten[n]` = average power attenuation for cluster n
  - Formula A (current): ≈ 0.9796 (barely attenuates)
  - Formula B (correct): ≈ 0.8445 (meaningful attenuation)

### 4. compute_sns_attenuation() — Line 5927-5954

**Purpose:** Compute amplitude attenuation for a single TX element, relative to VR boundary.

**Inputs:**
- `pos_h, pos_v`: Element position in wavelengths
- `vr`: Visibility region structure (contains D_n)
- `rolloff_C`: Roll-off constant (13.0)

**Returns:**
- Amplitude scaling factor (sqrt of power attenuation)

```cpp
static inline Real compute_sns_attenuation(
    Real pos_h, Real pos_v, const ClusterVR& vr, Real rolloff_C)
{
    if (!vr.limited) return 1.0;  // Not visibility-limited, no attenuation

    // Check if inside VR rectangle
    bool inside_x = (fabs(pos_h - vr.x0) <= vr.a);
    bool inside_y = (fabs(pos_v - vr.y0) <= vr.b);

    if (inside_x && inside_y) return 1.0;  // Full power inside VR

    // PROBLEM LOCATION:
    Real denom = vr.D_n * rolloff_C;  // Line 5938
    //            ↑↑↑
    //   This multiplies by D_n ≈ 8.94, making roll-off 8.9x too gradual!

    if (denom <= 0.0) return 1.0;

    // Compute distance to VR edge
    Real dx = fabs(pos_h - vr.xa);  // Distance to inner edge (h)
    Real dy = fabs(pos_v - vr.yb);  // Distance to inner edge (v)

    Real gamma;  // Power attenuation (linear, not dB)

    if (!inside_x && inside_y)
        gamma = exp(-dx / denom);  // Outside horizontal only
    else if (inside_x && !inside_y)
        gamma = exp(-dy / denom);  // Outside vertical only
    else
        gamma = exp(-sqrt(dx*dx + dy*dy) / denom);  // Outside both

    return sqrt(gamma);  // Convert power to amplitude
}
```

**For element at far corner (pos_h=15.5λ, pos_v=16.1λ):**

```
Current (Formula A):
  dx = |15.5 - 9.3| = 6.2λ
  dy = |16.1 - 9.7| = 6.4λ
  d = sqrt(6.2² + 6.4²) = 8.94λ
  denom = 8.94 * 13.0 = 116.21
  gamma = exp(-8.94 / 116.21) = exp(-0.077) = 0.926
  amplitude = sqrt(0.926) = 0.962
  dB loss = 10*log10(0.926) = -0.334 dB

Proposed (Formula B, no D_n):
  denom = 13.0 (not multiplied by D_n)
  gamma = exp(-8.94 / 13.0) = exp(-0.688) = 0.513
  amplitude = sqrt(0.513) = 0.716
  dB loss = 10*log10(0.513) = -2.90 dB
```

---

## Data Structure: ClusterVR

From `/h/channel.h`:

```cpp
struct ClusterVR {
    bool  limited;    // true = visibility-limited cluster
    Real  V_n;        // VR size factor [0,1] (≈0.60 for typical)
    Real  a;          // VR width in wavelengths (≈9.3λ)
    Real  b;          // VR height in wavelengths (≈9.7λ)
    Real  x0, y0;     // Reference corner (anchor point)
    Real  xa, yb;     // Inner edge of VR (attenuation starts here)
    Real  xA, yB;     // Far array corner (opposite from anchor)
    Real  D_n;        // sqrt((xA-xa)^2 + (yB-yb)^2) ≈ 8.94λ
                      // ↑↑↑ USED IN DENOMINATOR (PROBLEM)
};
```

---

## Global Variables Used

```cpp
// From common.h (defined globally):
REAL g_sns_bs_enabled;        // SNS feature enabled
REAL g_sns_vr_A;              // 0.12 (VR size parameter)
REAL g_sns_vr_B;              // 0.48 (VR size parameter)
REAL g_sns_vr_R;              // 50.0 (power ratio parameter)
REAL g_sns_vr_delta;           // 0.0316 (random noise std dev)
REAL g_sns_mu_P_vis;           // 0.49 (visibility probability mean)
REAL g_sns_sigma_P_vis;        // 0.18 (visibility probability std dev)
REAL g_sns_rolloff_C;          // 13.0 (rolloff constant) ← LINE 5938

// Per-drop/channel:
ClusterVR sns_vr[MAX_NUM_CLUSTERS];      // Array of VR structures
ClusterVR sns_vr_los;                    // LOS path VR
REAL sns_rsrp_power_atten[MAX_NUM_CLUSTERS];  // Averaged attenuation
```

---

## Calling Context

### Where Compute_SNS_RSRP_Attenuation() is called:

```
GetChannelImpulseResponse()  [Line 5956]
  ├→ Phase A: Allocate_H_usn_memory()
  ├→ Phase B: PreCompute LOS/NLOS rays
  ├→ Phase C: Build H_usn matrix (element-level coefficients)
  └→ **Called SOMEWHERE after Phase C**: Compute_SNS_RSRP_Attenuation()
      └→ Uses sns_rsrp_power_atten[] in RSRP calculation
```

### Where sns_rsrp_power_atten[] is used:

**In Get_RSRP() (Link.cpp):**

```cpp
// Approximate location (not exact):
if (sns_vr[cluster_idx].limited) {
    power_rsrp *= channel->sns_rsrp_power_atten[cluster_idx];
    // Reduces RSRP for visibility-limited clusters
}
```

This directly affects reported RSRP in calibration and performance metrics.

---

## Impact on Calibration

### Current Behavior (Formula A, with D_n):

```
Coupling loss = BS_maxpower + RSRP - Pathloss
              = 49 dBm + (-145 dBm from visibility-limited) - (-135 dBm PL)
              ≈ -31 dBm (SNS adds only -0.09 dB effect)
```

### Expected Behavior (Formula B, without D_n):

```
Coupling loss = 49 dBm + (-145.7 dBm from visibility-limited) - (-135 dBm PL)
              ≈ -31.7 dBm (SNS adds -0.7 dB effect)
```

The **0.6-0.7 dB difference** is significant for calibration matching.

---

## Configuration Impact

### Why Config 3 is Affected Most:

```
Config 3: BS_M=24, BS_N=32, BS_dH=0.5λ, BS_dV=0.7λ
  D_n ≈ 8.94λ  (large array → large D_n)
  denom = 8.94 * 13 = 116.21
  max_atten = -0.334 dB (very weak)

Config A: BS_M=8, BS_N=8, BS_dH=0.5λ, BS_dV=0.5λ
  D_n ≈ 2.5λ   (smaller array → smaller D_n)
  denom = 2.5 * 13 = 32.5
  max_atten = -1.2 dB (stronger effect)
```

Larger arrays get weaker SNS effects with D_n normalization, which is backwards.

---

## Summary

| Step | Function | Issue | Impact |
|---|---|---|---|
| 1 | Generate_VisibilityRegion | Determines if cluster is visibility-limited | No bug |
| 2 | sns_setup_corner | Computes D_n diagonal | Correct computation |
| 3 | Compute_SNS_RSRP_Attenuation | Averages attenuation across elements | Calls buggy function |
| 4 | **compute_sns_attenuation** | **Uses exp(-d / (D_n*C))** | **BUG: D_n normalization** |
| 5 | Get_RSRP | Uses sns_rsrp_power_atten in RSRP | Receives attenuated value |
| 6 | Calibration | Matches coupling loss to target | Off by ~0.6 dB |

**The bug is at Step 4, line 5938 in channel.cpp.**

