# Code Locations Summary: ZSD Circular Spread Analysis

## Files Modified / Analyzed

### 1. ORIGINAL Implementation (Correct Behavior)

#### File: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/Link.cpp`

**Location:** Angle transformation function (CORRECT)
```cpp
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
    else if (x < -180)  // ← KEY: This branch handles angles < -180
    {
        x = x + 360;
    }
    return x;
}
```
**Status:** ✓ CORRECT - Handles all angle ranges properly

---

#### File: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/channel.cpp`

**Location 1:** Per-ray angle precomputation (lines 4253-4282)
```cpp
void CHANNEL::Precompute_ray_angles()
{
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
    {
        int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
        for (int j = 0; j < num_rays; j++)
        {
            // RX (Arrival) angles
            Real rx_azimuth = AOA[i] + cluster_ASA * offset_angle[i][j];
            Real rx_zenith  = ZOA[i] + cluster_ZSA * offset_angle[i][j];

            // TX (Departure) angles
            Real tx_azimuth = AOD[i] + cluster_ASD * offset_angle_rand_coupling[i][j];
            Real tx_zenith  = ZOD[i] + (3.0 / 8.0) * pow(10, mu_ZSD) * offset_angle_rand_coupling[i][j];

            // Transform to proper ranges
            rx_azimuth = Transform_angle_minus_180_to_plus_180(rx_azimuth);
            rx_zenith  = Transform_angle_0_to_plus_180(rx_zenith);  // ← USES CORRECT FUNCTION
            tx_azimuth = Transform_angle_minus_180_to_plus_180(tx_azimuth);
            tx_zenith  = Transform_angle_0_to_plus_180(tx_zenith);  // ← USES CORRECT FUNCTION

            // Store
            ray_AOA[i][j][0] = rx_azimuth;  // azimuth
            ray_AOA[i][j][1] = rx_zenith;   // zenith
            ray_AOD[i][j][0] = tx_azimuth;  // azimuth
            ray_AOD[i][j][1] = tx_zenith;   // zenith
        }
    }
    ray_angles_precomputed = true;
}
```
**Purpose:** Compute per-ray angles from cluster centers + offsets, applies transformation
**Status:** ✓ Uses correct transformation function

---

**Location 2:** ZOD generation (lines 3880-3967)
```cpp
// Lines 3888-3895: Initial ZOD based on power distribution
if (power[i] == 0.)
{
    ZOD[i] = 0.;
}
else
{
    ZOD[i] = -1 * ZSD * log(power[i] / strongest_power) / C;
}

// Lines 3925-3967: ZOD transformation based on propagation condition
if (Propagation == 1) // LOS
{
    ZOD[i] = (X * ZOD[i] + Y_ZOD) - (X1_ZOD * first_ZOD + Y1_ZOD - LOS_ZOD_GCS);
}
else // NLOS, O2I
{
    if (Propagation == 2)
    {
        ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
    }
    else if (Propagation == 0)
    {
        ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
    }
}
```
**Purpose:** Generate cluster-level ZOD values (can range beyond [0, 180])
**Note:** These values are later combined with per-ray offsets and transformed
**Status:** ✓ Correct - results properly handled by Transform_angle_0_to_plus_180()

---

**Location 3:** Circular angle spread computation (lines 4412-4454)
```cpp
void CHANNEL::Set_circular_angle_spread()
{
    // Count total rays across all clusters/subclusters
    int total_rays = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
        total_rays += (int)NUM_RAY_per_ClusterNUM[i];

    if (total_rays <= 0) return;

    // Flatten per-ray angles and powers into temporary arrays
    Real* flat_aoa = new Real[total_rays];
    Real* flat_aod = new Real[total_rays];
    Real* flat_zoa = new Real[total_rays];
    Real* flat_zod = new Real[total_rays];
    Real* flat_pow = new Real[total_rays];

    int idx = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
    {
        int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
        Real ray_power = power[i] / num_rays;
        for (int j = 0; j < num_rays; j++)
        {
            flat_aoa[idx] = ray_AOA[i][j][0];  // azimuth
            flat_zoa[idx] = ray_AOA[i][j][1];  // zenith (TRANSFORMED)
            flat_aod[idx] = ray_AOD[i][j][0];  // azimuth
            flat_zod[idx] = ray_AOD[i][j][1];  // zenith (TRANSFORMED) ← After Transform_angle_0_to_plus_180()
            flat_pow[idx] = ray_power;
            idx++;
        }
    }

    circular_angle_spread_AOA = compute_circular_angle_spread(flat_aoa, flat_pow, total_rays);
    circular_angle_spread_AOD = compute_circular_angle_spread(flat_aod, flat_pow, total_rays);
    circular_angle_spread_ZOA = compute_circular_angle_spread(flat_zoa, flat_pow, total_rays);
    circular_angle_spread_ZOD = compute_circular_angle_spread(flat_zod, flat_pow, total_rays);

    delete[] flat_aoa;
    delete[] flat_aod;
    delete[] flat_zoa;
    delete[] flat_zod;
    delete[] flat_pow;
}
```
**Purpose:** Flatten per-ray angles and call circular spread helper
**Status:** ✓ Correct - receives properly transformed angles from Precompute_ray_angles()

---

**Location 4:** Circular spread helper (lines 4372-4410)
```cpp
static Real compute_circular_angle_spread(Real* angles, Real* powers, int N)
{
    Real power_sum = 0.0;
    for (int i = 0; i < N; i++) power_sum += powers[i];
    if (power_sum <= 0.0) return 0.0;

    Real min_spread = 1e30;
    for (int delta = 0; delta < 360; delta++)
    {
        Real d = (Real)delta;
        // Eq. A-6: weighted mean of shifted angles
        Real mu = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real shifted = fmod(angles[i] + d + 180.0, 360.0);
            if (shifted < 0.0) shifted += 360.0;
            shifted -= 180.0;
            mu += shifted * powers[i];
        }
        mu /= power_sum;

        // Eq. A-4,5: spread
        Real variance = 0.0;
        for (int i = 0; i < N; i++)
        {
            Real shifted = fmod(angles[i] + d + 180.0, 360.0);
            if (shifted < 0.0) shifted += 360.0;
            shifted -= 180.0;
            Real dev = fmod(shifted - mu + 180.0, 360.0);
            if (dev < 0.0) dev += 360.0;
            dev -= 180.0;
            variance += dev * dev * powers[i];
        }
        Real spread = sqrt(variance / power_sum);
        if (spread < min_spread) min_spread = spread;
    }
    return min_spread;
}
```
**Purpose:** 3GPP TR 25.996 Annex A circular spread with delta optimization
**Status:** ✓ Correct algorithm

---

### 2. v2 Implementation (BUGGY)

#### File: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelModel.cpp`

**Location 1:** Angle transformation function (lines 459-467) — **BUGGY**
```cpp
Real TransformAngle0To180(Real x)
{
    x = fmod(x, (Real)360.0);
    if (x > 180.0)
        x = 360.0 - x;
    if (x < 0.0)                         // ← BUG: Should be "else if", not separate "if"
        x = -x;
    return x;
}
```
**Status:** ✗ BUGGY - Missing `else if (x < -180.0)` branch
**Impact:** Returns invalid values > 180 for inputs < -180

---

#### File: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelSSP.cpp`

**Location 1:** Per-ray angle precomputation (lines 642-659)
```cpp
void ChannelModel::PrecomputeRayAngles()
{
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++) {
        int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
        for (int j = 0; j < num_rays; j++) {
            Real rx_azimuth = AOA[i] + cluster_ASA * offset_angle[i][j];
            Real rx_zenith  = ZOA[i] + cluster_ZSA * offset_angle[i][j];
            Real tx_azimuth = AOD[i] + cluster_ASD * offset_angle_rand_coupling[i][j];
            Real tx_zenith  = ZOD[i] + (3.0 / 8.0) * pow(10.0, mu_ZSD) * offset_angle_rand_coupling[i][j];

            ray_AOA[i][j][0] = TransformAngleMinus180To180(rx_azimuth);
            ray_AOA[i][j][1] = TransformAngle0To180(rx_zenith);      // ← USES BUGGY FUNCTION
            ray_AOD[i][j][0] = TransformAngleMinus180To180(tx_azimuth);
            ray_AOD[i][j][1] = TransformAngle0To180(tx_zenith);      // ← USES BUGGY FUNCTION
        }
    }
    ray_angles_precomputed = true;
}
```
**Purpose:** Same as original - compute per-ray angles
**Status:** ✗ Uses buggy transformation function, produces invalid angles

---

**Location 2:** ZOD generation (lines 485-506 in ChannelSSP.cpp)
```cpp
if (power[i] == 0.0) {
    ZOD[i] = 0.0;
    ZOA[i] = 0.0;
} else {
    ZOD[i] = -1.0 * ZSD * log(power[i] / strongest_power) / C;
    ZOA[i] = -1.0 * ZSA * log(power[i] / strongest_power) / C;
}

// ... (lines 494-506)
if (Propagation == LOS_propagation) {
    ZOA[i] = (X * ZOA[i] + Y_ZOA) - (X1_ZOA * first_ZOA + Y1_ZOA - LOS_ZOA_GCS);
    ZOD[i] = (X * ZOD[i] + Y_ZOD) - (X1_ZOD * first_ZOD + Y1_ZOD - LOS_ZOD_GCS);
} else if (Propagation == OUT2IN_propagation) {
    ZOA[i] = X * ZOA[i] + Y_ZOA + 90.0;
    ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
} else {
    ZOA[i] = X * ZOA[i] + Y_ZOA + LOS_ZOA_GCS;
    ZOD[i] = X * ZOD[i] + Y_ZOD + LOS_ZOD_GCS + mu_offset_ZOD;
}
```
**Purpose:** Same as original - generate cluster-level ZOD
**Status:** ✓ Correct generation logic, but results mishandled by buggy transformation

---

**Location 3:** Circular angle spread computation (lines 664-694 in ChannelSSP.cpp)
```cpp
void ChannelModel::ComputeCircularAngleSpread()
{
    int total_rays = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
        total_rays += (int)NUM_RAY_per_ClusterNUM[i];
    if (total_rays <= 0) return;

    Real* flat_aoa = new Real[total_rays];
    Real* flat_aod = new Real[total_rays];
    Real* flat_zoa = new Real[total_rays];
    Real* flat_zod = new Real[total_rays];
    Real* flat_pow = new Real[total_rays];

    int idx = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++) {
        int nr = (int)NUM_RAY_per_ClusterNUM[i];
        Real rp = power[i] / nr;
        for (int j = 0; j < nr; j++) {
            flat_aoa[idx] = ray_AOA[i][j][0];
            flat_zoa[idx] = ray_AOA[i][j][1];
            flat_aod[idx] = ray_AOD[i][j][0];
            flat_zod[idx] = ray_AOD[i][j][1];  // ← Gets CORRUPTED values from buggy transform
            flat_pow[idx] = rp;
            idx++;
        }
    }

    circular_angle_spread_AOA = ComputeCircularAngleSpreadHelper(flat_aoa, flat_pow, total_rays);
    circular_angle_spread_AOD = ComputeCircularAngleSpreadHelper(flat_aod, flat_pow, total_rays);
    circular_angle_spread_ZOA = ComputeCircularAngleSpreadHelper(flat_zoa, flat_pow, total_rays);
    circular_angle_spread_ZOD = ComputeCircularAngleSpreadHelper(flat_zod, flat_pow, total_rays);

    delete[] flat_aoa;
    delete[] flat_aod;
    delete[] flat_zoa;
    delete[] flat_zod;
    delete[] flat_pow;
}
```
**Purpose:** Same as original
**Status:** ✗ Receives corrupted angles from buggy PrecomputeRayAngles()

---

**Location 4:** Circular spread helper (lines 487-519 in ChannelModel.cpp)
```cpp
Real ComputeCircularAngleSpreadHelper(Real* angles, Real* powers, int N)
{
    Real power_sum = 0.0;
    for (int i = 0; i < N; i++) power_sum += powers[i];
    if (power_sum <= 0.0) return 0.0;

    Real min_spread = 1e30;
    for (int delta = 0; delta < 360; delta++) {
        Real d = (Real)delta;
        Real mu = 0.0;
        for (int i = 0; i < N; i++) {
            Real shifted = fmod(angles[i] + d + 180.0, 360.0);
            if (shifted < 0.0) shifted += 360.0;
            shifted -= 180.0;
            mu += shifted * powers[i];
        }
        mu /= power_sum;

        Real variance = 0.0;
        for (int i = 0; i < N; i++) {
            Real shifted = fmod(angles[i] + d + 180.0, 360.0);
            if (shifted < 0.0) shifted += 360.0;
            shifted -= 180.0;
            Real dev = fmod(shifted - mu + 180.0, 360.0);
            if (dev < 0.0) dev += 360.0;
            dev -= 180.0;
            variance += dev * dev * powers[i];
        }
        Real spread = sqrt(variance / power_sum);
        if (spread < min_spread) min_spread = spread;
    }
    return min_spread;
}
```
**Purpose:** Same as original
**Status:** ✓ Correct algorithm, but receives corrupted input data

---

## Data Flow Comparison

### ORIGINAL (Correct)
```
ZOD[] (cluster-level)
    ↓
Precompute_ray_angles():
  tx_zenith = ZOD[i] + offset
  ray_AOD[i][j][1] = Transform_angle_0_to_plus_180(tx_zenith)  ✓ CORRECT TRANSFORM
    ↓
Set_circular_angle_spread():
  flat_zod[idx] = ray_AOD[i][j][1]  (values in [0, 180])
    ↓
compute_circular_angle_spread(flat_zod)
  → Result: circular_angle_spread_ZOD ≈ 14.2° (at 95th percentile) ✓
```

### v2 (Buggy)
```
ZOD[] (cluster-level)
    ↓
PrecomputeRayAngles():
  tx_zenith = ZOD[i] + offset
  ray_AOD[i][j][1] = TransformAngle0To180(tx_zenith)  ✗ BUGGY TRANSFORM
                     (produces values > 180 when input < -180)
    ↓
ComputeCircularAngleSpread():
  flat_zod[idx] = ray_AOD[i][j][1]  (CORRUPTED: some values > 180!)
    ↓
ComputeCircularAngleSpreadHelper(flat_zod)
  → Result: circular_angle_spread_ZOD ≈ 1.4° (at 95th percentile) ✗
```

---

## Summary Table

| Aspect | ORIGINAL | v2 | Status |
|--------|----------|-----|--------|
| Transform function logic | if-else if-else if chain | Two separate if statements | ✗ v2 BUGGY |
| Handles (180, 360) | Yes | Yes | ✓ Both OK |
| Handles [-180, 0) | Yes | Yes | ✓ Both OK |
| Handles (-∞, -180) | Yes (line `else if (x < -180)`) | No (MISSING) | ✗ v2 BUGGY |
| Handles [0, 180] | Yes | Yes | ✓ Both OK |
| Output range guaranteed | [0, 180] | Can exceed 180 | ✗ v2 BUGGY |
| Ray angle precomputation | Identical logic | Identical logic | ✓ Both same |
| Circular spread helper | Identical code | Identical code | ✓ Both same |
| Circular spread flattening | Identical logic | Identical logic | ✓ Both same |
| ZSD statistics at 95% | 14.2° | 1.4° | ✗ v2 Wrong |

