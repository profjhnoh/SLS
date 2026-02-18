# ZSD Circular Angle Spread Discrepancy: Root Cause Analysis

## Executive Summary

The 10x reduction in ZSD (Zenith Spread of Departure) circular angle spread values at the upper tail in v2 is **NOT** caused by differences in the circular angle spread algorithm itself (which is identical in both implementations). Instead, the root cause is a **subtle but critical difference in the angle transformation function** that converts zenith angles to the [0, 180] degree range.

**Critical Finding:** The v2 `TransformAngle0To180()` function has a bug that produces INCORRECT results for negative zenith angles in the range [-190, -180). This causes per-ray zenith angles near the boundaries to be systematically shifted AWAY from the full [0, 180] range, artificially compressing the spread.

---

## 1. Circular Angle Spread Algorithm Analysis

### Both Implementations Use Per-Ray Angles

Contrary to the initial hypothesis, **both implementations use per-ray angles** (flattened from all clusters/rays) to compute circular angle spread:

**ORIGINAL** (`src/channel.cpp:4414-4454`):
```cpp
void CHANNEL::Set_circular_angle_spread()
{
    int total_rays = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
        total_rays += (int)NUM_RAY_per_ClusterNUM[i];

    Real* flat_zod = new Real[total_rays];
    int idx = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++) {
        int num_rays = (int)NUM_RAY_per_ClusterNUM[i];
        Real ray_power = power[i] / num_rays;
        for (int j = 0; j < num_rays; j++) {
            flat_zod[idx] = ray_AOD[i][j][1];  // Per-ray ZOD
            flat_pow[idx] = ray_power;
            idx++;
        }
    }
    circular_angle_spread_ZOD = compute_circular_angle_spread(flat_zod, flat_pow, total_rays);
}
```

**v2** (`v2/src/ChannelSSP.cpp:664-694`):
```cpp
void ChannelModel::ComputeCircularAngleSpread()
{
    int total_rays = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++)
        total_rays += (int)NUM_RAY_per_ClusterNUM[i];

    Real* flat_zod = new Real[total_rays];
    int idx = 0;
    for (int i = 0; i < NUM_PATH_for_channelcoeff; i++) {
        int nr = (int)NUM_RAY_per_ClusterNUM[i];
        Real rp = power[i] / nr;
        for (int j = 0; j < nr; j++) {
            flat_zod[idx] = ray_AOD[i][j][1];  // Per-ray ZOD
            flat_pow[idx] = rp;
            idx++;
        }
    }
    circular_angle_spread_ZOD = ComputeCircularAngleSpreadHelper(flat_zod, flat_pow, total_rays);
}
```

**Conclusion:** Both implementations flatten per-ray angles identically. The algorithm difference hypothesis is **INCORRECT**.

---

## 2. Circular Angle Spread Helper Function Comparison

The `ComputeCircularAngleSpreadHelper()` function is **IDENTICAL** in both implementations (modulo whitespace):

**Both versions** (`src/channel.cpp:4373-4410` and `v2/src/ChannelModel.cpp:487-519`):
```cpp
Real compute_circular_angle_spread(Real* angles, Real* powers, int N)
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
        // ... variance calculation ...
        Real spread = sqrt(variance / power_sum);
        if (spread < min_spread) min_spread = spread;
    }
    return min_spread;
}
```

The algorithm itself does NOT have bugs for zenith angles. The delta optimization searches over all 360-degree rotations to find the optimal mean direction that minimizes variance, which is mathematically sound.

**Conclusion:** The circular spread helper function is **NOT the problem**.

---

## 3. CRITICAL BUG: Transform Angle Function Difference

This is where the root cause lies. The per-ray ZOD angles are computed before being passed to the circular spread function. They are transformed using `Transform_angle_0_to_plus_180()` / `TransformAngle0To180()`.

### ORIGINAL Implementation
**File:** `src/Link.cpp` (implementation) and `src/channel.cpp:4270` (usage)

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
    else if (x < -180)
    {
        x = x + 360;
    }
    return x;
}
```

### v2 Implementation
**File:** `v2/src/ChannelModel.cpp:459-467`

```cpp
Real TransformAngle0To180(Real x)
{
    x = fmod(x, (Real)360.0);
    if (x > 180.0)
        x = 360.0 - x;
    if (x < 0.0)
        x = -x;
    return x;
}
```

### The Difference

The **v2 version is missing the `else if (x < -180)` case**. This causes **incorrect behavior for angles in the range [-190, -180)**.

### Proof via Test Cases

```
Input: -190°
  Original: 170°  [CORRECT: 360 + (-190) = 170]
  v2:       190°  [INCORRECT: fmod(-190, 360) = -190, then -(-190) = 190]

Input: -185°
  Original: 175°  [CORRECT: 360 + (-185) = 175]
  v2:       185°  [INCORRECT: fmod(-185, 360) = -185, then -(-185) = 185]
```

The v2 function returns values outside the valid [0, 180] range for input angles < -180.

### Impact on ZSD Calculation

When per-ray ZOD angles are computed:
```cpp
Real tx_zenith = ZOD[i] + (3.0 / 8.0) * pow(10.0, mu_ZSD) * offset_angle_rand_coupling[i][j];
ray_AOD[i][j][1] = TransformAngle0To180(tx_zenith);
```

For some rays, `tx_zenith` may fall below -180 (since cluster center ZOD is centered around some baseline, and the offset can be negative and large). When this happens:

- **Original:** Correctly maps to [0, 180] (e.g., -190 → 170)
- **v2:** Returns an invalid value outside [0, 180] (e.g., -190 → 190, which when interpreted as magnitude becomes 190)

This causes v2's ray angles to be pushed toward the extremes (180°) when they should be distributed across the full [0, 180] range.

---

## 4. Impact on Circular Spread Calculation

### Why the Spread Becomes Artificially Small

The circular spread algorithm finds the optimal mean direction that minimizes variance. When ray angles are incorrectly transformed and cluster near boundaries or invalid ranges:

1. **v2 bug creates an artificial concentration** of angles near 180° (instead of having them properly distributed across [0, 180])
2. The circular spread algorithm then finds an optimal mean that minimizes variance around this artificially concentrated distribution
3. **Result:** Much smaller variance, hence much smaller spread (~1.4° vs 14.2° at 95th percentile)

### The Algorithm is Correct

The circular spread algorithm itself works fine. It correctly identifies that the (incorrect) v2 angles have lower variance because they're artificially compressed into a smaller effective range.

---

## 5. ZOD Array Format Verification

Both implementations use the **same ZOD generation process** with the same format:

**File:** `src/channel.cpp:3893-3965` and `v2/src/ChannelSSP.cpp:485-506`

The ZOD[] arrays are **ABSOLUTE ZENITH ANGLES IN GCS (Global Coordinate System)**, NOT offsets. They are centered around values computed from:

```cpp
LOS_ZOD_GCS = 180.0 - LOS_ZOA_GCS;
LOS_ZOA_GCS = 90.0 - atan((bs_height - ms_height) / distance) * (180 / pi);
```

So ZOD[] values typically range from 0° to ~180°, but intermediate calculations can produce values outside this range that need transformation.

---

## 6. Root Cause Summary

| Aspect | Finding |
|--------|---------|
| **Algorithm** | Identical in both versions ✓ |
| **Per-ray angle usage** | Identical in both versions ✓ |
| **Circular spread helper** | Identical in both versions ✓ |
| **Angle transformation** | **BUG in v2: Missing `else if (x < -180)` case** ✗ |
| **ZOD format** | Identical (absolute GCS angles) ✓ |

---

## 7. Recommended Fix

The v2 `TransformAngle0To180()` function in `v2/src/ChannelModel.cpp:459-467` should be corrected to match the original:

```cpp
Real TransformAngle0To180(Real x)
{
    x = fmod(x, (Real)360.0);
    if (x > 180.0)
        x = 360.0 - x;
    else if ((x < 0.0) && (x >= -180.0))
        x = -x;
    else if (x < -180.0)
        x = x + 360.0;
    return x;
}
```

This ensures that all intermediate angle values are correctly mapped to the valid [0, 180] range before being passed to the circular spread calculation.

---

## 8. Verification Steps

To confirm this fix resolves the discrepancy:

1. Apply the fix to `v2/src/ChannelModel.cpp:459-467`
2. Recompile v2
3. Re-run the ZSD circular spread statistics
4. Verify that v2 ZSD values now match the original at the 95th percentile (should be ~14.2° instead of ~1.4°)

---

## Files Involved

| File | Role |
|------|------|
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/Link.cpp` | ORIGINAL: Correct `Transform_angle_0_to_plus_180()` implementation |
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/channel.cpp:4373-4410` | ORIGINAL: Identical circular spread helper |
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/channel.cpp:4412-4454` | ORIGINAL: Per-ray flattening for circular spread |
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelModel.cpp:459-467` | v2: **BUGGY** `TransformAngle0To180()` implementation |
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelModel.cpp:487-519` | v2: Identical circular spread helper |
| `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelSSP.cpp:664-694` | v2: Per-ray flattening for circular spread |

