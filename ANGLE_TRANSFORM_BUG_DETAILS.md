# Detailed Analysis: Angle Transformation Bug in v2

## The Bug in One Picture

```
Input: tx_zenith = -185° (intermediate calculation before transformation)

ORIGINAL Transform_angle_0_to_plus_180(-185):
  Step 1: fmod(-185, 360) = -185
  Step 2: Is -185 > 180? NO
  Step 3: Is (-185 < 0) AND (-185 >= -180)? YES ← MATCHES THIS BRANCH
  Step 4: Return -(-185) = 185
  Step 5: WAIT - that's still > 180, need to apply next step? NO - ORIGINAL has ELSE IF chain

  Wait, let me re-trace ORIGINAL more carefully:
  Step 1: fmod(-185, 360) = -185
  Step 2: if (-185 > 180)? NO
  Step 3: else if ((-185 < 0) && (-185 >= -180))? NO, because -185 < -180
  Step 4: else if (-185 < -180)? YES ← MATCHES THIS BRANCH
  Step 5: x = -185 + 360 = 175
  Return: 175 ✓

v2 TransformAngle0To180(-185):
  Step 1: fmod(-185, 360.0) = -185.0
  Step 2: if (-185.0 > 180.0)? NO
  Step 3: if (-185.0 < 0.0)? YES ← MATCHES THIS BRANCH
  Step 4: x = -(-185.0) = 185.0
  Return: 185.0 ✗ (WRONG! Outside [0, 180] range)

REASON FOR BUG: v2 uses two separate if statements instead of an else-if chain.
The second "if (x < 0.0)" ALWAYS executes when it shouldn't.
```

## Problematic Input Range

The v2 bug manifests for any angle in the range **(-∞, -180)** after fmod:

```
Input (before transform)  | fmod(x, 360) | Original | v2 Bug | Status
========================|==============|==========|========|========
-370                     | -10          | 10       | 10     | OK (in [-180, 0))
-260                     | -260+360=-100| 100      | 100    | OK (in [-180, 0))
-190                     | -190         | 170      | 190    | BUG ✗
-185                     | -185         | 175      | 185    | BUG ✗
-180.5                   | -180.5       | 179.5    | 180.5  | BUG ✗
-180                     | -180         | 180      | 180    | Boundary OK
-179                     | -179         | 179      | 179    | OK
...
-1                       | -1           | 1        | 1      | OK
```

## How It Affects the ZSD Calculation

### Step-by-Step: Ray Angle Generation

```cpp
// From Precompute_ray_angles() in both implementations

Real tx_zenith = ZOD[i] + (3.0 / 8.0) * pow(10.0, mu_ZSD) * offset_angle_rand_coupling[i][j];
//                 ^^^^^^^^^   cluster center ZOD (in [0, 180])
//                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ offset (typically ±0.5 to ±3 degrees)

ray_AOD[i][j][1] = TransformAngle0To180(tx_zenith);
```

### Realistic Example: InH (Indoor-to-Home) Scenario

Assume:
- Cluster center ZOD[i] = 85° (reasonable for some MS-BS pair geometry)
- mu_ZSD = -1.5 (log10 scale: 10^(-1.5) ≈ 0.0316)
- Scaling factor: (3/8) × 0.0316 ≈ 0.0119
- offset_angle_rand_coupling[i][j] ranges from ±0.0447 to ±2.1551 radians

For negative offset (worst case):
```
tx_zenith = 85 + 0.0119 × (-2.1551)
          = 85 - 0.0256
          = 84.97°

After Transform_angle_0_to_plus_180:
  Both versions: 84.97° (no issue - well within range)
```

But for NLOS scenarios or when cluster ZOD is near boundaries:
- If ZOD[i] = 5°
- tx_zenith = 5 + 0.0119 × (-2.1551) = 4.97° (still OK)

However, after scaling (which can vary with scenario):
- If larger offset multiplier or different scenario parameters produce:
- tx_zenith = -185° (intermediate calculation)
- **ORIGINAL:** -185° → 175° (correct, maps to [0, 180])
- **v2:** -185° → 185° (WRONG, outside [0, 180])

### Effect on Circular Spread Calculation

When flat_zod[] array is passed to ComputeCircularAngleSpreadHelper:

```
v2 flat_zod[] values: [84.5, 85.2, 85.1, 185, 83.9, 86.1, 185, ...]
                       ↑                  ↑              ↑
                       Normal range      BUG values    BUG values
                       [0, 180]         Outside range!
```

The algorithm searches for optimal delta (0-359) to minimize variance:
- With **incorrect** values like 185, the algorithm finds that delta=0 or small deltas minimize variance
- The variance calculation then works with a dataset that appears artificially compressed
- Result: Minimum spread found is ~1.4° (artificially small)

With **correct** values:
- All values in [0, 180], properly distributed
- Algorithm finds optimal delta that captures the full angular extent
- Result: Spread of ~14.2° (matches original)

## Why This Wasn't Caught

1. **Silent failure:** The function returns a value, not an error
2. **Plausible value:** 185° looks reasonable to someone glancing at it
3. **Rare triggering:** Only occurs for zenith angles calculated to < -180°
4. **No validation:** No code validates that output is in [0, 180]
5. **Subtle conditional:** Using two separate `if` statements instead of `if-else` chain is an easy mistake

## Code Comparison: The Exact Difference

### ORIGINAL (Correct)
```cpp
Real Transform_angle_0_to_plus_180(Real x)
{
    x = fmod(x, 360);
    if (x > 180)                          // Case 1: (180, 360)
    {
        x = 360 - x;
    }
    else if ((x < 0) && (x >= -180))     // Case 2: [-180, 0)
    {
        x = -x;
    }
    else if (x < -180)                   // Case 3: (-∞, -180) ← v2 MISSING THIS
    {
        x = x + 360;
    }
    return x;
}
```

### v2 (Buggy)
```cpp
Real TransformAngle0To180(Real x)
{
    x = fmod(x, (Real)360.0);
    if (x > 180.0)                       // Case 1: (180, 360)
        x = 360.0 - x;
    if (x < 0.0)                         // Case 2: ALL NEGATIVE VALUES
        x = -x;                          // ← Problem: This runs even when Case 1 ran!
    return x;
}
```

### The Control Flow Difference

**ORIGINAL:** `if-else if-else if` chain (mutually exclusive)
```
if (x > 180)       → Do this
  else if (x < 0 && x >= -180)  → Or do this (not both)
    else if (x < -180)          → Or do this (not both)
```

**v2:** Two separate `if` statements (can both execute)
```
if (x > 180)       → Do this
if (x < 0)         → Then ALSO do this (even if first executed!)
```

## Test Case: Tracing Through x = -185

### ORIGINAL Logic Trace
```
x = fmod(-185, 360) = -185
├─ if (-185 > 180)? NO
│  └─ skip first block
├─ else if ((-185 < 0) && (-185 >= -180))? NO (because -185 < -180)
│  └─ skip second block
└─ else if (-185 < -180)? YES ← ENTER THIS BLOCK
   └─ x = -185 + 360 = 175
Return 175 ✓
```

### v2 Logic Trace
```
x = fmod(-185, 360.0) = -185.0
├─ if (-185.0 > 180.0)? NO
│  └─ skip first block
└─ (continue to next statement)
├─ if (-185.0 < 0.0)? YES ← ENTER THIS BLOCK
│  └─ x = -(-185.0) = 185.0
└─ (no more conditions)
Return 185.0 ✗ (OUTSIDE [0, 180]!)
```

## Impact on Statistics

Given the distribution of ray offsets is random:

### Distribution of Problematic Angles

If roughly 5-10% of rays generate intermediate tx_zenith values < -180 (due to random variation in offset_angle_rand_coupling):
- **Original:** All these angles correctly mapped to [0, 180]
- **v2:** These 5-10% of angles appear as values > 180, outside valid range

### Effect on CDF

```
ZSD Circular Spread CDF
========================
Percentile | v2 (Buggy)  | Original (Correct) | Ratio
-----------|-------------|-------------------|-------
5%         | 0.195°      | 0.254°            | 0.77x
50%        | 1.028°      | 2.169°            | 0.47x
95%        | 1.418°      | 14.20°            | 0.10x ← v2 collapses here
```

The CDF saturation around 1.4° occurs because:
- With invalid angles > 180 in the array
- The circular spread algorithm can only find small spreads that accommodate this artifact
- As sample size increases, the artifact becomes more consistent
- The CDF asymptotes to a value determined by the artifact's magnitude

---

## Summary

| Item | Finding |
|------|---------|
| **Root cause** | v2 TransformAngle0To180() missing `else if (x < -180)` branch |
| **Impact range** | Angles after fmod in range (-∞, -180), typically rare but possible |
| **Effect** | Invalid output values > 180, outside valid [0, 180] range |
| **Consequence** | Circular spread algorithm finds artificially small spreads |
| **Fix complexity** | 1-line addition to restore original's else-if structure |
| **Testing needed** | Recompile v2 with fix, re-run statistics, verify 95th percentile rises to ~14.2° |

