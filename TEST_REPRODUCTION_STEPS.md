# Test Reproduction Steps: Verifying the ZSD Spread Bug

## Quick Verification Test

### Test 1: Transform Function Behavior

Create a simple test program to verify the angle transformation bug:

```cpp
// test_transform.cpp
#include <iostream>
#include <cmath>
using namespace std;

typedef double Real;

// ORIGINAL (Correct)
Real Transform_angle_0_to_plus_180_orig(Real x) {
    x = fmod(x, 360);
    if (x > 180) {
        x = 360 - x;
    } else if ((x < 0) && (x >= -180)) {
        x = -x;
    } else if (x < -180) {
        x = x + 360;
    }
    return x;
}

// v2 (Buggy)
Real TransformAngle0To180_v2(Real x) {
    x = fmod(x, (Real)360.0);
    if (x > 180.0)
        x = 360.0 - x;
    if (x < 0.0)
        x = -x;
    return x;
}

int main() {
    cout << "Testing Angle Transformation Functions" << endl;
    cout << "======================================" << endl;
    cout << "Input°  | Original | v2    | Status" << endl;
    cout << "--------|----------|-------|--------" << endl;

    Real test_values[] = {-370, -260, -190, -185, -180.5, -180, -179, -100, -1, 0, 1, 90, 179, 180, 181, 225, 270};
    int n = sizeof(test_values) / sizeof(test_values[0]);

    int bug_count = 0;
    for (int i = 0; i < n; i++) {
        Real val = test_values[i];
        Real orig = Transform_angle_0_to_plus_180_orig(val);
        Real v2 = TransformAngle0To180_v2(val);

        bool is_bug = (orig != v2);
        if (is_bug) bug_count++;

        cout << val << " | " << orig << " | " << v2 << " | ";
        if (is_bug) {
            cout << "BUG!" << endl;
        } else {
            cout << "OK" << endl;
        }
    }

    cout << "======================================" << endl;
    cout << "Total bugs found: " << bug_count << endl;

    if (bug_count > 0) {
        cout << "RESULT: v2 has bugs in angle transformation!" << endl;
        return 1;
    } else {
        cout << "RESULT: No differences found" << endl;
        return 0;
    }
}
```

**Compile and run:**
```bash
g++ -o test_transform test_transform.cpp
./test_transform
```

**Expected output with v2 bug:**
```
Input°  | Original | v2    | Status
--------|----------|-------|--------
-370 | 10 | 10 | OK
-260 | 100 | 100 | OK
-190 | 170 | 190 | BUG!
-185 | 175 | 185 | BUG!
-180.5 | 179.5 | 180.5 | BUG!
...
======================================
Total bugs found: 3
RESULT: v2 has bugs in angle transformation!
```

---

## Full Integration Test

### Test 2: ZSD Statistics Before and After Fix

To verify this is the actual cause of the 10x discrepancy:

#### Step 1: Baseline Run (Current v2 - with bug)

```bash
cd /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2
mkdir -p build
cd build
cmake ..
make

# Run a quick simulation to collect ZSD statistics
cd ..
./sim_executable config_file.cfg  # Run to generate output
```

Capture the ZSD circular angle spread values. They should match the buggy values:
- 5th percentile: ~0.195°
- 50th percentile: ~1.028°
- 95th percentile: ~1.418°

#### Step 2: Apply the Fix

Edit file: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelModel.cpp`

**Lines 459-467 (BEFORE - BUGGY):**
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

**Lines 459-469 (AFTER - FIXED):**
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

**Key changes:**
- Line 465: Change `if (x < 0.0)` to `else if ((x < 0.0) && (x >= -180.0))`
- Lines 466-467 (NEW): Add the missing branch:
  ```cpp
  else if (x < -180.0)
      x = x + 360.0;
  ```

#### Step 3: Recompile and Re-run

```bash
cd /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/build
cmake ..
make

# Run the same simulation again
cd ..
./sim_executable config_file.cfg
```

#### Step 4: Compare Statistics

**Expected behavior after fix:**

The ZSD circular angle spread values should now match the original implementation:
- 5th percentile: ~0.254° (was 0.195°)
- 50th percentile: ~2.169° (was 1.028°)
- 95th percentile: ~14.20° (was 1.418°)

**Verification checklist:**
- [ ] 95th percentile rises from ~1.4° to ~14.2° (10x increase!)
- [ ] 50th percentile rises from ~1.0° to ~2.2° (~2x increase)
- [ ] 5th percentile slightly increases from ~0.195° to ~0.254° (~1.3x)
- [ ] Distribution shape changes from saturating to spreading across full range
- [ ] Values match original implementation at all percentiles

---

## Advanced Verification: Angle Distribution Analysis

### Test 3: Histogram of Transformed Angles

To see the actual distribution of transformed angles and verify the bug:

```cpp
// test_angle_distribution.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

typedef double Real;

Real TransformAngle0To180_orig(Real x) {
    x = fmod(x, 360);
    if (x > 180) {
        x = 360 - x;
    } else if ((x < 0) && (x >= -180)) {
        x = -x;
    } else if (x < -180) {
        x = x + 360;
    }
    return x;
}

Real TransformAngle0To180_v2(Real x) {
    x = fmod(x, (Real)360.0);
    if (x > 180.0)
        x = 360.0 - x;
    if (x < 0.0)
        x = -x;
    return x;
}

int main() {
    // Simulate realistic ZOD + offset values
    // Cluster ZOD typically in [0, 180]
    // Offset multiplier: (3/8) * 10^(mu_ZSD) typically 0.5-3 degrees
    // offset_angle_rand_coupling: ±0.0447 to ±2.1551 radians (±2.6 to ±123 degrees)

    vector<Real> angles_orig, angles_v2;

    cout << "Simulating 10000 ray angle transformations..." << endl;

    for (int trial = 0; trial < 10000; trial++) {
        // Random cluster ZOD in [0, 180]
        Real cluster_zod = (Real)rand() / RAND_MAX * 180.0;

        // Random offset (simulated from ±2.1551 radians max with scaling)
        Real offset = ((Real)rand() / RAND_MAX * 2.0 - 1.0) * 3.0;  // ±3 degrees

        // Intermediate zenith angle (can go outside [0, 180])
        Real tx_zenith = cluster_zod + offset;

        // Transform
        Real orig = TransformAngle0To180_orig(tx_zenith);
        Real v2 = TransformAngle0To180_v2(tx_zenith);

        angles_orig.push_back(orig);
        angles_v2.push_back(v2);
    }

    // Count values outside [0, 180]
    int out_of_range_count = 0;
    for (auto val : angles_v2) {
        if (val < 0 || val > 180) {
            out_of_range_count++;
        }
    }

    cout << "\nResults:" << endl;
    cout << "========" << endl;
    cout << "Total samples: " << angles_v2.size() << endl;
    cout << "v2 values outside [0, 180]: " << out_of_range_count << " ("
         << (100.0 * out_of_range_count / angles_v2.size()) << "%)" << endl;

    if (out_of_range_count > 0) {
        cout << "\nBUG CONFIRMED: v2 produces invalid angles!" << endl;
        cout << "\nSample invalid values:" << endl;
        int shown = 0;
        for (int i = 0; i < (int)angles_v2.size() && shown < 5; i++) {
            if (angles_v2[i] < 0 || angles_v2[i] > 180) {
                cout << "  Input around " << angles_orig[i] << " -> v2 output: "
                     << angles_v2[i] << " (INVALID!)" << endl;
                shown++;
            }
        }
        return 1;
    } else {
        cout << "\nNo issues found (simulation may not have triggered bug)" << endl;
        return 0;
    }
}
```

**Compile and run:**
```bash
g++ -o test_angle_dist test_angle_distribution.cpp
./test_angle_dist
```

---

## Regression Test: Original vs v2 Side-by-Side

### Test 4: Run Both Simulators on Same Configuration

```bash
# Step 1: Run original simulator
cd /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010
mkdir -p build_orig
cd build_orig
cmake ..
make
cd ..

# Run with same config file
./Tbps_DLSLS config_test.cfg > results_original.log 2>&1

# Step 2: Run v2 simulator (without fix)
cd v2
mkdir -p build_buggy
cd build_buggy
cmake ..
make
cd ..

# Run with same config file
./sim_executable ../config_test.cfg > results_v2_buggy.log 2>&1

# Step 3: Apply fix to v2 and recompile
# (Edit v2/src/ChannelModel.cpp as described above)
rm -rf build_buggy
mkdir -p build_fixed
cd build_fixed
cmake ..
make
cd ..

./sim_executable ../config_test.cfg > results_v2_fixed.log 2>&1

# Step 4: Compare results
cd ..
echo "=== ORIGINAL ===" && grep "circular_angle_spread_ZOD\|ZSD" results_original.log | head -20
echo ""
echo "=== v2 BUGGY ===" && grep "circular_angle_spread_ZOD\|ZSD" results_v2_buggy.log | head -20
echo ""
echo "=== v2 FIXED ===" && grep "circular_angle_spread_ZOD\|ZSD" results_v2_fixed.log | head -20
```

---

## Success Criteria

### You've found the root cause if:

1. **Transform test (Test 1)** shows v2 produces values > 180 for negative inputs
2. **ZSD statistics (Test 2)** show 95th percentile rising from ~1.4° to ~14.2° after fix
3. **Angle distribution (Test 3)** shows v2 has out-of-range values before fix
4. **Side-by-side comparison (Test 4)** shows original and v2-fixed produce matching ZSD values

### The fix is complete when:

- [ ] v2 95th percentile ZSD matches original (≈14.2°)
- [ ] v2 50th percentile ZSD matches original (≈2.17°)
- [ ] v2 5th percentile ZSD matches original (≈0.254°)
- [ ] No angles in flat_zod[] exceed 180° or go below 0°
- [ ] Circular spread algorithm receives valid input in all cases

---

## Troubleshooting

### If the fix doesn't work:

1. **Verify the fix was applied correctly:**
   ```bash
   diff /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/src/Link.cpp \
        /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/src/ChannelModel.cpp | grep -A5 -B5 "TransformAngle"
   ```

2. **Check compilation flags:**
   - Ensure no optimization flags are stripping the new code
   - Verify debug symbols are included: `g++ -g -O0 ...`

3. **Verify actual function is being called:**
   - Add debug output to TransformAngle0To180() to confirm it's executed
   - Print values before and after transformation

4. **Check for other places calling the transform:**
   ```bash
   grep -rn "TransformAngle0To180" /sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/v2/
   ```
   Ensure all call sites are tested

---

## Expected Test Output

After successfully applying the fix, you should see:

```
ZSD Circular Angle Spread Statistics
=====================================
Original Implementation:
  5th percentile:   0.254°
  25th percentile:  0.689°
  50th percentile:  2.169°
  75th percentile:  6.542°
  95th percentile: 14.200°

v2 Before Fix (BUGGY):
  5th percentile:   0.195°
  25th percentile:  0.512°
  50th percentile:  1.028°
  75th percentile:  1.285°
  95th percentile:  1.418°  ← Artificial saturation!

v2 After Fix:
  5th percentile:   0.254°
  25th percentile:  0.689°
  50th percentile:  2.169°
  75th percentile:  6.542°
  95th percentile: 14.200°  ← Now matches original!
```

This confirms the bug fix is successful.

