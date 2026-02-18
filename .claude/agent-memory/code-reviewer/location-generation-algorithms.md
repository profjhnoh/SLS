# Location Generation Algorithms - Original Simulator

## Overview
This document provides the exact algorithms used for BS and MS location generation in the original 5G DL SLS simulator, to be used as reference for v2 implementation.

## BS Location Generation

### InH (Indoor Hotspot) - TYPE 11
**Fixed grid layout**

#### Standard InH (g_mTRP_mode < 2):
- 12 BSs in 2 rows × 6 columns
- **Row 1 (Top)**: y = 35m, x = {-50, -30, -10, 10, 30, 50}
- **Row 2 (Bottom)**: y = 15m, x = {-50, -30, -10, 10, 30, 50}
- Spacing: 20m horizontal, 20m vertical

#### InH with mTRP Mode 1:
- 12 BSs + 5 mTRPs
- mTRPs positioned at: y = 25m, x = {-40, -20, 0, 20, 40}

#### InH with mTRP Mode 2:
- 24 BSs in 3 rows × 8 columns
- **Row 1**: y = 40m, x = {-52.5, -37.5, -22.5, -7.5, 7.5, 22.5, 37.5, 52.5}
- **Row 2**: y = 25m, x = {-52.5, -37.5, -22.5, -7.5, 7.5, 22.5, 37.5, 52.5}
- **Row 3**: y = 10m, x = {-52.5, -37.5, -22.5, -7.5, 7.5, 22.5, 37.5, 52.5}
- Spacing: 15m horizontal, 15m vertical

#### Scaling Factor (InH):
```cpp
bool a = 0;  // false => apply scaling
Real n = 20/5. = 4.0;
if (a == 0) {
    bs[bs_idx].loc.x /= n;
    bs[bs_idx].loc.y /= n;
    ms[ms_idx].loc.x /= n;
    ms[ms_idx].loc.y /= n;
}
```
**Net result**: All InH coordinates scaled down by factor of 4.

### Hexagonal Grid (UMa, UMi, RMa) - TYPE 12, 13

#### 7-Site Configuration (simple_num_BS = 7):
```
cell_radius = ISD * (1 / sqrt(3))

BS locations:
BS[0]: (0, 0)
BS[1]: (0, ISD)
BS[2]: (-1.5*cell_radius, 0.5*ISD)
BS[3]: (-1.5*cell_radius, -0.5*ISD)
BS[4]: (0, -ISD)
BS[5]: (1.5*cell_radius, -0.5*ISD)
BS[6]: (1.5*cell_radius, 0.5*ISD)
```

#### 19-Site Configuration (simple_num_BS = 19):
7-site layout + additional 12 BSs in outer ring:
```
BS[7]:  (0, 2*ISD)
BS[8]:  (-1.5*cell_radius, 1.5*ISD)
BS[9]:  (-3*cell_radius, ISD)
BS[10]: (-3*cell_radius, 0)
BS[11]: (-3*cell_radius, -ISD)
BS[12]: (-1.5*cell_radius, -1.5*ISD)
BS[13]: (0, -2*ISD)
BS[14]: (1.5*cell_radius, -1.5*ISD)
BS[15]: (3*cell_radius, -ISD)
BS[16]: (3*cell_radius, 0)
BS[17]: (3*cell_radius, ISD)
BS[18]: (1.5*cell_radius, 1.5*ISD)
```

#### 1-Site Configuration (simple_num_BS = 1):
- **Not implemented in original code**
- Comment indicates central BS at (0,0) with 3 interfering replicas (actual num_BS = 4)

#### mTRP Mode 1 (Hexagonal):
One mTRP per sector, positioned at ISD/3 from BS center:
```
dist_macro_btw_micro = ISD/3
mTRP_sector0: (dist*cos(30°), dist*sin(30°))
mTRP_sector1: (dist*cos(150°), dist*sin(150°))
mTRP_sector2: (dist*cos(270°), dist*sin(270°))
```

#### mTRP Mode 2 (Hexagonal):
Three mTRPs per sector forming equilateral triangle:
```
dist_macro_btw_sector_center = ISD/3
center_x = dist*cos(30°)
center_y = dist*sin(30°)

dist_micro_wrt_center = (ISD/3)*sqrt(3)/3

Sector 0 mTRPs (relative to sector center):
  mic[0]: center + dist_micro*(cos(60°), sin(60°))
  mic[1]: center + dist_micro*(cos(180°), sin(180°))
  mic[2]: center + dist_micro*(cos(-60°), sin(-60°))

Sectors 1,2: Rotate mic[0-2] by 120° and 240° respectively
```

## MS Location Generation

### InH - Random Uniform in 20m × 25m (or 20m × 15m) Rectangles

#### Standard InH (g_mTRP_mode < 2):
```cpp
// For each BS
int num_bs_with_users = (single_cell_mode == 1) ? 1 : num_BS;
for bs_idx = 0 to num_bs_with_users-1:
    for ms_idx = 0 to num_Indoor_TRxP*num_MS_persector-1:
        do {
            ms.loc.x = (bs.loc.x - 10) + 20 * U(0,1)  // U(0,1) = uniform random
            if (bs_idx < 6):  // Top row
                ms.loc.y = (bs.loc.y - 10) + 25 * U(0,1)
            else:  // Bottom row
                ms.loc.y = (bs.loc.y - 15) + 25 * U(0,1)
        } while (distance(ms, bs) < min_distance)

        ms.Indoor = true
        ms.moving_direction = U(-180, 180)  // degrees
        ms.speed = user_speed * 1000/3600   // m/s
        ms.nearest_bs_idx = bs_idx
```

#### InH with mTRP Mode 2:
```cpp
// Random in entire 120m × 50m area
for ms_idx = 0 to num_BS*num_Indoor_TRxP*num_MS_persector-1:
    ms.loc.x = 120*U(0,1) - 60  // [-60, 60]
    ms.loc.y = 50*U(0,1)         // [0, 50]

    // Find nearest BS (including mTRPs)
    min_dist = infinity
    for bs_idx = 0 to num_BS+num_mTRP-1:
        dist = euclidean_distance(bs[bs_idx].loc, ms.loc)
        if dist < min_dist:
            ms.nearest_bs_idx = bs_idx
            min_dist = dist
```

### Hexagonal Cell - Sector-Based Random Drop

#### Algorithm (MS.cpp: Set_MS_Location(), lines 34-150):
```cpp
cell_radius = ISD * (1 / sqrt(3))

// Step 1: Generate random point in 60° sector wedge
do {
    x_tmp = cell_radius * U(0,1)
    y_tmp = cell_radius * sqrt(3)/2 * U(0,1)

    // Reflect points above diagonal boundary
    if (y_tmp > -sqrt(3)*(x_tmp - cell_radius)) {
        y_tmp2 = cell_radius - x_tmp
        x_tmp2 = y_tmp

        loc.x = x_tmp2*cos(-π/6) - y_tmp2*sin(-π/6)
        loc.y = x_tmp2*sin(-π/6) + y_tmp2*cos(-π/6)
    }
    else {
        loc.x = x_tmp
        loc.y = y_tmp
    }
} while (sqrt(loc.x^2 + loc.y^2) < min_distance)

// Step 2: Rotate to sector orientation
angle = atan(loc.y / loc.x) * 180/π
R = sqrt(loc.x^2 + loc.y^2)

if (sector_idx == 1):
    angle += 120
    loc.x = R*cos(angle*π/180)
    loc.y = R*sin(angle*π/180)
elif (sector_idx == 2):
    angle += 240
    loc.x = R*cos(angle*π/180)
    loc.y = R*sin(angle*π/180)

// Step 3: Translate to BS location (done in Initialdrop)
ms.loc.x += bs[bs_idx].loc.x
ms.loc.y += bs[bs_idx].loc.y
```

**Sector orientation**:
- Sector 0: 0° (pointing right/east)
- Sector 1: 120° (pointing upper-left)
- Sector 2: 240° (pointing lower-left)

### Indoor/Outdoor Decision

#### Dense Urban (TYPE 12):
```cpp
if (carrier_freq > 20 GHz):
    Indoor = false  // All outdoor for 30 GHz+
else:
    if U(0,1) < 0.8:
        Indoor = true   // 80% indoor
    else:
        Indoor = false  // 20% outdoor
```

#### Rural (TYPE 13):
```cpp
if (Configuration_Type == 2):  // Config C
    if U(0,1) < 0.4:
        Indoor = true   // 40% indoor
    else:
        Indoor = false  // 60% outdoor
else:
    if U(0,1) < 0.5:
        Indoor = true   // 50% indoor
    else:
        Indoor = false  // 50% outdoor
```

#### InH (TYPE 11):
```cpp
Indoor = false  // All outdoor (commented code had indoor option)
```

### UE Speed Assignment

#### Dense Urban (TYPE 12):
```cpp
if (Indoor == false):
    speed = 30 km/h  // Outdoor pedestrian
else:
    speed = user_speed  // Indoor (typically 3 km/h)
```

#### Rural (TYPE 13):
```cpp
if (Indoor == false):
    if (Configuration_Type == 2):  // Config C
        if U(0,1) < 0.5:
            speed = user_speed  // 3 km/h outdoor pedestrian (50%)
        else:
            speed = 30 km/h     // Outdoor in-car (50%)
    else:
        speed = 120 km/h  // All outdoor high-speed
else:
    speed = user_speed  // Indoor
```

#### InH (TYPE 11):
```cpp
speed = user_speed * 1000/3600  // m/s
```

## UE Height Assignment

### Dense Urban (TYPE 12)

#### Outdoor UEs:
```cpp
ms.MS_HEIGHT_FINAL = ms_height_out  // 1.5m
ms.floor_idx = 0
ms.ms_d_in = 0
```

#### Indoor UEs (UMi - Config C):
```cpp
// No floor selection - height assigned in Set_LOS_Prob
ms.MS_HEIGHT_FINAL = ms_height_in  // 0 (placeholder, computed later)
ms.floor_idx = 0
```

#### Indoor UEs (UMa - Configs A/B):
```cpp
// Multi-floor building
N = floor(5*U(0,1) + 4)  // Number of floors: 4-8
n = floor(N*U(0,1) + 1)   // Floor index: 1 to N

ms_height_in = 3*(n-1) + 1.5  // Floor height formula

d_in = 25*U(0,1)  // Indoor distance (not MIN for UMa!)
d_out = distance - d_in

ms.MS_HEIGHT_FINAL = ms_height_in
ms.floor_idx = n - 1
ms.ms_d_in = d_in
```

### Rural (TYPE 13)

#### Outdoor UEs:
```cpp
ms.MS_HEIGHT_FINAL = ms_height_out  // 1.5m
ms.floor_idx = 0
```

#### Indoor UEs:
```cpp
d_in = MIN(10*U(0,1), 10*U(0,1))  // Min of two samples

ms.MS_HEIGHT_FINAL = ms_height_in  // 1.5m
ms.floor_idx = 0
```

### InH (TYPE 11)

#### Outdoor UEs (UMi):
```cpp
ms.MS_HEIGHT_FINAL = ms_height_out  // 1.5m
ms.floor_idx = 0
```

#### Indoor UEs (UMi):
```cpp
// Multi-floor building for InH
if (self_bs_idx == 0):  // Only for first BS
    N = floor(5*U(0,1) + 4)  // 4-8 floors
    n = floor(N*U(0,1) + 1)   // Floor 1 to N

    ms_height_in = 3*(n-1) + 1.5

    d_in = MIN(25*U(0,1), 25*U(0,1))  // Min of two samples

    ms.MS_HEIGHT_FINAL = ms_height_in
    ms.floor_idx = 0  // Note: Always 0 in InH despite multi-floor
```

## O2I Loss Calculation

### Material Losses (frequency-dependent):
```cpp
glass_L    = 2 + 0.2*fc_GHz        // Standard multi-pane glass
IRRglass_L = 23 + 0.3*fc_GHz       // Infrared Reflective (IRR) glass
concrete_L = 5 + 4*fc_GHz          // Concrete
wood_L     = 4.85 + 0.12*fc_GHz    // Wood
```

### Dense Urban (TYPE 12)

#### Outdoor:
```cpp
Otoi_loss = 0
```

#### Indoor - Model B (Channel_Model_Type == 1):
```cpp
if U(0,1) > 0.2:  // 80% Low Loss
    Otoi_loss = 5 - 10*log10(0.3*10^(-glass_L/10) + 0.7*10^(-concrete_L/10))
    indoor_L = 0.5 * MIN(25*U(0,1), 25*U(0,1))
    Otoi_loss += indoor_L + 4.4*N(0,1)
    high_loss_flag = 0
else:  // 20% High Loss
    Otoi_loss = 5 - 10*log10(0.7*10^(-IRRglass_L/10) + 0.3*10^(-concrete_L/10))
    indoor_L = 0.5 * MIN(25*U(0,1), 25*U(0,1))
    Otoi_loss += indoor_L + 6.5*N(0,1)
    high_loss_flag = 1
```

#### Indoor - Model A (Channel_Model_Type == 0 or 2):
```cpp
if fc >= 6 GHz:
    if U(0,1) > 0.2:  // 80% Low Loss (ITU-R M.2412), 50% per 3GPP
        Otoi_loss = 5 - 10*log10(0.3*10^(-glass_L/10) + 0.7*10^(-concrete_L/10))
        indoor_L = 0.5 * MIN(25*U(0,1), 25*U(0,1))
        Otoi_loss += indoor_L + 4.4*N(0,1)
        high_loss_flag = 0
    else:  // 20% High Loss
        Otoi_loss = 5 - 10*log10(0.7*10^(-IRRglass_L/10) + 0.3*10^(-concrete_L/10))
        indoor_L = 0.5 * MIN(25*U(0,1), 25*U(0,1))
        Otoi_loss += indoor_L + 6.5*N(0,1)
        high_loss_flag = 1
else:  // fc < 6 GHz
    Otoi_loss = 20
    indoor_L = 0.5 * ms.ms_d_in
    Otoi_loss += indoor_L
```

### Rural (TYPE 13)

#### Outdoor:
```cpp
Otoi_loss = 0
```

#### Indoor - Model B (Channel_Model_Type == 1):
```cpp
// 100% Low Loss
Otoi_loss = 5 - 10*log10(0.3*10^(-glass_L/10) + 0.7*10^(-concrete_L/10))
indoor_L = 0.5 * MIN(10*U(0,1), 10*U(0,1))
Otoi_loss += indoor_L + 4.4*N(0,1)
```

#### Indoor - Model A (Channel_Model_Type == 0):
```cpp
if fc > 6 GHz:
    // 100% Low Loss
    Otoi_loss = 5 - 10*log10(0.3*10^(-glass_L/10) + 0.7*10^(-concrete_L/10))
    indoor_L = 0.5 * MIN(10*U(0,1), 10*U(0,1))
    Otoi_loss += indoor_L + 4.4*N(0,1)
else:  // fc <= 6 GHz
    Otoi_loss = 10
    indoor_L = 0.5 * 10*U(0,1)  // No MIN for <6GHz
    Otoi_loss += indoor_L
```

### InH (TYPE 11):
```cpp
Otoi_loss = 0  // No O2I loss for indoor scenario
```

## In-Car Penetration Loss

### Dense Urban (TYPE 12):
```cpp
incar_loss = 0  // All pedestrians (per 3GPP calibration)
```

### Rural (TYPE 13):
```cpp
if (Indoor == true):
    incar_loss = 0
else:  // Outdoor
    if (Configuration_Type == 2):  // Config C
        if U(0,1) < (0.4/0.6):  // Pedestrian (40% of 60% outdoor)
            incar_loss = 0
        else:  // In-car (20% of 60% outdoor)
            incar_loss = 9 + 5*N(0,1)
    else:  // Config A/B
        incar_loss = 9 + 5*N(0,1)  // 50% all outdoor in-car
```

### InH (TYPE 11):
```cpp
incar_loss = 0
```

## Key Constants by Scenario

### InH_eMBB (TYPE 11):
```
ISD = 20m
min_distance = 0m
ms_height = 1.5m
bs_height = 3m
num_floor = 1
Indoor = false (all outdoor per current code)
```

### Dense_Urban_eMBB (TYPE 12):
```
ISD = 500m (Config A/B), 200m (Config C)
min_distance = 35m (Config A/B), 10m (Config C), 0m (Config D)
ms_height_out = 1.5m
ms_height_in = 0 (computed dynamically)
macro_bs_height = 25m
micro_bs_height = 10m
num_floor = 8
Indoor_probability = 0.8 (below 20 GHz), 0.0 (above 20 GHz)
```

### Rural_eMBB (TYPE 13):
```
ISD = 1732m (Config A/B), 6000m (Config C)
min_distance = 10m
ms_height_out = 1.5m
ms_height_in = 1.5m
macro_bs_height = 35m
micro_bs_height = 10m
num_floor = 1
Indoor_probability = 0.5 (Config A/B), 0.4 (Config C)
```

## Known Bugs in Original Code

### 1. d_2D_in MIN Logic (BUG #5 from MEMORY.md)
```cpp
// Link.cpp:1689, 1711, 1749
indoor_L = 0.5 * MIN(25*U(0,1), 25*U(0,1))
```
**Issue**: Uses minimum of two uniform samples instead of single uniform sample. Per 3GPP TR 38.901 Table 7.4.3-2, d_2D_in should be uniform U(0, d_2D_in_max), not min(U1, U2).

**Correct formula**:
- UMa/UMi: `d_2D_in = 25*U(0,1)` (not MIN)
- RMa: `d_2D_in = 10*U(0,1)` (not MIN)

**Exception**: Multi-floor buildings require special handling - check 3GPP spec for multi-floor d_2D_in distribution.

### 2. UMa Indoor d_in (channel.cpp:756)
```cpp
d_in = 25. * randnum.u();  // Correct (single sample)
```
This is correct for UMa, but contradicts the MIN logic in Link.cpp O2I loss calculation. The O2I loss should use the d_in computed here, not recompute with MIN.

### 3. InH floor_idx Always 0
Despite computing multi-floor heights (n = 1 to N), `floor_idx` is always set to 0 in InH. This may be intentional or a copy-paste error.

## Implementation Notes for v2

1. **Use cell_radius formula**: `cell_radius = ISD / sqrt(3)` everywhere
2. **Sector angles**: 0°, 120°, 240° (not 30°, 150°, 270°)
3. **d_2D_in**: Fix MIN bug - use single uniform sample
4. **O2I loss**: Separate logic for Model A vs Model B, frequency-dependent
5. **Floor height**: `h_UT = 3*(n-1) + 1.5` where n ∈ [1, N], N ∈ [4, 8]
6. **Wraparound**: Original code uses simple_wrap_mat for interference calculation (not documented here)
