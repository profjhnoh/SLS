# 3GPP 38.901 (TR 38.901 v18.1.0) - Extracted Sections

## Source Document
- File: `/sessions/laughing-pensive-brahmagupta/mnt/Tbps_DLSLS_20251010/doc/38901-j10.docx`
- Sections: 7.3, 7.4, 7.8
- Extraction Date: 2026-02-19

---

# SECTION 7.3: ANTENNA MODELLING

## Paragraphs 532-626

### Key Antenna Concepts

[532] (Heading 2): 7.3 Antenna modelling

[533] (Heading 3): 7.3.0 Antenna array structure

[534] (Normal): This clause captures the antenna array structures considered in this SI for calibration.

[535] (Normal): BS antenna model:

[536] (Normal): The BS antenna is modelled by a uniform rectangular panel array, comprising MgNg panels, as illustrated in Figure 7.3-1 with Mg being the number of panels in a column and Ng being the number of panels in a row. Furthermore the following properties apply:

[537] (B1): - Antenna panels are uniformly spaced in the horizontal direction with a spacing of dg,H and in the vertical direction with a spacing of dg,V.

[538] (B1): - On each antenna panel, antenna elements are placed in the vertical and horizontal direction, where N is the number of columns, M is the number of antenna elements with the same polarization in each column.

[539] (B2): - Antenna numbering on the panel illustrated in Figure 7.3-1 assumes observation of the antenna array from the front (with x-axis pointing towards broad-side and increasing y-coordinate for increasing column number).

[540] (B2): - The antenna elements are uniformly spaced in the horizontal direction with a spacing of dH and in the vertical direction with a spacing of dV.

[541] (B2): - The antenna panel is either single polarized (P =1) or dual polarized (P =2).

[542] (Normal): The rectangular panel array antenna can be described by the following tuple.

[544] (TF): Figure 7.3-1: Cross-polarized panel array antenna model

[545] (Normal): The antenna radiation power pattern of each antenna element is generated according to Table 7.3-1.

[546] (TH): Table 7.3-1: Radiation power pattern of a single antenna element

[548] (Normal): UT antenna model:

[549] (Normal): The UT antenna can be modelled by placing the antenna element/array at the placement candidate locations relative to the device center.

[550] (Normal): For handheld UT devices, the following additional properties apply:

[551] (B1): - Four corners of the rectangle and the center of edges of the rectangle that reflect UT form factor are identified as placement candidate locations.

[552] (B1): - Each antenna is assumed to be oriented along the direction determined by the vector connecting the centre of the rectangle to the candidate location.

[553] (B1): - Reference UT orientation vector of the handheld UT is perpendicular to the plane of the flat UT handheld device, and reflects the way the user typically holds the device.

[554] (B1): - Device dimensions for UT antenna modelling are (15 cm, 7 cm, 0 cm).

[556] (TF): Figure 7.3-2: Handheld UT antenna placement candidate locations relative to centre of device (top down view)

[557] (B1): - Reference handheld UT orientation is illustrated in Figure 7.3-3. The representation the handheld UT orientation will be given as a time-varying Euler angle (ΩUT,α, ΩUT,β, ΩUT,γ) in the UT body-fixed Coordinate System (LCS).

[559] (TF): Figure 7.3-3: Illustration of reference handheld UT orientation in UT LCS

[560] (B2): - For calibration of handheld UT, the following UT rotation based on reference UT orientation is assumed.

[561] (B3): - For calibration with spatial non-stationarity modelling:

[562] (B4): - For one-hand blockage, ΩUT,α = 0 – 360 degrees, ΩUT,β = 45 degrees, ΩUT,γ = 0 degrees,

[563] (B4): - For dual-hand blockage, ΩUT,α = 0 – 360 degrees, ΩUT,β = 0 degrees, ΩUT,γ = 45 degrees,

[564] (B4): - For hand and head blockage, ΩUT,α = 0 – 360 degrees, ΩUT,β = 90 degrees, ΩUT,γ = 0 degrees,

[565] (B3): - For all other cases:

[566] (B4): - ΩUT,α = 0 – 360 degrees, ΩUT,β = 45 degrees, ΩUT,γ = 0 degrees,

[567] (B3): - Example of UT rotations are illustrated in Figure 7.3-4 and 7.3-5.

[569] (TF): Figure 7.3-4: Example of UT rotation of ΩUT,α = 0 - 360 degrees, ΩUT,β = 90 degrees, ΩUT,γ = 0 degrees

[571] (TF): Figure 7.3-5: Example of UT rotation of ΩUT,α = 0 - 360 degrees, ΩUT,β = 0 degrees, ΩUT,γ = 90 degrees

[572] (B1): - The antenna radiation power pattern of each antenna element for handheld UT is generated according to Table 7.3-2, which uses formulas specified in [4].

[573] (TH): Table 7.3-2: Radiation power pattern of a single antenna element for handheld UT

[575] (Normal): For consumer premise equipment (CPE) devices, the following properties apply:

[576] (B1): - Total 9 candidate antenna locations in 4 corners of the vertical square plane, 4 center of the edges of the vertical square plane and the center of the device.

[577] (B1): - Device dimensions for UT antenna modelling are (0 cm, 20 cm, 20 cm)

[579] (TF): Figure 7.3-6: CPE UT antenna placement candidate locations relative to centre of device (side view)

[581] (Heading 3): 7.3.1 Antenna port mapping

[582] (Normal): Legacy BS array antennas, i.e. uniform linear arrays with fix phase shifts between its M elements to obtain a beamtilt in the vertical plane, can be represented by a single antenna port.

[583] (EQ): (7.3-1)

[584] (Normal): where m=1, …, M, θs is the electrical vertical steering angle defined between 0° and 180° (90° represents perpendicular to the antenna array), λ = c/fc is the wavelength, c = 3.0108 m/s is the speed of light.

[585] (Heading 3): 7.3.2 Polarized antenna modelling

[586] (Normal): In general the relationship between radiation field and power pattern is given by:

[587] (EQ): . (7.3-2)

[588] (Normal): The following two models represent two options on how to determine the radiation field patterns based on a defined radiation power pattern.

[589] (Normal): Model-1:

[590] (Normal): In case of polarized antenna elements assume Ψ is the polarization slant angle where 0 degrees corresponds to a purely vertical polarization.

[591] (EQ): , (7.3-3)

[593] (Normal): where

[594] (EQ): , (7.3-3a)

[595] (EQ): . (7.3-3b)

[596] (Normal): Note that the zenith and the azimuth field components Eθ, Eϕ, Hθ and Hϕ are defined in terms of the spherical basis vectors of a spherical coordinate system.

[597] (Normal): Model-2:

[598] (Normal): In case of polarized antennas, the polarization is modelled as angle-independent in both azimuth and elevation, in an LCS-to-GCS transformation.

[599] (EQ): (7.3-4)

[600] (Normal): and

[601] (EQ): , (7.3-5)

[602] (Normal): respectively, where Ψ is the polarization slant angle and F(θ, ϕ) is the 3D antenna element power pattern as a function of azimuth and zenith angles.

[603] (Normal): Handheld UT Model:

[604] (Normal): For cases when a candidate antenna placement location is used for one antenna field pattern:

[605] (B1): - Reference radiation pattern of the UT antenna model is vertically polarized with all the gain in the theta field component.

[606] (Normal): For cases when a candidate antenna placement location is used for two distinct antenna polarization field patterns:

[607] (B1): - Reference radiation pattern of the UT antenna model is,

[608] (B2): - For first antenna field pattern: Ψ = 90° and Etx,θ = 0.

[609] (B2): - For second antenna field pattern: Ψ = 0° and Etx,ϕ = 0.

[610] (Normal): Each polarized field component of the reference radiation pattern Eθ and Eϕ should be rotated according to the orientation and polarization direction of each UT antenna candidate location.

[611] (EQ): , (7.3-6)

[612] (EQ): , (7.3-7)

[613] (EQ): , (7.3-8)

[614] (Normal): where the 3D-rotation angles, Ω'UT,α, Ω'UT,β and Ω'UT,γ, are obtained according to the orientation and polarization direction of each UT antenna candidate location.

[615] (Normal): For cases when a candidate antenna placement location is used for one antenna field pattern (e.g., single polarization):

[616] (B1): - The polarization direction is indicated by the arrow in Figure 7.3-7, which is parallel with the plane of the handheld UT device.

[618] (TF): Figure 7.3-7: Handheld UT antenna polarization directions for one antenna field pattern (top down view)

[619] (Normal): For cases when a candidate antenna placement location is used for two antenna field patterns (e.g., dual polarization) (see note in the figure):

[620] (B1): - For the first antenna field pattern, the polarization direction is indicated by the arrows in Figure 7.3-7 but additionally perpendicular to the plane of the UT device.

[621] (B1): - For the second antenna field pattern, the polarization direction is perpendicular to the polarization direction of the first antenna field pattern.

[622] (B1): - An example for candidate antenna location (6) is given in Figure 7.3-8.

[624] (TF): Figure 7.3-8: Handheld UT polarization direction for two antenna filed pattern (side view)

---

# SECTION 7.4: PATHLOSS, LOS PROBABILITY AND PENETRATION MODELLING

## Paragraphs 626-670

[626] (Heading 2): 7.4 Pathloss, LOS probability and penetration modelling

[627] (Heading 3): 7.4.1 Pathloss

[628] (Normal): The pathloss models are summarized in Table 7.4.1-1 and the distance definitions are indicated in Figure 7.4.1-1 and Figure 7.4.1-2. Note that the distribution of the shadow fading is log-normal, and its standard deviation for each scenario is given in Table 7.4.1-1.

[630] (Normal): Note that

[631] (EQ): (7.4-1)

[633] (TH): Table 7.4.1-1: Pathloss models

[635] (Heading 3): 7.4.2 LOS probability

[636] (Normal): The Line-Of-Sight (LOS) probabilities are given in Table 7.4.2-1.

[637] (TH): Table 7.4.2-1 LOS probability

[639] (Heading 3): 7.4.3 O2I penetration loss

[640] (Heading 4): 7.4.3.1 O2I building penetration loss

[641] (Normal): The pathloss incorporating O2I building penetration loss is modelled as in the following:

[642] (EQ): (7.4-2)

[643] (Normal): where PLb is the basic outdoor path loss given in Clause 7.4.1, where d is replaced by d'O2I. L_w is the building penetration loss through the external wall, L_in is the inside loss dependent on the depth into the building, and σP is the standard deviation for the penetration loss.

[644] (Normal): L_w is characterized as:

[645] (EQ): (7.4-3)

[646] (Normal): L_mat is an additional loss is added to the external wall loss to account for non-perpendicular incidence; L_mat,i is the penetration loss of material i, example values of which can be found in Table 7.4.3-1; p_i is proportion of i-th materials, where sum(p_i)=1; and N is the number of materials.

[647] (TH): Table 7.4.3-1: Material penetration losses

[649] (Normal): Table 7.4.3-2 gives μ_w, σ_in and σP for two O2I penetration loss models. The O2I penetration is UT-specifically generated, and is added to the SF realization in the log domain.

[650] (TH): Table 7.4.3-2: O2I building penetration loss model

[652] (Normal): d'O2I is minimum of two independently generated uniformly distributed variables between 0 and 25 m for UMa and UMi-Street Canyon, between 0 and 25 m for SMa UTs in commercials buildings, between 0 and 10 m for SMa UTs in residential buildings, and between 0 and 10 m for RMa. d'O2I shall be UT-specifically generated.

[653] (Normal): Both low-loss and high-loss models are applicable to UMa and UMi-Street Canyon.

[654] (Normal): Only the low-loss model is applicable to RMa.

[655] (Normal): Only the high-loss model is applicable to InF.

[656] (Normal): Low-loss, high-loss, and low-loss A models are applicable to SMa.

[657] (Normal): The composition of low and high loss is a simulation parameter that should be determined by the user of the channel models, and is dependent on the use of metal-coated glass in buildings and the deployment scenarios. Such use is expected to differ in different markets and regions of the world and also may increase over years to new regulations and energy saving initiatives. Furthermore, the use of such high-loss glass currently appears to be more predominant in commercial buildings than in residential buildings in some regions of the world (see note).

[658] (NO): Note: One example survey for the US market can be found in [5]. The survey does not necessarily be representative for all the scenarios. Other ratios outside of the survey should not be precluded.

[659] (Normal): For backwards compatibility with TR 36.873 [3], the following building penetration model should be used for UMa and UMi single-frequency simulations at frequencies below 6 GHz.

[660] (TH): Table 7.4.3-3. O2I building penetration loss model for single-frequency simulations <6 GHz

[662] (Heading 4): 7.4.3.2 O2I car penetration loss

[663] (Normal): The pathloss incorporating O2I car penetration loss is modelled as in the following:

[664] (EQ): (7.4-4)

[665] (Normal): where PLb is the basic outdoor path loss given in Clause 7.4.1. μ = 9, and σP = 5. The car penetration loss shall be UT-specifically generated. Optionally, for metallized car windows, μ = 20 can be used. The O2I car penetration loss models are applicable for at least 0.6-60 GHz.

[666] (Heading 3): 7.4.4 Autocorrelation of shadow fading

[667] (Normal): The long-term (log-normal) fading in the logarithmic scale around the mean path loss PL (dB) is characterized by a Gaussian distribution with zero mean and standard deviation. Due to the slow fading process versus distance x (x is in the horizontal plane), adjacent fading values are correlated. Its normalized autocorrelation function R(x) can be described with sufficient accuracy by the exponential function ITU-R Rec. P.1816 [18]

[668] (EQ): (7.4-5)

[669] (Normal): with the correlation length dcor being dependent on the environment, see the correlation parameters for shadowing and other large scale parameters in Table 7.5-6 (Channel model parameters). In a spatial consistency procedure in Clause 7.6.3, the cluster specific random variables are also correlated following the exponential function with respect to correlation distances in the two dimensional horizontal plane.

---

## TABLE 7.4.1-1: PATHLOSS MODELS

| Scenario | LOS/NLOS | Pathloss Formula | Shadow Fading Std [dB] | Applicability Range, Antenna Height Default Values |
|----------|----------|-----------------|------------------------|-----------------------------------------------------|
| RMa | LOS | PL = 20 log10(fc) + 32.45 + 20 log10(d), see note 5 | For 10 m ≤ d2D < dBP, For dBP ≤ d2D < 10 km, | h = avg. building height W = avg. street width |
| RMa | NLOS | PL = 20 log10(fc) + 32.45 + 20 log10(d), for d > 10 m | | h = avg. building height W = avg. street width |
| UMa | LOS | PL = 28.0 + 22 log10(d) + 20 log10(fc), see note 1 | | See note 7 |
| UMa | NLOS | PL = 13.54 + 39.08 log10(d) + 20 log10(fc), for d > 10 m | | Explanations: see note 3 |
| UMa | NLOS | Optional | | |
| UMi - Street Canyon | LOS | PL = 32.4 + 21 log10(d) + 20 log10(fc), see note 1 | | See note 7 |
| UMi - Street Canyon | NLOS | PL = 35.46 + 38.08 log10(d) + 20 log10(fc), for d > 10 m | | Explanations: see note 4 |
| UMi - Street Canyon | NLOS | Optional | | |
| InH - Office | LOS | PL = 17.3 + 38 log10(d) + 20 log10(fc) | | |
| InH - Office | NLOS | PL = 20.05 + 46.3 log10(d) + 20 log10(fc) | | |
| InH - Office | NLOS | Optional | | |
| InF | LOS | PL = 20 log10(d) + 32.45 + 20 log10(fc) | | |
| InF | NLOS | InF-SL: PL = 38.5 + 20 log10(d) + 20 log10(fc) | | |
| InF | NLOS | InF-DL: PL = 38.3 + 20 log10(d) + 20 log10(fc) | | |
| InF | NLOS | InF-SH: PL = 35.7 + 20 log10(d) + 20 log10(fc) | | |
| InF | NLOS | InF-DH: PL = 35.5 + 20 log10(d) + 20 log10(fc) | | |
| SMa | LOS | PL = 20 log10(fc) + 32.45 + 20 log10(d), PL = 20 log10(fc) + 32.45 + 40 log10(d) | For 10 m ≤ d2D < dBP, For dBP ≤ d2D < 5000 m, | 10 m ≤ d2D < dBP dBP ≤ d2D < 5000 m 25 m < hBS < 35 m h = 10 m |
| SMa | NLOS | PL = 20 log10(fc) + 32.45 + 20 log10(d) + Lclutter | | 10 m ≤ d2D < 5000 m 1 m < hUT < 14 m h = 10 m, W = 10 m |

**NOTE 1:** Effective breakpoint distance d'BP = 4 h'BS h'UT fc/c, where fc is the centre frequency in Hz, c = 3.0108 m/s is the propagation velocity in f...

---

## TABLE 7.4.2-1: LOS PROBABILITY

| Scenario | LOS Probability Formula |
|----------|------------------------|
| RMa | PrLOS = min(18/d, 1) * (1 - exp(-d/36)) + exp(-d/36) |
| UMi - Street canyon | PrLOS = (18/d) * (1 - exp(-d/36)) + exp(-d/36), where d is 2D distance in meters |
| UMa | PrLOS = (18/d) * (1 - exp(-d/63)) + exp(-d/63), where d is 2D distance in meters |
| Indoor - Mixed office | PrLOS = 1 |
| Indoor - Open office | PrLOS = 1 |
| InF-SL, InF-SH, InF-DL, InF-DH | Where parameters, and are defined in Table 7.2-4 |
| InF-HH | PrLOS = 1 |
| SMa | Where parameters = 20 m, 8 m, = 15 m with vegetation options (0%, 10%, 20%) |

---

# SECTION 7.8: CHANNEL MODEL CALIBRATION

## Paragraphs 1481-1530

[1481] (Heading 2): 7.8 Channel model calibration

[1482] (Heading 3): 7.8.1 Large scale calibration

[1483] (Normal): For large scale calibration, fast fading is not modeled. The calibration parameters can be found in Table 7.8-1. The calibration results can be found in R1-2306406.

[1484] (TH): Table 7.8-1: Simulation assumptions for large scale calibration

[1486] (Normal): Additional calibration parameters can be found in Table 7.8-1A. It is assumed that parameters from Table 7.8-1 is used in the simulations unless specified elsewhere.

[1487] (TH): Table 7.8-1A: Simulation assumptions for large scale calibration

[1489] (Heading 3): 7.8.2 Full calibration

[1490] (Normal): The calibration parameters for full calibration including the fast fading modelling can be found in Table 7.8-2. Unspecified parameters are assumed to be similar to the large scale calibration parameters from Table 7.8-1 unless otherwise specified.

[1491] (TH): Table 7.8-2: Simulation assumptions for full calibration

[1493] (Normal): Additional full calibration parameters can be found in Table 7.8-2A. It is assumed that parameters from Table 7.8-2 is used unless:

[1494] (B1): - SCS of 15 kHz

[1495] (B1): - UT attachment is based on RSRP (formula) from BS port 0

[1496] (B1): - BS antenna configuration 1 and 2 both apply

[1497] (B1): - UT antenna configuration, pattern, and polarization modeling are labeled as UT antenna config A

[1498] (Normal): The calibration results based on additional calibration parameters can be found in R1-2506406.

[1499] (TH): Table 7.8-2A: Simulation assumptions for full calibration

[1501] (Heading 3): 7.8.3 Calibration of additional features

[1502] (Normal): The calibration parameters for the calibration of oxygen absorption, large bandwidth and large antenna array, spatial consistency, blockage, and indoor factory scenario are defined below.

[1503] (Normal): The additional calibration parameters for the calibration of near field channel modeling, BS side spatial non-stationarity, and UT side spatial non-stationarity, are provided in Tables 7.8-6A, 7.8-6B, and 7.8-6C, respectively, where near field distance d_nf is:

[1504] (B1): - Up to 1.5 m for UMa with maximum antenna elements in the array is 5k for single Polarization.

[1505] (B1): - Up to 1 m for UMi with maximum antenna elements in the array is 2.22k for single Polarization.

[1506] (B1): - Up to 0.71 m for Indoor factory with maximum antenna elements in the array is 1.12k for single Polarization.

[1507] (B1): - Up to 0.25 (for rectangular antenna array), 0.5 (for linear antenna array) m for Indoor office with maximum antenna elements in the array is 500 for single Polarization.

[1508] (Normal): The additional calibration results can be found in R1-2506406.

[1509] (TH): Table 7.8-3: Simulation assumptions for calibration for oxygen absorption

[1511] (TH): Table 7.8-4: Simulation assumptions for calibration for large bandwidth and large antenna array

[1513] (TH): Table 7.8-5: Simulation assumptions for calibration for spatial consistency

[1515] (TH): Table 7.8-6: Simulation assumptions for calibration for blockage

[1517] (TH): Table 7.8-6A: Simulation assumptions for calibration for near field channel modeling

[1519] (TH): Table 7.8-6B: Simulation assumptions for calibration for BS side spatial non-stationarity

[1521] (TH): Table 7.8-6C: Simulation assumptions for calibration for UT side spatial non-stationarity

[1523] (Heading 3): 7.8.4 Calibration of the indoor factory scenario

[1524] (Normal): For the InF, the calibration parameters can be found in Table 7.8-7. The calibration results can be found in R1-1909704.

[1525] (Normal): It should be noted absolute delay model had not been agreed by the deadline, so companies were not able to submit CDF of absolute delay for InF scenario.

[1526] (TH): Table 7.8-7: Simulation assumptions for large scale calibration for the indoor factory scenario

[1528] (Heading 2): 7.9 Channel model(s) for ISAC

[1529] (Heading 3): 7.9.0 Introduction

---

## KEY CALIBRATION PARAMETERS SUMMARY

### Near Field Distances (Section 7.8.3)
- **UMa**: Up to 1.5 m (max 5k antenna elements, single polarization)
- **UMi**: Up to 1 m (max 2.22k antenna elements, single polarization)
- **Indoor Factory**: Up to 0.71 m (max 1.12k antenna elements, single polarization)
- **Indoor Office**: Up to 0.25 m (rectangular) or 0.5 m (linear) arrays (max 500 antenna elements, single polarization)

### Calibration Scope
The section covers calibration parameters for:
1. Large scale calibration (Table 7.8-1)
2. Full calibration with fast fading (Table 7.8-2)
3. Oxygen absorption calibration (Table 7.8-3)
4. Large bandwidth and antenna array calibration (Table 7.8-4)
5. Spatial consistency calibration (Table 7.8-5)
6. Blockage calibration (Table 7.8-6)
7. Near field channel modeling (Table 7.8-6A)
8. BS side spatial non-stationarity (Table 7.8-6B)
9. UT side spatial non-stationarity (Table 7.8-6C)
10. Indoor factory scenario calibration (Table 7.8-7)

### Default Calibration Assumptions
- SCS (Subcarrier Spacing): 15 kHz
- UT attachment: Based on RSRP from BS port 0
- BS antenna configurations: 1 and 2
- UT antenna configuration: Labeled as "UT antenna config A"

---

# KEY FORMULAS REFERENCED

## Pathloss Calculation
PL(d) = L₀ + 10·n·log₁₀(d) + 20·log₁₀(fc) + Xσ

Where:
- L₀: Reference pathloss at 1 meter
- n: Pathloss exponent
- d: Distance in meters
- fc: Center frequency in GHz
- Xσ: Shadow fading (log-normal distribution)

## Breakpoint Distance (RMa, SMa scenarios)
d'BP = 4·h'BS·h'UT·fc/c

Where:
- h'BS: BS antenna height (m)
- h'UT: UT antenna height (m)
- fc: Center frequency (Hz)
- c: Speed of light = 3.0108 m/s

## Autocorrelation of Shadow Fading
R(x) = exp(-x/dcor)

Where:
- x: Distance in horizontal plane
- dcor: Correlation distance (scenario-dependent)

## O2I Building Penetration
PL_O2I = PLb + L_w + L_in + Xσ

Where:
- PLb: Basic outdoor pathloss
- L_w: External wall penetration loss
- L_in: Inside loss (depth-dependent)
- Xσ: Penetration loss standard deviation

---

# DOCUMENT REFERENCES

**3GPP Technical Report:** 3GPP TR 38.901 v18.1.0
**Title:** "Study on channel model for frequency spectrum above 6 GHz"
**Related Standards:**
- 3GPP TR 36.873: "Study on 3D channel model for LTE"
- 3GPP TS 38.104: "Base Station (BS) radio transmission and reception"
- ITU-R Rec. P.1816: "Spatial profile and time variation of broadband land mobile services"
- ITU-R Rec. P.2040-1: "Effects of building materials and structures on radiowave propagation"

