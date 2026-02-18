# 3GPP Spec Interpreter Memory

## Pathloss Verification (2026-02-14)
See [pathloss-verification.md](pathloss-verification.md) for full details.

### Critical Bugs Found in v2
1. **InH B NLOS wrong PL1**: `ChannelPathloss.cpp:92-98` uses sub-6 NLOS formula (43.3*log10+11.5) instead of LOS formula (32.4+17.3*log10) for max(PL_LOS, PL'_NLOS)
2. **UMi B NLOS sigma_SF=4.0**: `PathlossParams.cpp:35` should be 7.82 per spec Table 7.4.1-1. Same bug in original code.
3. **RMa LMLC missing -12dB**: `RMa_B_NLOS_LMLC` should be max(PL_LOS, PL_NLOS-12) but -12 offset not implemented

### Calibration-only Limitations
- UMa hE hardcoded to 1.0 (valid for h_UT<=13m only)
- UMa LOS prob omits C(d,h_UT) for h_UT>13m
- UMi/UMa LOS prob doesn't clamp min(18/d,1) (functionally OK)

## Key File Locations (v2)
- Pathloss formulas: `v2/src/ChannelPathloss.cpp`
- Pathloss params/sigma_SF: `v2/tables/PathlossParams.cpp`
- LOS probability: `v2/src/ChannelPathloss.cpp:ComputeLOSProbability()`
- LOS prob configs: `v2/tables/LOSProbParams.cpp`

## Key File Locations (original)
- Pathloss: `src/channel.cpp` Set_PATHLOSS() line 902
- LOS prob: `src/channel.cpp` Set_LOS_Prob() line 557
- Channel params: `src/channel.cpp` Set_Channel_Parameters() line 1558

## ZSA Verification (2026-02-14)
See [zsa-verification.md](zsa-verification.md) for full details.

### Bugs Found
1. **v2 log10(1+fc) vs log10(fc) for UMa NLOS ZSA**: `ChannelSSP.cpp:65` uses `log10(1+fc_GHz)` but UMa NLOS spec uses `log10(fc_GHz)`. ~4% ZSA underestimate at 6GHz.
2. **Questionable c_ZSA=3 for UMa O2I**: `channel.cpp:2610`, `ChannelParams.cpp:410`. Spec may intend c_ZSA=7 (inherited from base table). This narrows ray spread for ~80% of UEs.

### Key Mappings
- ZSA LSP: `generateLSP.cpp:757,766,776` (original), `SimEngine.cpp:670,683,696` (v2)
- ZOA generation: `channel.cpp:3613` Set_ZOAZOD() (original), `ChannelSSP.cpp:456` GenerateZOAZOD() (v2)
- Circular AS: `channel.cpp:4373` compute_circular_angle_spread(), `ChannelSSP.cpp:672` (v2)
- Angle scaling: `channel.cpp:3679-3803` (original), `AngleScalingParams.cpp` (v2)

### UMa ZSA Spec Values (Table 7.5-6)
- LOS: mu=0.96(sub6)/0.95(high), sigma=0.15(sub6)/0.16(high), c_ZSA=7
- NLOS: mu=-0.2856*log10(fc)+1.445(sub6)/-0.3236*log10(fc)+1.512(high), sigma=0.17(sub6)/0.16(high), c_ZSA=7
- O2I: mu=1.01, sigma=0.43, c_ZSA=3(code)/7(spec unclear)

## Spec Reference
- Pathloss: TR 38.901 Table 7.4.1-1 (Section 7.4.1)
- LOS prob: TR 38.901 Table 7.4.2-1 (Section 7.4.2)
- sigma_SF for UMi NLOS (standard): 7.82 dB
- sigma_SF for UMi NLOS (optional): 8.2 dB
- ZSA params: TR 38.901 Table 7.5-6 Part 1 (UMa section)
- ZOA scaling: TR 38.901 Table 7.5-4
- Ray offsets: TR 38.901 Table 7.5-3
