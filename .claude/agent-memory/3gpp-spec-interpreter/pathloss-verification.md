# Pathloss Verification Details (2026-02-14)

## Bug #1: InH Model B NLOS wrong max() formula
- File: `v2/src/ChannelPathloss.cpp:92-98`
- Code uses `43.3*log10(d3D)+11.5+20*log10(fc)` (sub-6 NLOS) as PL1
- Should use `32.4+17.3*log10(d3D)+20*log10(fc)` (LOS formula) per spec
- Original code (channel.cpp:960-963) correctly uses LOS formula
- Fix: Separate InH_B_NLOS from InH_NLOS_high case

## Bug #2: UMi Model B NLOS sigma_SF
- File: `v2/tables/PathlossParams.cpp:35`
- Current: `sigma_SF = 4.0` (copied from LOS)
- Spec: `sigma_SF = 7.82` (Table 7.4.1-1 UMi NLOS)
- Original code also has this bug (sigma_SF not set for UMi B NLOS)
- Fix: Change to 7.82

## Bug #3: RMa Config C LMLC
- File: `v2/src/ChannelPathloss.cpp` + `v2/tables/PathlossParams.cpp`
- RMa_B_NLOS_LMLC formula identical to RMa_B_NLOS
- Original: `max(PL_LOS, PL_NLOS - 12)` for Config C
- v2: `max(PL_LOS, PL_NLOS)` without -12
- Fix: Either subtract 12 in formula or add LMLC offset handling in ComputePathloss()

## Verified Correct Items
- All LOS formulas (InH, UMi, UMa, RMa)
- UMa NLOS formula and sigma_SF
- RMa NLOS formula (non-LMLC) and sigma_SF
- All breakpoint distance formulas
- All LOS probability formulas (functionally correct)
- max(LOS,NLOS) logic in ComputePathloss() for NLOS
- O2I sigma_SF values (UMa=7, RMa=8)
