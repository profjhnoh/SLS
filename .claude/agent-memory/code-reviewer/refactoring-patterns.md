# Refactoring Patterns - Get_signal_interference v2

## Overview
Reviewed v2 refactoring (2026-02-12): Get_signal_interference() split into 11 functions.

## Critical Findings

### 1. Missing link_antgain Assignment (CRITICAL BUG)
**Location**: select_serving_cell_v2 (Line.cpp:3892-3940)
**Issue**: Original code stores `link_antgain = RSRP_antgain` at lines 347, 384, 547, 585, but v2 does NOT set this field.
**Impact**: link_antgain left at 0 (default value), may break downstream code that depends on it.
**Fix Required**: Add `link_antgain = srv.max_rsrp_dB;` in select_serving_cell_v2() around line 3917.

### 2. Typo Bug FIXED in v2
**Location**: Original line 333
**Original Bug**: `analog_beam_selection[bs_idx].sector_a = sector_z;` (assigns z to a field)
**v2 Fix**: Lines 3880-3881 correctly assign both fields separately.
**Verdict**: v2 FIXES this bug!

### 3. InH 1TRxP Config A Interference Logic
**Original**: Lines 368-410 accumulate interference inline during serving cell selection
**v2**: Lines 3951-3963 compute post-hoc from candidates vector
**Analysis**: EQUIVALENT - both sum all non-serving cells with max beam. Order differs but result identical.

### 4. InH 3TRxP g_comp_mode Code
**Original**: Lines 741-759 handle g_comp_mode if enabled
**v2**: Lines 4010-4030 handle same logic
**Issue**: Original line 743 uses `signal` variable (last loop value from line 792), v2 line 4012 also uses `signal`.
**Concern**: In v2, `signal` is not in scope at line 4012 - this will NOT compile!
**Status**: CRITICAL - v2 code references undefined variable `signal` in comp_mode block.

### 5. InH 3TRxP comp_interf_strength[0]
**Original**: Lines 754-759 unconditionally compute comp_interf_strength[0]
**v2**: Lines 4024-4029 compute same, but only if `cfg.sectors_per_bs == 3`
**Analysis**: CORRECT - comp variables only relevant for 3TRxP scenario.

## Design Quality

### Strengths
- Clean separation of concerns (config, pathloss, beam search, selection, interference)
- Reusable CandidateCell struct avoids repeated channel lookups
- ScenarioConfig encapsulates branching logic
- Better testability (each function can be unit tested)
- Fixes original sector_a/sector_z typo bug

### Concerns
- **Backwards compatibility**: link_antgain field not set
- **Compilation error**: `signal` undefined in g_comp_mode block (line 4012)
- **RNG order changed**: randnum.u() call order differs → different per-UE results (statistically OK, but not bit-exact)
- **Memory**: std::vector allocation overhead (minor, ~10KB for 300 candidates)

## Behavioral Equivalence

| Scenario | Status | Notes |
|----------|--------|-------|
| InH 1TRxP Config A | ✓ PASS | Post-hoc interference sum equivalent |
| InH 1TRxP Config B | ✓ PASS | Random beam logic preserved |
| InH 3TRxP Config A | ✓ PASS | Matches original |
| InH 3TRxP Config B | ✓ PASS | Matches original |
| InH 3TRxP comp_mode | ✗ FAIL | Undefined `signal` variable |
| Dense Urban/Rural | ✓ PASS | RSRP table lookup preserved |
| mTRP interference | ✓ PASS | Correct logic |
| single_cell_mode | ✓ PASS | Correct filtering |

## Recommendations

### MUST FIX (Blocking Issues)
1. **Line 4012**: Replace `signal` with computed value from comp_sector
2. **Line 3917**: Add `link_antgain = srv.max_rsrp_dB;`

### SHOULD FIX (Compatibility)
3. Add regression test comparing original vs v2 outputs for all 3 scenarios
4. Document that RNG order differs (not bit-exact but statistically equivalent)

### NICE TO HAVE
5. Add const correctness to helper functions (build_scenario_config_v2, etc.)
6. Consider using std::array instead of std::vector (fixed size known at compile time)
7. Add comments explaining Config A vs Config B interference methods

## Testing Strategy
- Run both versions side-by-side with fixed RNG seed
- Compare: self_bs_idx, self_sector_idx, link_pathloss, signal, interference, geometry
- Test all 3 scenarios: InH 1TRxP, InH 3TRxP, Dense Urban
- Test edge cases: empty candidate list (shouldn't happen), single_cell_mode enabled
