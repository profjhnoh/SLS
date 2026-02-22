# Code Refactorer Memory - Tbps_DLSLS

## Project Build
- MinGW GCC 15.1.0, build dir: `build/`, command: `cmake --build . -j8`
- Executable: `Tbps_DLSLS.exe` in project root

## Key File Locations
- Pathloss: `channel.cpp` `Set_PATHLOSS()` ~line 889
- LOS probability: `channel.cpp` `Set_LOS_Prob()` ~line 544
- Channel parameters: `channel.cpp` `Set_Channel_Parameters()` ~line 1545
- Signal/RSRP: `Link.cpp` lines 833-950 (Dense Urban serving cell)
- Antenna pattern: `Link.cpp` `Get_BS_antenna_pattern()` ~line 2571
- Field pattern: `Link.cpp` `Get_BS_antenna_field_pattern()` ~line 2751
- O2I loss: `Link.cpp` ~lines 1692-1716

## Dead Code Patterns
- `g_comp_mode` is always 0 — all `if (g_comp_mode)` blocks are dead code
- Fully removed from (2026-02-20):
  - `generateSSP.cpp:75` — if/else; kept else branch, dropped comp_sector_idx branch
  - `channel_update.cpp:594` — entire `if (comp_mode && g_comp_mode)` block
  - `receive_downlink.cpp:1168` — if/else; kept else body, removed inner duplicate decls
  - `scheduling.cpp:169,179` — whole outer `if (comp_mode==1)` block (both inner ifs dead)
  - `measure.cpp:57` — `comp_ue_pct` replaced with hardcoded `5`; trailing comments cleaned
  - `PerdropStatistics.cpp:206-208` — 3 output lines for dead globals removed
- When removing if/else where if-branch is dead: promote else-body, remove duplicate
  inner variable declarations that would then shadow already-declared outer variables

## Indentation Convention
- Link.cpp uses tabs (not spaces) for indentation
- When editing with Python string replacement, use `\t` explicitly
- Edit tool fails on tab-indented files — use Python script replacements instead

## Known Bugs (from project memory)
- `(5/4)` integer division in `channel.cpp:1246,713,763` → should be `5.0/4.0`
- Double shadow fading O2I in `Link.cpp:845`
- `d_2D_in = min(U1,U2)` is CORRECT per 3GPP TR 38.901
- SV dB: amplitude → use `20*log10` not `10*log10`

## Refactoring Patterns
- Always use Python for exact string replacement in tab-indented C++ files
- Read exact line repr() first to confirm tab vs space indentation before editing
- After removing a block, check for leftover double blank lines and clean them
- Verify with grep that no references remain after removal
