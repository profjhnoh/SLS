# Cfg Parameter Integration Memory

## Configuration File Parsing Pattern

### File Structure
- **cfg parsing**: `src/setSimulParam.cpp` using `Get_parameter(infile, "param_name", default_value)`
- **Variable definitions**: `h/Initiallization.h` (actual variables, no `extern`)
- **Extern declarations**: `h/common.h` (all `extern` declarations)
- **Scenario defaults**: `src/Initiallization.cpp` in `Set_Parameter()` function
- **Override location**: `src/Initiallization.cpp` after wavelength calculations (line ~522), before print statements (line ~523)

### Variable Declaration Pattern
1. Define in `h/Initiallization.h`: `Real variable_name;`
2. Declare extern in `h/common.h`: `extern Real variable_name;`
3. Parse in `src/setSimulParam.cpp` near line 127
4. Override in `src/Initiallization.cpp` after line 522

### Sentinel Value Pattern
- Use `-9999` as sentinel for "not set in cfg file"
- Check with `> -9998` to allow for floating point comparison

### Override Pattern Example
```cpp
if (cfg_BS_Tx_Power > -9998) {
    Total_BS_Tx_power = cfg_BS_Tx_Power;
    bs_maxpower = cfg_BS_Tx_Power;
    cout << "** BS Tx Power overridden by cfg: " << cfg_BS_Tx_Power << " dBm **" << endl;
}
```

### Important Dependencies
- When overriding noise figure, must also recalculate noise: `noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);`
- When overriding BS power, update both `Total_BS_Tx_power` and `bs_maxpower`

## Build System
- CMake-based build in `build/` directory
- Build command: `cd build && cmake .. && make -j8`
- Project compiles with warnings but warnings are expected
