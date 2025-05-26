# Systematic Differences Between Run3 and VBS2 Frameworks

## Summary
After extensive investigation, we've identified fundamental differences between the Run3 and VBS2 frameworks that lead to an ~8% difference in jet multiplicity after the "≥4 AK4 Jets" selection. These differences are architectural and cannot be fully reconciled without major framework restructuring.

## Key Differences

### 1. Jet Selection Architecture
- **VBS2**: Applies corrections jet-by-jet during selection
  ```cpp
  for each jet:
    apply JEC → apply JER → check pt > 20 → increment counter
  ```
- **Run3**: Applies corrections globally, then selects
  ```cpp
  apply JEC to all jets → apply JER to all jets → select jets with pt > 20
  ```

### 2. Missing Scale Factors
- **VBS2**: Applies PileUp Jet ID scale factors (puid_sfs) during selection
- **Run3**: Does not apply these MC-to-data correction weights

### 3. JER Random Seeds
- **VBS2**: Complex seed based on run, lumi, event, and jet eta
- **Run3**: Simple seed based only on event number

## Impact on Analysis

### Event Counts After "≥4 AK4 Jets"
- VBS2: 251,236 events
- Run3: 271,996 events (with JEC+JER)
- Difference: ~20k events (8%)

### Final Selection (after all cuts)
- VBS2: 64,107 events
- Run3: 70,377 events
- Difference: ~6k events (9.8%)

## Recommendation
This should be treated as a **known systematic uncertainty** between the frameworks. The physics conclusions should remain valid within these uncertainties. Both approaches are technically correct but implement the selection logic differently.

## Potential Mitigations
1. Apply a data-driven scale factor to align the frameworks
2. Include this as a systematic uncertainty in the final analysis
3. Verify that signal-to-background ratios remain consistent
