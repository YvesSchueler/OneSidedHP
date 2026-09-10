# Adjusted one-sided Hodrick–Prescott filter

MATLAB, Python, and R implementations of the adjusted one-sided HP filter (HP-1s*) from

> Wolf, E., Mokinski, F., and Schüler, Y. (2026). On adjusting the one-sided Hodrick–Prescott filter. *Journal of Money, Credit and Banking*, 58, 919–931. [doi:10.1111/jmcb.13240](https://doi.org/10.1111/jmcb.13240)

Please cite the paper when using the code.

## What the filter does

The one-sided HP filter (HP-1s) runs the two-sided filter on an expanding sample and keeps the last observation, so it uses only past data and is suitable for real-time work. Applied with the standard smoothing parameter (e.g. λ = 1,600 for quarterly data), however, it does not extract the same cyclical frequencies as the two-sided filter, and it dampens the amplitude of the cycle.

The adjustment corrects both. Given the two-sided smoothing parameter λ₂, it uses a smaller smoothing parameter λ₁ inside the one-sided filter and rescales the extracted cycle by a factor κ, with both chosen so that the power transfer function of HP-1s matches that of HP-2s. For λ₂ = 1,600 this gives λ₁ ≈ 650 and κ ≈ 1.15. For other values of λ₂ (6.25 to 1,000,000), the code supplies λ₁ and κ from fitted polynomials, or recomputes them by direct optimization.

## Files

| file | language |
|---|---|
| `adj1s_hpfilter.m` | MATLAB |
| `adj1s_hpfilter.py` | Python (NumPy, SciPy) |
| `adj1s_hpfilter.R` | R (base) |

## Usage

Input in all three versions: `y`, a single time series without missing values; `lm_2`, the two-sided smoothing parameter.
Output: the adjusted cyclical component, and the parameters `lm_1` and `kappa` used.

**MATLAB**
```matlab
[ycycle_adj, lm_1, kappa] = adj1s_hpfilter(y, 1600);        % polynomial adjustment (fast)
[ycycle_adj, lm_1, kappa] = adj1s_hpfilter(y, 1600, 1);     % adjustment by optimization
```

**Python**
```python
from adj1s_hpfilter import adj1s_hpfilter
ycycle_adj, lm_1, kappa = adj1s_hpfilter(y, 1600)              # polynomial adjustment (fast)
ycycle_adj, lm_1, kappa = adj1s_hpfilter(y, 1600, opt=True)    # adjustment by optimization
```

**R**
```r
source("adj1s_hpfilter.R")
res <- adj1s_hpfilter(y, 1600)             # polynomial adjustment (fast)
res$ycycle_adj; res$lm_1; res$kappa
```

Options:

- `opt`: if true, λ₁ and κ are found by minimizing the distance between the power transfer functions of HP-1s and HP-2s (evaluated on a 1,000-observation filter weight vector); otherwise they are read from the fitted polynomials. The two routes give nearly identical values for λ₂ in the range covered by the polynomials. Optimization takes roughly a minute.
- `sample` (Python and R only): if true, the optimization evaluates the transfer function at the actual sample length instead of 1,000. Useful for large λ₂ combined with short samples (e.g. λ₂ = 400,000 and T < 100).

All three versions return a cyclical component of length T; the first two entries are zero, since the filter needs at least three observations. Defaults are `lm_2 = 1600` and `opt` off in all three versions.

## Contact

Yves Schüler, Deutsche Bundesbank — yves.schueler (at) bundesbank.de

The views expressed are those of the authors and do not necessarily coincide with the views of the Deutsche Bundesbank or the Eurosystem.
