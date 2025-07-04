# Resistance vs Temperature

## Procedure

1. Define temperature fluctuation cutoff $\Delta T_{cutoff}$ e.g. $30 \mathrm{mK}$
2. Obtain final $R_s\ \mathrm{vs}\ T$ dataset; two options:
    - Probably, Gonnelli's original idea
        1. 1. Extract two datasets of $R\ \mathrm{vs}\ T$: one for front, and one for back
        2. For both datasets individually, remove points for which temperature fluctuation is more than cutoff
        3. Define a list of points from lowest to highest temperature
        4. For each sample temperature, compute $R_f$ and $R_b$ as the linear interpolation between the two points with the nearest temperature, taken from the respective datasets
        5. Solve Van der Pauw equation on these interpolated values
    - Filter + interpolation:
        1. Extract two datasets of $R\ \mathrm{vs}\ T$: one for front, and one for back
        2. For both datasets individually, remove points for which temperature fluctuation is more than cutoff
        3. Apply filter (e.g. moving average, temperature bin average) optionally with (simple linear) interpolation to obtain two filtered datasets with same temperature values
        4. Solve Van der Pauw equation for pair of front and back resistances
    - Dataline-wise VdP
        1. Select an temperature fluctuation cutoff $\Delta T_{inter} \leq \Delta T_{cutoff} $ between front and back temperature 
        2. Remove all rows for which $\Delta T_{f-f}, \Delta T_{b-b} > \Delta T_{cutoff}$ and $\Delta T_{f-b} > \Delta T_{inter}$
        3. For each row solve Van der Pauw equation, and take as $T$ the average between front and back temperatures
3. Perform least square fit of Bloch-Gruneisen model on resulting dataset; two options:
    - Direct fit of the model: 
    - Simplify the fit by:
        1. Performing linear fit of constant function on the first points of the dataset (where the resistance remains reasonable flat) to find $R_0$
        2. Perform fit on $\log(R - R_0)$


## Observations
- Using Gonnelli's original idea, we sample and interpolate noise which is a bit meaningless. One way to counteract this issue is to a have sample temperature almost as dense as the original values: in this way the original noise is reasonably? preserved in the derived sample, and thus the lsq fit parameters' uncertainties have direct connection with original noise.
- Using both filter and interpolation, we can have sample temperature much less dense than the original dataset (moreover they can be equispaced). Multiple choices of filtering can be adopted. However, we expect lsq fit parameters' uncertainty to be affected by the original noise only very slightly, since the contribution from the latter is greatly suppressed with initial filtering.
- Performing dataline-wide VdP preserves preserves noise on the resistance more realistically. However, it introduces some minor error on the temperature since front and back measurements are not performed at the same temperature. The only way to handle this error is to make sure that the error is fully covered by the uncertainty i.e. have a reasonable overlap of the uncertainty range of $T_f$ and $T_b$ (for each individual dataline)
- The linear region of resistance vs temperature includes a very big portion of points, which could make the flat and polynomial regions underrepresented (depends on relative density of points). The best solution may be truncating the temperatures used for the fit to a lower value.
- Performing the lsq fit of the logarithm is easier and likely has better convergence, but logarithm of zero or negative values will surely occur: how do we handle this case? Discarding these points may be detrimental as this tecnique would halve the datapoint at low temperature, reducing the "weight" of this region in the eyes of the fit. Moreover, it would skew the fit a little towards more positive values. 
