# Nb-Characterization

## DFT
--Insert the name of the file associated to this phase of the characterization--
--Insert the required steps below and mark the ones that have been already completed--
--Insert the data extrapolated from the simulations--

## AFM
--Insert the name of the file associated to this phase of the characterization--
--Insert the required steps below and mark the ones that have been already completed--
--Insert the data extrapolated from the analysis--

## Resistivity
### Resistivity and Mean Free Path
1. $\rho (T)$ via VDP technique.
  - Iterative solution of the VDP equation from $R_{front}$ and $R_{back}$
  - Filter data (we can try both methods):
    - Discard $R_{front}$ and $R_{back}$ that are too different
    - Discard $R_{front}$ and $R_{back}$ corresponding to high T fluctuation
  - Plot $\rho (T)$
    - Merge proper datasets to obtain complete curve from 300k to 2.9k
2. Condcutivity $\sigma (T)$
  - Through $n$ and $m^{*}$ from DFT
  - Through $v_f$ and DOS from DFT
    - In both cases one parameter is exact (n and DOS) and the other is approximated $m^{*}$ and $v_f$
  - Find $\tau$.
    - Not easy to compare $\tau$ with literature, so we find and discuss $l$ instead
  - Ioffe-Regel parameter calculation: we want it greater than 1 amap.
    - Calculation explicitely requested for low T above Tc, suggested also for higher T.
### Resistivity vs Temperature
