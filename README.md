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

1. $\rho (T)$ via VDP technique. MORE DETAILED IN ./R-vs-T/readme.md
    - Iterative solution of the VDP equation from $R_{front}$ and $R_{back}$.
        - Filter data (we can try both methods):
        - Discard high $\Delta$T for $R_{front}$ and $R_{back}$.
        - Discard $R_{front}$ and $R_{back}$ corresponding to high T fluctuation.
    - Plot $\rho (T)$.
        - Merge proper datasets to obtain complete curve from 300k to 2.9k.
2. Condcutivity $\sigma (T)$
    - Through $n$ and $m^{*}$ from DFT.
    - Through $v_f$ and DOS from DFT.
        - In both cases one parameter is exact (n and DOS) and the other is approximated ($m^{*}$ and $v_f$).
    - Find $\tau$.
        - Not easy to compare $\tau$ with literature, so we find and discuss $l$ instead.
    - Ioffe-Regel parameter calculation: we want it greater than 1 amap.
        - Calculation explicitely requested for low T above Tc, suggested also for higher T.

### Resistivity vs Temperature

1. RRR calculation
    - Discussion of RRR wrt thickness and mean free path. Comparison with literature.
2. Bloch Gruneisen Model FIT
    - Fundamental parameters: $\rho_0$, n, $\Theta_{BG}$.
        - We can let them free, but it's not easy to obtain proper fit.
        - We can set $\rho_0$ since we have it.
        - A way to determine n is to plot $ln(\rho - \rho_0)$ vs $ln(T)$. This plot will hae initial angular coefficient equal to n! But $\rho_0$ must be really precise.
    - The fit is always linear for high T.
    - Great focus on the scattering crossover (where the curve changes power)
        - If from the fit we get n=5 there are no big problems, material is quite good.
        - If we obtain something smaller than 5, we need to distinguish among different scattering mechanism, identifying which is the dominant one.
        - If $n \sim 2$, the material is strongly correlated.
    - From the A coefficient we determine $\lambda_{TR}$, to compare with DFT and literature.
