# DFT calculation

## Convergenge w.r.t. cutoff energy

```bash
mpirun -n 3 abinit conv-Ecut.abi 1> out-Ecut.log 2> err.log
```

![Convergence plot: Etot vs Ecut](./fig/conv-Ecut.svg)

Fitting an exponential we get:

$$ y = \exp(12.89144 - 0.585344x) - 62.282 \quad [\mathrm{Ha}] $$

Taking $E_{cut} = 40.5\ \textrm{Ha}$, we get an error of ca. $20\ \mathrm{\mu Ha}$

## Convergence w.r.t. $k$-points and $t_{smear}$

```bash
mpirun -n 3 abinit conv-kpt-smear.abi 1> out-kt.log 2> err.log
```
![Convergence plot: Etot vs tsmear for different ngkpt](./fig/conv-kpt-smear.svg)

From plot, best option with reasonable amount of $k$-points is:
```bash
ngkpt 7 7 7
tsmear 0.006 # 6mHa
```
