# DFT calculation

## Convergenge w.r.t. cutoff energy

```bash
mpirun -n 4 abinit conv-Ecut.abi 1> out-Ecut.log 2> err.log
```

![Convergence plot: Etot vs Ecut](./fig/conv-Ecut.svg)

Fitting an exponential we get:

$$ y = \exp(12.23956 - 0.5668441x) - 62.28018 \quad [\mathrm{Ha}] $$

Taking $E_{cut} = 42\ \textrm{Ha}$, we get an error of ca. $9.4\ \mathrm{\mu Ha}$

## Convergence w.r.t. $k$-points and $t_{smear}$

```bash
# Inside tmux or equivalent
mpirun -n 4 abinit conv-kpt-smear.abi 1> out-kt.log 2> err.log

# To follow iterations on other terminal
tail -f out-kt.log | grep -E "DATASET|ITER STEP"
```
![Convergence plot: Etot vs tsmear for different ngkpt](./fig/conv-kpt-smear.svg)

From plot, best option with reasonable amount of $k$-points is:
```bash
ngkpt 7 7 7
tsmear 0.006 # 6mHa
```

## Relaxation

```bash
# Lauch relaxation on tmux or equivalent
mpirun -n 4 abinit relaxation.abi 1> relaxation.log 2> err.log

# See progress
tail -f relaxation.log | grep -E "DATASET|ITER STEP"

# Get final lattice constant
cat relaxation.abo | grep "acell"
```

We perform BFGS structural relaxation to optimize cell size. We get the relaxed lattice constant of $6.2576431670\ \mathrm{Bohr} = 0.33114196111\ \mathrm{nm}$
which differs only slightly from the empirical value of $0.33004\ \mathrm{nm}$

## Electronic properties
```bash
mpirun -n 4 abinit electronic.abi 1> electronic.log 2> err.log
```
The relevant outputs are:
- Bands in `out/electronic_DS2_EBANDS.data`
- Density of states in `out/electronic_DS3_DOS`
- Fermi surface in `out/electronic_DS3_BXSF`

![Bands and density of states](./fig/ebands.svg)

## Phonon properties

1. Create derivative database for response function of phonons
    ```bash
    # Inside tmux or equivalent
    mpirun -n 4 abinit ddb.abi 1> ddb.log 2> err.log

    # To follow iterations on other terminal
    tail -f ddb.log | grep -E "== DATASET|ITER STEP|Perturbation"
    ```
2. Merge all the databases together:
    ```bash
    mrgddb < merge-ddb.txt
    ```
3. Compute phononic bands and DOS:
    ```bash
    anaddb phononic.abi
    ```
The relevant outputs are:
- Bands in `out/phononic_PHBANDS.data`
- Density of states in `out/phononic_PHDOS`

![Phonon bands and DOS](./fig/phbands.svg)