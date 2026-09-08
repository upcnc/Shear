# Shear — Wang-style IBM biofilm detachment

Discrete-unit immersed-boundary model of biofilm **partial tearing and wash-off** in a straight microchannel, after

> Jiankun Wang, Yumeng Fu, Jin Wu, Jin Li, Xiaoling Wang,  
> Hydrodynamic modelling of biofilm detachment via the Immersed Boundary Method,  
> *International Biodeterioration & Biodegradation* 207 (2026) 106245.

This repository is the shear / detachment test bed. Growth, pore-throat FSI and honeycomb transport live in `upcnc/Biofilm`.

## What the paper does, and what this code keeps

| Paper | This code |
|---|---|
| IB2d + circular units of 12 Lagrangian points | D2Q9 LBM + Peskin 4-point IBM, 8-point rings |
| Unit–unit springs, form if d < Tc, break if |L-L0| > Tf = 0.6 r | same rule |
| Bottom units tethered to the wall | three wall springs per basal unit |
| Cases: semicircle layered / random, rectangle layered | case 1 / 2 / 3 |
| Critical modulus: soft films slough, stiff films only deform | sweep E |
| Grid h = 0.1 um | coarsened to 2 um so it runs on a laptop |

You should see:

- low E (e.g. 5e5-1e6 Pa): inter-unit bonds break, the top lifts, basal units stay — partial tear, not one-piece drop;
- high E (e.g. 4e6 Pa): the body leans downstream and saturates;
- case 1 tears in layers; case 2 tears more irregularly; case 3 loses a large mid-block and leaves a residual carpet.

## Build

```bash
g++ -O3 -std=c++17 -fopenmp wang_ibm_detach.cpp -o wang_detach
```

Windows MinGW: add `-lstdc++fs` if filesystem is missing. The source defines M_PI.

## Run

```bash
./wang_detach result_case1_soft  1  1.0e6
./wang_detach result_case1_stiff 1  4.0e6
./wang_detach result_case2_soft  2  1.0e6
./wang_detach result_case3_soft  3  1.0e6
```

Arguments: outdir, case (1/2/3), E in Pa.

## Output

- history.csv — living inter-unit bonds, still-anchored units, detached fraction
- state_XXXXXX.dat — Tecplot field + unit centres + living bonds
- summary.txt

## Parameters you should sweep first

1. E = 5e5, 1e6, 2e6, 4e6 on case 1 — find the jump from continuous wash-off to saturated lean (paper Ec ~ 3.5e6 Pa on their grid).
2. Same E on cases 1-3 — morphology changes the critical modulus.
3. Leave Tf = 0.6 r fixed; that is their experimental calibration.

Maxwell viscoelastic rest-length evolution is not in this first commit.
