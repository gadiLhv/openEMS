# Mode-file sampling: node point-values vs Yee edge integrals

Record of the benchmark behind commit `f010ec3`, kept because the result is not
obvious and the failure mode is silent.

## What is being compared

openEMS's DOFs are line integrals -- `Volt_amp` is `int(E.dl)` on a primal edge,
`Curr_amp` is `int(H.dl)` on a dual edge -- but mode files have historically been
written by point-evaluating the FEM solution at grid nodes. The two agree to
`O(h^2 d2E)` for a smooth mode and diverge to `O(1)` at a re-entrant metal wedge,
where `E ~ rho^-1/3` is unbounded and the tangential component is two-valued:
zero along the conductor, singular through the gap. One nodal number cannot
carry both, so the solver's point evaluation returns whichever element it lands
in, and those corner samples end up the largest values in the file.

Three arms per case, all at identical settings:

- `stock`  -- the committed mode file (rect and coax ship 2x oversampled;
              see "Open" below, the oversampling is deliberate and load-bearing)
- `point`  -- node point-sampling on the simulation's own port grid
- `edge`   -- Gauss-Legendre integration along each Yee edge, contracted to the
              edge average, then collapsed onto the primal nodes the way
              `Engine_Interface_FDTD`'s `NODE_INTERPOLATE` does
              (`engine_interface_fdtd.cpp:125-149` for E, `:213-231` for H)

`point` vs `edge` is the clean A/B: same `.geo`, same solve, same normalisation,
same grid -- only the sampling differs. `stock` is the reference.

RectWG is a null *by construction*, not merely by measurement: `Ey = sin(pi x/a)`
is constant along y, so its Yee edge average equals its point value and the two E
files are bit-identical. Only `Hx` differs, by a cell-average of the sine
(1.000000 -> 0.999286). If rect moves, something other than the sampling moved it.

## Results

15,000 iterations -- nothing converges, so every return loss is VOID:

    case  arm       steps   conv?  final E dB   dB/1k TS   |Gam| dB      resid   balance
    rect  stock     15000      no      -40.95     -0.872       VOID          -   0.996-1.014
    rect  point     15000      no      -39.76     -0.830       VOID          -   0.996-1.014
    rect  edge      15000      no      -39.77     -0.967       VOID          -   0.998-1.017
    coax  stock     15000      no      -31.83     -4.746       VOID   2.78e-01   0.875-92.349
    coax  point     15000      no      -18.72     -2.621       VOID   2.73e-01   0.860-90.030
    coax  edge      15000      no      -23.11     -3.273       VOID   2.39e-01   0.789-90.125
    cpw   point     15000      no       -0.00     +1.357       VOID   6.47e-02   0.352-1.175
    cpw   edge      15000      no       -0.00      0.000       VOID   6.04e-02   0.384-1.154

30,000 iterations:

    case  arm       steps   conv?  final E dB   dB/1k TS   |Gam| dB      resid   balance
    rect  stock     16305     YES      -50.46     -1.509     -41.88          -   0.997-1.015
    rect  point     21648     YES      -53.64     -0.840     -43.06          -   0.999-1.014
    rect  edge      17463     YES      -51.39     -1.263     -42.24          -   0.999-1.017
    coax  stock     30000      no      -30.49     -0.999       VOID   2.25e-01   0.790-1.356
    coax  point     30000      no      -33.01     -1.068       VOID   1.90e-01   0.804-1.295
    coax  edge      26122     YES      -67.67     -3.356     -16.73   1.79e-01   0.935-1.302
    cpw   point     30000      no      -44.87     -3.281       VOID   1.55e-03   1.169-1.251
    cpw   edge      28490     YES      -52.64     -3.682     -21.52   1.32e-03   1.129-1.173

`steps` is the true iteration count (`Time for N iterations`), not the last
printed sample. `conv?` is the absence of openEMS's max-timesteps warning.

## What it says

- **`edge` is the only arm that converges** on both coax and CPW within 30k.
  Coax reaches -60 dB at 26122 steps and lands -16.73 dB; CPW at 28490 steps
  lands -21.52 dB, reproducing the -21.50 dB of the full-length run.
- **A fixed-step energy reading is not a convergence measurement.** At 15k the
  coax `stock` arm looks best (-31.83 vs -23.11) -- it was caught mid-transient.
  Over 30k the decay *rates* reverse it: `edge` sustains -3.36 dB/1k while
  `stock` and `point` stall at ~-1.0.
- **Residual and energy balance improve monotonically** wherever the scheme is
  applied: coax resid 2.25e-01 -> 1.79e-01, balance 0.790-1.356 -> 0.935-1.302;
  CPW 1.55e-03 -> 1.32e-03 and 1.169-1.251 -> 1.129-1.173.
- **Rect stays null.** All three converge within 1.2 dB, and `point` is nominally
  best on |Gamma| while being slowest to converge -- i.e. no consistent winner,
  the correct answer for a mode with no wedge.
- The coax gains ~2.5 dB (-14.18 -> -16.73) on a geometry with **no re-entrant
  corner at all**, which was not expected.

## Open

- **Dual-grid sampling.** The `stock` rect and coax files are 2x oversampled by
  design, because openEMS also samples on the dual grid: `ProcessModeMatch` takes
  its disc lines with `dualMesh` (`processmodematch.cpp:185-189`) and the
  excitation samples `E_n` at `GetYeeCoords(..., false)`, dual along the
  component's own axis. A primal-only file forces openEMS to average two stored
  samples to answer those, which smooths the very thing the file specifies.
  Every `edge` arm above paid that penalty and still won. The fix is to emit on
  the interleaved primal-union-dual grid, storing the edge integral at the edge
  midpoint and the node-collapsed value at the primal node, so every consumer
  lands on a stored sample. Untested.
- `ngauss = 5` is unverified for convergence; regenerate at 9 and diff.
- The coax return loss is unimpressive in absolute terms. Likely post-processing:
  no broadband behaviour was put into the TEM cases, which matters at the low end.
- `MODE='MUR'` diverges on the CPW geometry for old and new mode files alike --
  separate, pre-existing, untouched here.

## Reproducing

The generator lives outside this repo, in `blit_port_mode_solver/`:
`yee_mode_export.py` (the sampling module) and `gen_cpw_modes_yee.py`, whose
`--point` flag regenerates the node-sampled baseline for A/B.

The benchmark used copies of the three test scripts with `NrTS` overridden to
15000 and 30000, swapping mode files between arms. CPW must be run with
`CPW_ABS_MODE=SCALAR`. Note a mode file is tied to its mesh: the pre-`f010ec3`
`CPW_E.csv` is 28x19 and does not align with the 25x19 the test now builds, so
`git checkout` is not a valid baseline -- regenerate with `--point` instead.
