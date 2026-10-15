# Part A2: MATLAB simulation and visualization

Run from this directory, as the solver expects relative paths:

```matlab
cd('/Users/giannis/Documents/repos/comp-em-fdtd/src/partA')
a2_a
```

The assignment's A2(a) asks for **dt > dt_c**. The default is now
`courant_multiplier = 1.1`, deliberately unstable. Your original value `0.9`
is the stable A2(b) experiment. Edit that single setting to compare them;
the geometry stays fixed. `make_video = false` skips video export.

## Review of your implementation

The three update functions correctly implement the vacuum TM curls:

```text
dHx/dt = -(1/mu0) dEz/dy
dHy/dt = +(1/mu0) dEz/dx
dEz/dt = (1/eps0) (dHy/dx - dHx/dy)
```

The Yee arrays are indexed `[x,y]`: Ez has `(nx+1,ny+1)` nodes, Hx has
`(nx+1,ny)` nodes at half y cells, and Hy has `(nx,ny+1)` nodes at half x
cells. At entry to iteration n, Ez is at `(n-1)*dt` and H is at
`(n-3/2)*dt`. Updating H first advances it to `(n-1/2)*dt`; the interior
electric update then advances Ez to `n*dt`. Replacing the source using
`times(n+1)` is therefore correct. Initially H is zero at `-dt/2` and the
electric source is imposed at t=0.

Updating only the electric interior and excluding the source's y endpoints
preserves all four zero PEC edges. The source uses amplitude 1 V/m and the
convention `Re{phasor * exp(+j*omega*t)}`. Its cosine turns on abruptly,
so startup transients accompany the intended 3.75 GHz carrier.

The source position needed correction: index i represents `(i-1)*dx`,
so `nx/2 + 1 = 31`, rather than 30, is the electric node at `x=a/2=3*lambda0`.
The grid has 60 x 7 cells; dx = 7.99447 mm, dy = 7.13792 mm, and
dt_c = 17.7604 ps. Both spacings appear correctly in the Courant bound.
Integer nx is explicitly rounded, and dx is computed from the exact length.

## Output and interpretation

`plot_a2.m` saves three figures under the project's ignored `output/partA/`:

- `*_amplitude.png`: max|Ez| on a logarithmic axis, plus an off-source probe.
- `*_snapshots.png`: signed 2D fields near 0.5, 2, 5, and 12 source periods.
  Each panel has its own physical color scale, labelled in V/m.
- `*_spacetime.png`: an early-time and a full-run x-time diagram. Each panel
  has its own fixed physical scale. Late unstable growth can overwhelm the
  full-run panel; use the early panel to inspect startup propagation.

The x-time diagram samples the nearest existing y node to b/2, and labels
its actual y/b. With seven y cells there is no electric node exactly at b/2.
Field images transpose the stored `[x,y]` slice so x is horizontal and y is
vertical. The dotted vertical line marks the source.

`animate_a2.m` writes `*.mp4` with **every** `ez_tracker` row as a frame,
including t=0. At 20 fps, playback duration is `size(ez_tracker,1)/20`
seconds; this is slowed playback of nanosecond simulation time.
For dt/dt_c > 1, colors show `Ez/max|Ez|` **for each frame**. The frame title
and lower logarithmic plot show the actual amplitude in V/m, so normalized
colors must not be read as evidence of bounded fields. Stable runs use a
single fixed physical scale in V/m throughout the video.

The helpers can also be called directly on an existing history without
rerunning the solver, using the calls at the end of `a2_a.m`.

The video renderer was tested in MATLAB R2025b on macOS. It exports fixed
1024 x 656 PNG frames before encoding; dimensions divisible by 16 avoid an
MPEG-4 encoder artifact observed on this installation. No additional toolbox
or external video encoder is required.

## Validation performed in MATLAB

Both 12-period cases were executed from `src/partA`. Independent linear-field
tests confirmed the three analytic curl increments on unequal dx/dy, and
the electric update preserved PEC edges. Full histories had finite fields,
exact source values, and exactly zero left/right symmetry error.

| dt/dt_c | Updated steps | Stored frames | Peak absolute Ez [V/m] |
| --- | ---: | ---: | ---: |
| 0.9 | 201 | 202 | 2.33124 |
| 1.1 | 164 | 165 | 1.90919e39 |

The stable history's residual from the prescribed transverse sine profile
was 2.0e-15 V/m. The large unstable amplitude is expected for A2(a);
its precise value and onset depend on roundoff seeding, and are not an
analytical physical prediction. Video frame counts were checked by decoding
the completed files, and representative frames were visually inspected.

## Your interpretation tasks

1. Use the growth curve and late snapshots to distinguish numerical
   instability from physical interference or resonant buildup in a lossless
   cavity. What spatial structure accompanies the large amplitude?
2. Starting from your three discrete updates, insert a Fourier mode and find
   the condition for a real temporal frequency. Explain why the worst-case
   spatial mode gives the bound containing both dx and dy, and how tiny
   numerical perturbations can excite modes absent from the intended source.
3. Compare S=0.9 and S=1.1 at equal physical times. Identify outgoing and
   reflected traces, and discuss why an abrupt turn-on can produce arrivals
   before the carrier group-delay estimate. Use the actual color scales
   when comparing amplitudes.
