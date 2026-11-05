# AVBMC detailed-balance reference

Run `python avbmc_reference.py`. The script uses nine particles in a periodic
L=6 box. The ideal case has no interactions; the interacting case has a
hard core at 0.7 and an attractive square well of depth 1 from 0.7 to 1.25,
with kT=1. It saves `summary.json`. The old acceptance factors are also
implemented in their K=1 form (`n_neigh=1`, `nbrd=N_in`,
`nnonassoc=N_out`), so the ideal-gas binomial check can detect their bias.

The Fortran `attempt_avbmc_pair_moves` corresponds term by term to
`avbmc_step`: fixed-type target selection; current `near`/`far` counts; a
uniform source-particle choice; K position and uniform quaternion destination
trials; old state plus K-1 source trials; Boltzmann Rosenbluth selection; the
same IN/OUT volume, count, and direction factors; and a Metropolis draw.
The Fortran host `calc_single_particle_energy` supplies the interaction energy
that `particle_energy` supplies here. `update_single_particle_cell` and the
final checkerboard rebuild synchronize accepted positions and quaternions.
The reference square well is deliberately simpler than the production SSP/SB
potentials, so these tests validate the algorithm, not GPU potential parity.

See `summary.json` for the measured chi-square and KS p-values. Samples are
thinned and the first quarter of each run is discarded. GPU validation decks
are in the campaign validation v6-v10 directories.
