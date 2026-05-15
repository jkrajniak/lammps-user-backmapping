# Theoretical Background

This page explains the backmapping method implemented in the package.

!!! cite "References"
    Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of
    Polymers from Coarse-Grained to Atomistic Descriptions",
    *J. Chem. Theory Comput.* 2016, 12, 5549--5562.
    [DOI: 10.1021/acs.jctc.6b00595](https://doi.org/10.1021/acs.jctc.6b00595)

    Krajniak, Zhang, Pandiyan, Nies, Samaey, "Reverse Mapping Method for Complex
    Polymer Systems", *J. Comput. Chem.* **2018**, *39*, 648--664 (first
    published online **6 December 2017**). Extends the adaptive-resolution
    reverse mapping idea to **topologies whose connectivity is not fixed
    a priori**, including **epoxy** and **trimethylol melamine** networks,
    a **hyperbranched** polymer, and **PET** from step-growth chemistry, using
    the same single control-parameter construction as in the linear-polymer
    formulation.
    [DOI: 10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)

## Overview

Backmapping (reverse mapping) is the process of reintroducing atomistic detail
into a coarse-grained simulation. The method implemented here uses a
time-dependent, smooth transition controlled by a resolution parameter
**&lambda;** that is ramped from 0 (fully CG) to 1 (fully AT) during the
simulation.

For **linear polymer melts**, the underlying formulation is described in the
**JCTC 2016** reference above. **Polymer networks and other complex
architectures** whose CG connectivity arises from **chemical reactions** or
step-growth processes were treated with the same reverse-mapping construction
in the **J. Comput. Chem.** work (online **2017**, volume **2018**). This
repository focuses on delivering that methodology through **LAMMPS** and
**`backmap-prep`** so users are not tied to legacy **ESPResSo++** bonded
implementations; fully worked **network** examples are a natural extension of
the current **melt** tutorials.

Unlike instantaneous mapping approaches that place atoms geometrically and then
energy-minimize, this method gradually introduces atomistic interactions while
fading out coarse-grained ones. The system relaxes dynamically, avoiding the
need for aggressive minimization steps.

## The Lambda Parameter

Each atom carries a per-atom resolution parameter &lambda;:

- **&lambda; = 0**: the atom is at purely coarse-grained resolution
- **&lambda; = 1**: the atom is at fully atomistic resolution
- **0 &lt; &lambda; &lt; 1**: the atom is in the transition region

Lambda increases linearly at each timestep:

\[
\lambda(t + \Delta t) = \min\bigl(\lambda(t) + \alpha,\ 1\bigr)
\]

where **&alpha;** is the ramp rate. A smaller &alpha; produces a slower,
gentler transition. Typical values are 10<sup>-4</sup> to 10<sup>-3</sup>.

### Uniform vs Non-uniform Lambda

In **uniform** mode (default), all atoms start at the same initial &lambda;
value (&lambda;<sub>0</sub>) and increase together.

In **non-uniform** mode, each atom receives a random offset to its initial
lambda, creating a staggered transition where different parts of the system
reach full AT resolution at different times. This can help reduce artifacts
from a sudden global transition.

## Force Weighting

All interactions (pair, bond, angle) are weighted by lambda. The weighting
scheme uses the product of the lambda values of the interacting atoms.

### Non-bonded Pair Interactions

For a pair of atoms *i* and *j*:

- **AT interactions**: weighted by \( w_\text{AT} = \lambda_i \times \lambda_j \)
- **CG interactions**: weighted by \( w_\text{CG} = 1 - \lambda_i \times \lambda_j \)

This ensures:

- At &lambda; = 0: only CG interactions are active (\( w_\text{CG} = 1 \), \( w_\text{AT} = 0 \))
- At &lambda; = 1: only AT interactions are active (\( w_\text{CG} = 0 \), \( w_\text{AT} = 1 \))
- During the transition: both sets of interactions contribute with complementary weights

### Bonded Interactions (Bonds and Angles)

Cross-CG bonded interactions (bonds and angles that span CG bead boundaries)
use the same weighting scheme:

- **AT cross bonds/angles**: \( w = \lambda_i \times \lambda_j \)
- **CG cross bonds/angles**: \( w = 1 - \lambda_i \times \lambda_j \)

For angles, the weight is computed from the first and last atoms of the
angle triplet (*i*-*j*-*k*), using \( \lambda_i \) and \( \lambda_k \).

Intra-bead AT bonds and angles (within a single CG bead) use standard LAMMPS
styles without lambda weighting, since they exist at both resolutions.

## Cross Interactions

Cross interactions are bonded terms that span CG bead boundaries. They are
essential for maintaining structural integrity during the transition.

### CG Cross Bonds

These are the bonds that connect neighboring CG beads. They are typically
derived from iterative Boltzmann inversion (IBI) or other coarse-graining
methods, and are provided as tabulated potentials. They fade out as &lambda;
increases:

\[
F_\text{CG bond} = (1 - \lambda_i \lambda_j) \times F_\text{table}(r)
\]

### AT Cross Bonds

These are the atomistic bonds between atoms in different CG beads. They fade
in as &lambda; increases:

\[
F_\text{AT bond} = \lambda_i \lambda_j \times (-k)(r - r_0)
\]

### AT Cross Angles

Same weighting as AT cross bonds, applied to angle potentials between atoms
in different CG beads.

## CG Bead Management

During the backmapping simulation, each CG bead coexists with the AT atoms
it represents. The `fix backmap` style manages the relationship:

1. **Position tracking**: After each timestep, the CG bead position is updated
   to the center-of-mass (COM) of its constituent AT atoms.

2. **Force distribution**: Forces on CG beads (from CG pair interactions and
   CG cross bonds) are redistributed to AT atoms proportional to their mass
   fraction:

    \[
    \mathbf{F}_i^\text{AT} \mathrel{+}= \frac{m_i}{M_\text{CG}} \mathbf{F}_\text{CG}
    \]

    where \( m_i \) is the mass of AT atom *i* and \( M_\text{CG} \) is the
    total mass of all AT atoms in the bead.

3. **Velocity zeroing**: CG bead velocities are set to zero at each step to
   prevent them from drifting independently of their AT atoms.

## Simulation Protocol

When AT atoms are first placed inside CG beads, their positions are
approximate and may overlap. A robust multi-phase protocol prevents the
large forces from these overlaps from destabilising the simulation
("Bond atoms missing" errors). The protocol has been validated on
dodecane and polyethylene systems from 180 to 3 500 atoms.

### Phase 0a: Energy Minimisation

CG beads are frozen (`fix setforce 0.0 0.0 0.0`) and the system is
energy-minimised to relieve the worst AT overlaps without moving CG
positions.

### Phase 0b: AT Relaxation

AT atoms are integrated with `fix nve/limit` (max displacement
&asymp; 0.01 &Aring;/step) at a very small timestep (e.g. 0.01 fs) for
&sim;10 000 steps. CG beads remain frozen. This gently relaxes remaining
steric clashes.

### Phase 1: Lambda Ramp

All atoms are integrated with `fix nve/limit` (max displacement
&asymp; 0.05 &Aring;/step) at a small timestep (e.g. 0.10 fs). The
lambda ramp is activated (`fix_modify bm active yes`). Lambda increases
by &alpha; each step; CG interactions fade out while AT interactions
fade in. Using `nve/limit` instead of full NVE prevents force spikes
from launching atoms out of the box.

### Phase 2: NVT Equilibration

Once &lambda; reaches 1.0, the system is fully atomistic. The integrator
is switched to NVT and the timestep is ramped gradually:

| Segment | Timestep | Steps | Purpose |
|---------|----------|-------|---------|
| 1 | 0.25 fs | 10 000 | Gentle thermalisation |
| 2 | 0.50 fs | 5 000 | Intermediate relaxation |
| 3 | 1.00 fs | 5 000 | Production timestep |

The gradual ramp avoids sudden energy injection that can blow up a
freshly backmapped system.

!!! tip "Small systems"
    For small test systems (&lt; 500 atoms), the simple protocol of
    NVE + Langevin at `dt = 1.0 fs` may work. The multi-phase protocol
    is essential for production-scale systems (&gt; 1 000 atoms).
