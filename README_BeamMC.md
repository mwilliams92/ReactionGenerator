# BeamMC

Standalone C++11 Monte Carlo for beam energy loss, straggling, approximate multiple scattering, and two-body reaction recoils.

## Build and run

```bash
g++ -O2 -std=c++11 BeamMC.cpp -o BeamMC
./BeamMC beam
./BeamMC reaction
```

The executable reads `beam_input.txt`, `target_input.dat`, `degrader_input.txt`, and, in reaction mode, `input_reaction.txt`.

## Stopping tables

The code uses one SRIM stopping table for the complete material rather than separate tables for each element. This is appropriate for compound-target SRIM calculations such as `Ga64_in_CH2.srim`.

A material file contains:

```text
thickness_mm = 1.0
density_g_cm3 = 0.95
effective_Z = 6
effective_A = 12
stopping_table = Ga64_in_CH2.srim
srim_energy_unit = keV
srim_stopping_unit = MeV_cm2_g
```

`effective_Z` and `effective_A` are used only for approximate multiple scattering. The stopping table must be generated for the actual projectile ion in the complete compound material. The parser takes the first numeric column as energy and the last numeric column as stopping power; headers are ignored.

## Reaction stopping tables

Reaction mode reads these additional entries from `input_reaction.txt`:

```text
beam_srim_target = Ga64_in_CH2.srim
recoil_srim_target = Zn61_in_CH2.srim
recoil_srim_degrader = Zn61_in_degrader.srim
```

`beam_srim_target` is used for the incoming beam before the reaction. `recoil_srim_target` is used for the recoil while it traverses the remaining target. `recoil_srim_degrader` is used after the target. The degrader is skipped when `enabled = 0`; its reaction-specific table is then not needed.

## Input files

`beam_input.txt` specifies beam species, energy, energy spread, emittance, spot size, particle count, step size, random seed, output file, and straggling coefficient.

`target_input.dat` specifies target thickness, density, effective scattering parameters, and the beam stopping table.

`degrader_input.txt` specifies whether the degrader is enabled, its thickness and density, effective scattering parameters, and its stopping table for ordinary beam mode. In reaction mode, `recoil_srim_degrader` overrides this stopping table for the recoil.

`input_reaction.txt` specifies beam, target, ejectile, and recoil Z/A values, recoil excitation energy, Q value, and reaction-specific stopping tables.

## Output

Energies are in MeV, positions in mm, and `theta_x_mrad`/`theta_y_mrad` are horizontal and vertical angles in mrad relative to the nominal beam axis.

For each reaction event, the reaction depth is sampled uniformly between the target entrance and exit. The incoming beam is transported through that depth using `beam_srim_target`; the two-body recoil is then generated at that location and transported through the remaining target thickness using `recoil_srim_target`, followed by the optional degrader using `recoil_srim_degrader`. Thus, recoil energy loss in the residual target is included before degrader transport. The model samples an isotropic centre-of-mass angle and uses non-relativistic two-body kinematics.

This is a baseline model using linear stopping-power interpolation, Gaussian straggling, and approximate Highland multiple scattering. Validate it against SRIM/TRIM or experimental data before production use.
