# BeamMC

Standalone C++11 Monte Carlo for energy loss, Gaussian energy straggling, approximate Highland multiple scattering, and optional two-body reaction recoils.

## Build and run

```bash
g++ -O2 -std=c++11 BeamMC.cpp -o BeamMC
./BeamMC beam
./BeamMC reaction
```

The executable reads these files from the current directory:

- `beam_input.txt`
- `target_input.dat`
- `degrader_input.txt`
- `input_reaction.txt` (only in `reaction` mode)

Outputs are CSV files. Beam-mode output is controlled by `output` in `beam_input.txt`; reaction mode writes `reaction_output.csv`.

## beam_input.txt

```text
Z = 31
A = 64
energy_MeV = 100.0
energy_spread_MeV = 0.5
emittance_x_mm_mrad = 2.0
emittance_y_mm_mrad = 2.0
beam_sigma_x_mm = 1.0
beam_sigma_y_mm = 1.0
particles = 100000
step_mm = 0.01
seed = 12345
output = beam_output.csv
straggling_coefficient = 0.02
```

Emittances are geometric rms emittances in mm mrad. The initial rms angular spreads are calculated as `emittance / beam_sigma`, converted to radians. `energy_spread_MeV` is a Gaussian sigma.

## target_input.dat

Fractions are interpreted as **mass fractions** and are normalized automatically.

```text
thickness_mm = 1.0
density_g_cm3 = 1.0
srim_energy_unit = keV
srim_stopping_unit = MeV_cm2_g
element = 1 1 1.00794 1.0 srim_H.dat
```

For a compound or mixture, add one `element = ...` line per element:

```text
element = 6 12 12.011 0.50 srim_C.dat
element = 8 16 15.999 0.50 srim_O.dat
```

Each SRIM file must contain numeric rows whose first numeric column is energy and whose last numeric column is stopping power. Header text is ignored. The default stopping-power unit is `MeV_cm2_g`; `MeV_cm2_mg` and `keV_um` are also accepted. Set the unit explicitly for your SRIM export.

## degrader_input.txt

Set `enabled = 0` to bypass the degrader. When enabled, use the same material format as the target:

```text
enabled = 1
thickness_mm = 0.05
density_g_cm3 = 2.0
srim_energy_unit = keV
srim_stopping_unit = MeV_cm2_g
element = 6 12 12.011 1.0 srim_C.dat
```

## input_reaction.txt

```text
beam_Z = 31
beam_A = 64
target_Z = 1
target_A = 1
ejectile_Z = 2
ejectile_A = 4
recoil_Z = 30
recoil_A = 61
recoil_excitation_MeV = 0.0
q_value_MeV = 0.0
```

`q_value_MeV` is optional but should be supplied for a physical reaction. The requested input species and excitation energy alone do not determine Q; nuclear mass excesses are required. The reaction model samples an isotropic centre-of-mass angle, calculates non-relativistic two-body kinematics, then transports the recoil through the remaining target and optional degrader.

## Output angles

`theta_x_mrad` and `theta_y_mrad` are horizontal and vertical angles relative to the nominal beam axis, in mrad. Positions are in mm and energies are in MeV.

## Physics scope

This is a practical baseline implementation, not a replacement for a validated transport code. It uses linear interpolation of SRIM stopping tables, a Gaussian straggling model, and the Highland multiple-scattering approximation. Validate against SRIM/TRIM or experimental data before using it for a production analysis. Reaction target recoil energy loss is simulated; energy loss of the ejectile is not needed for the recoil output and is therefore not generated.
