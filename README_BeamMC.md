# BeamMC

BeamMC is a standalone C++11 Monte Carlo simulation for transporting an ion beam through a target and optional degrader. It calculates energy loss, Bohr energy straggling, approximate multiple scattering, and—optionally—two-body reaction recoils.

The code is intended for compact target/degrader studies where SRIM stopping tables are available for the specific projectile and complete material being traversed.

## Build and run

### Beam mode

```bash
g++ -O2 -std=c++11 BeamMC.cpp -o BeamMC
./BeamMC beam
```

### Reaction mode

```bash
./BeamMC reaction
```

The program reads:

- `beam_input.txt`
- `target_input.dat`
- `degrader_input.txt`
- `input_reaction.txt` in reaction mode

All files should normally be in the same directory as the executable. SRIM files may be kept in a separate directory; specify the relative path directly in the relevant input file, for example:

```text
STOPPING_TABLE              SRIM_Files/Strontium_in_SiHe.txt
```

The directory must exist relative to the directory from which `BeamMC` is run.

## ROOT output

ROOT output is optional. Compile with ROOT support using:

```bash
g++ -O2 -std=c++11 -DBEAMMC_WITH_ROOT BeamMC.cpp -o BeamMC \\
    $(root-config --cflags --libs)
```

Then enable it in `beam_input.txt` if supported by the current source version:

```text
ROOT_OUTPUT                 1
ROOT_FILE                   BeamMC.root
```

The ROOT file contains outgoing-energy histograms, a two-dimensional vertical-versus-horizontal angle histogram, and a `TTree` with incident and outgoing particle quantities. The CSV output remains available independently of ROOT output.

If ROOT output is requested, the executable must have been compiled with `-DBEAMMC_WITH_ROOT` and the ROOT compiler/linker flags.

## Input-file format

Input variables use uppercase names with underscores. The value follows the variable name after spaces or tabs. Comments begin with `#` and blank lines are ignored.

For example:

```text
BEAM_ENERGY                279.0
BEAM_ENERGY_UNIT           MeV
ENERGY_SPREAD              0.0
ENERGY_SPREAD_MODE         sigma_abs
```

The parser also accepts `=` between a variable and its value in versions that retain the older key/value syntax, but the aligned whitespace format above is recommended.

## `beam_input.txt`

This file contains the incident beam and simulation settings. A typical file is:

```text
# ============================================================================
# BEAM SPECIES AND ENERGY
# ============================================================================
BEAM_Z                     38
BEAM_A                     93
BEAM_ENERGY                279.0
BEAM_ENERGY_UNIT           MeV
ENERGY_SPREAD              0.1
ENERGY_SPREAD_MODE         sigma_percent

# ============================================================================
# TRANSVERSE BEAM PARAMETERS
# ============================================================================
EMITTANCE_X_MM_MRAD        2.0
EMITTANCE_Y_MM_MRAD        2.0
BEAM_SIGMA_X_MM            1.0
BEAM_SIGMA_Y_MM            1.0

# ============================================================================
# MONTE CARLO AND TRANSPORT
# ============================================================================
PARTICLES                  100000
STEP_UM                    0.1
BOHR_STRAGGLING_SCALE      1.0
SEED                       12345
OUTPUT                     beam_output.csv

# ============================================================================
# OPTIONAL ROOT OUTPUT
# ============================================================================
ROOT_OUTPUT                0
ROOT_FILE                  beam_output.root
```

`STEP_UM` is the transport step in micrometres. The default is `0.1` micrometres, appropriate for targets with thicknesses of order micrometres. Smaller steps can be used when greater resolution is needed.

### Beam energy units

`BEAM_ENERGY_UNIT` may be:

- `MeV`: `BEAM_ENERGY` is the total kinetic energy in MeV.
- `MeV/u`: `BEAM_ENERGY` is the energy per nucleon. The code converts it using

$$
E_{\mathrm{MeV}} = E_{\mathrm{MeV/u}} A,
$$

where $A$ is the beam mass number.

The absolute energy spread follows the same convention as the beam energy. Thus, for `MeV/u`, an absolute spread in `MeV/u` is multiplied by $A$ before transport.

### Energy-spread modes

`ENERGY_SPREAD_MODE` may be:

- `sigma_percent`: RMS spread as a percentage of the nominal beam energy.
- `fwhm_percent`: FWHM spread as a percentage of the nominal beam energy.
- `sigma_abs`: absolute RMS spread in the selected beam-energy unit.
- `fwhm_abs`: absolute FWHM spread in the selected beam-energy unit.

For FWHM input, the Gaussian RMS width is calculated as

$$
\sigma = \frac{\mathrm{FWHM}}{2\sqrt{2\ln 2}}.
$$

A zero spread produces a monoenergetic incident beam.

## `target_input.dat`

The target file specifies the geometry, density, approximate scattering material, and the SRIM table for the incoming beam in the complete target material.

Example:

```text
# ============================================================================
# TARGET GEOMETRY AND MATERIAL
# ============================================================================
THICKNESS_MM               0.001
DENSITY_G_CM3              1.5894
EFFECTIVE_Z                10
EFFECTIVE_A                20
STOPPING_TABLE              SRIM_Files/Strontium_in_SiHe.txt
SRIM_ENERGY_UNIT            auto
SRIM_STOPPING_UNIT          auto
```

`THICKNESS_MM` is the physical target thickness in mm. `DENSITY_G_CM3` is used for areal-density calculations, Bohr straggling, and approximate scattering. `EFFECTIVE_Z` and `EFFECTIVE_A` are effective values for the compound and are used only by the approximate multiple-scattering model; stopping power comes from the SRIM table.

For native SRIM output, use `auto` for both unit settings. The parser reads the SRIM header, identifies the data table, reads the ion energy, and sums the electronic and nuclear columns:

$$
\left(\frac{dE}{dx}\right)_{\mathrm{total}}
=
\left(\frac{dE}{dx}\right)_{\mathrm{electronic}}
+
\left(\frac{dE}{dx}\right)_{\mathrm{nuclear}}.
$$

The native SRIM format is recognised even when the header appears as `Ion dE/dx dE/dx Projected ...` rather than `Ion Energy ...`.

For the supplied `Strontium_in_SiHe.txt` file, SRIM reports `keV / micron`. For example, near 279 MeV, the electronic and nuclear stopping values are added and the result is interpreted as MeV per micrometre after conversion. A stopping value of approximately `7.5E3 keV/micron` therefore corresponds to approximately `7.5 MeV/micron`.

## `degrader_input.txt`

The degrader is optional. Set `ENABLED` to `0` to disable it.

```text
# ============================================================================
# DEGRADER
# ============================================================================
ENABLED                    0
THICKNESS_MM               0.001
DENSITY_G_CM3              1.0
EFFECTIVE_Z                6
EFFECTIVE_A                12
STOPPING_TABLE              SRIM_Files/Alpha_in_degrader.txt
SRIM_ENERGY_UNIT            auto
SRIM_STOPPING_UNIT          auto
```

In ordinary beam mode, `STOPPING_TABLE` is used for the transported beam in the degrader. In reaction mode, the recoil-specific table from `input_reaction.txt` is used instead when supplied.

## `input_reaction.txt`

Reaction mode performs a two-body reaction at a random depth in the target. The input file specifies the beam, target, ejectile, recoil, reaction energetics, and projectile-dependent SRIM tables.

For the current example reaction,

```text
93Sr(alpha,n)96Zr
```

the incoming beam is `93Sr`, the target is an alpha particle, the ejectile is a neutron, and the heavy recoil is `96Zr`:

```text
# ============================================================================
# REACTION PARTICLES
# ============================================================================
BEAM_Z                     38
BEAM_A                     93
TARGET_Z                    2
TARGET_A                    4
EJECTILE_Z                  0
EJECTILE_A                  1
RECOIL_Z                    40
RECOIL_A                    96

# ============================================================================
# REACTION ENERGETICS
# ============================================================================
RECOIL_EXCITATION_MEV      0.0
Q_VALUE_MEV                0.0

# ============================================================================
# REACTION-SPECIFIC SRIM TABLES
# ============================================================================
BEAM_SRIM_TARGET            SRIM_Files/Strontium_in_alpha_target.txt
RECOIL_SRIM_TARGET          SRIM_Files/Zirconium96_in_alpha_target.txt
RECOIL_SRIM_DEGRADER        SRIM_Files/Zirconium96_in_degrader.txt

# ============================================================================
# OUTPUT
# ============================================================================
OUTPUT                      reaction_output.csv
```

The numerical `Q_VALUE_MEV` must be supplied from the adopted nuclear masses or reaction data. It is not calculated from the Z/A values alone. `RECOIL_EXCITATION_MEV` reduces the available kinetic energy for the final state.

The table assignments are:

- `BEAM_SRIM_TARGET`: beam transport from the target entrance to the reaction point.
- `RECOIL_SRIM_TARGET`: heavy-recoil transport through the remaining target thickness.
- `RECOIL_SRIM_DEGRADER`: heavy-recoil transport through the degrader.

All referenced files may be placed in `SRIM_Files`, as shown above.

## Reaction-event sequence

For every reaction event, BeamMC:

1. Generates an incident beam particle with the requested energy spread and transverse phase-space distribution.
2. Samples the reaction depth uniformly between the target entrance and exit.
3. Transports the beam through that depth using `BEAM_SRIM_TARGET`.
4. Generates the two-body recoil using the transported beam energy and direction.
5. Transports the recoil through the remaining target thickness using `RECOIL_SRIM_TARGET`.
6. Transports the recoil through the optional degrader using `RECOIL_SRIM_DEGRADER`.
7. Writes the final recoil energy and angles.

For a forward-going recoil, the remaining target path is increased by the emission angle. A backward-going recoil travels towards the upstream target surface and does not subsequently pass through the downstream degrader.

The current two-body model is non-relativistic and samples an isotropic centre-of-mass emission direction. It should be validated against a dedicated kinematics calculation before precision use.

## Energy loss and straggling

SRIM supplies the mean stopping power. BeamMC interpolates the total stopping power linearly in energy and applies the deterministic loss at each transport step.

The current stochastic energy-loss model is Bohr straggling. The Bohr variance is calculated from the step areal density and the effective target composition, with the configured scale factor applied. In the standard convention:

- `BOHR_STRAGGLING_SCALE = 1.0` uses the full Bohr variance.
- `BOHR_STRAGGLING_SCALE = 0.01` scales the Bohr variance by 1%.
- `BOHR_STRAGGLING_SCALE = 0.0` disables stochastic straggling.

The transported kinetic energy is constrained to remain non-negative and passive material cannot increase the particle energy above its value before the step.

Approximate Highland multiple scattering is applied independently using the effective target parameters. The reported angles are horizontal and vertical projected angles relative to the nominal beam axis, in mrad.

## Output quantities

CSV and ROOT outputs use the following conventions:

- Energy: MeV.
- Position: mm.
- Horizontal angle: `theta_x_mrad`.
- Vertical angle: `theta_y_mrad`.
- Angles: mrad relative to the nominal beam axis.

Reaction output contains the final heavy-recoil quantities after the residual target and optional degrader. ROOT output additionally stores incident and outgoing energy and angular quantities in histograms and a `TTree`.

## Important checks

- Ensure every SRIM path is correct relative to the directory from which `BeamMC` is launched.
- Use a stopping table generated for the actual ion and complete material being traversed.
- In reaction mode, provide separate beam-target and recoil-target tables when their stopping powers differ.
- Confirm that `Q_VALUE_MEV` and recoil excitation energy are physically appropriate for the selected reaction.
- Check that the SRIM energy range covers the energies encountered during transport.
- Use `STEP_UM = 0.1` as a starting point for micron-scale targets and reduce it if convergence tests require it.
- Compare energy-loss and angular distributions with SRIM/TRIM or experimental data before using the simulation for final results.