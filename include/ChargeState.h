#ifndef BEAMMC_CHARGE_STATE_H
#define BEAMMC_CHARGE_STATE_H

#include <random>

namespace BeamMC {

struct ChargeStateResult {
    double meanQ;    // <q>
    double sigmaQ;   // σ_q (RMS width in charge units)
    int    sampledQ; // integer charge state (after sampling)
};

class ChargeStateCalculator {
public:
    ChargeStateCalculator() {}

    // Zp, Ap: projectile (beam or recoil)
    // energyMeV: kinetic energy (MeV)
    // targetZ: atomic number (or effective Z) of the medium at exit
    // reducedWidthW: Schiwietz & Grande reduced width w (user input)
    ChargeStateResult sample(int Zp, int Ap,
                             double energyMeV,
                             double targetZ,
                             double reducedWidthW,
                             std::mt19937_64& random) const;

private:
    ChargeStateResult computeMoments(int Zp, int Ap,
                                     double energyMeV,
                                     double targetZ,
                                     double reducedWidthW) const;
};

} // namespace BeamMC

#endif // BEAMMC_CHARGE_STATE_H