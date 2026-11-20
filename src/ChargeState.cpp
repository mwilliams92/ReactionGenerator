#include "ChargeState.h"

#include <algorithm>
#include <cmath>

namespace BeamMC {

namespace {
constexpr double AMU_MEV = 931.49410242;
constexpr double C_CM_PER_S    = 2.99792458e10;     // speed of light [cm/s]
constexpr double X_DENOM_SCALE = 3.6e8;             // from Shima definition

double safePow(double base, double exp)
{
    if (base <= 0.0) return 0.0;
    return std::pow(base, exp);
}

} // namespace

ChargeStateResult ChargeStateCalculator::computeMoments(
    int Zp, int Ap, double energyMeV, double targetZ, double reducedWidthW) const
{
    ChargeStateResult result;
    result.meanQ = 0.0;
    result.sigmaQ = 0.0;
    result.sampledQ = 0;

    if (Zp <= 0 || Ap <= 0 || energyMeV <= 0.0 || targetZ <= 0.0) {
        return result;
    }

	const double Z1 = static_cast<double>(Zp);
    const double Z2 = std::max(1.0, targetZ);

	// Relativistic velocity from kinetic energy E [MeV] and mass A*amu.
    const double massMeV = static_cast<double>(Ap) * AMU_MEV;
    const double gamma   = 1.0 + energyMeV / massMeV;
    const double beta2   = 1.0 - 1.0 / (gamma * gamma);

    if (beta2 <= 0.0) {
        return result;
    }

	const double beta = std::sqrt(beta2);
    const double v_cm_s = beta * C_CM_PER_S;  // beam velocity in cm/s

	// Shima reduced velocity:
    //   X = v / [ (3.6×10^8) * Z_1^0.45 ]
    const double denom = X_DENOM_SCALE * std::pow(Z1, 0.45);
    if (denom <= 0.0) {
        return result;
    }

	const double X = v_cm_s / denom;
    if (X <= 0.0) {
        return result;
    }

	// Shima mean charge state formula:
    //
    // q̄ / Z_1 = [ 1 - exp( -1.25 X + 0.32 X^2 - 0.11 X^3 ) ]
    //           × [ 1 - 0.0019 (Z_2 - 6) sqrt(X) + 1e-5 X (Z_2 - 6)^2 ].
    const double X2 = X * X;
    const double X3 = X2 * X;

	const double poly = -1.25 * X + 0.32 * X2 - 0.11 * X3;
    double bracket1 = 1.0 - std::exp(poly);
	bracket1 = std::max(0.0, std::min(1.0, bracket1));

	const double dZ  = Z2 - 6.0;
    double bracket2 =
        1.0
        - 0.0019 * dZ * std::sqrt(X)
        + 1.0e-5 * X * dZ * dZ;

	// Physically, the second factor should not drive q̄/Z1 negative.
    if (bracket2 < 0.0) bracket2 = 0.0;


	const double qFrac = bracket1 * bracket2;  // q̄ / Z_1
    double qbar = Z1 * qFrac;

	if (qbar < 0.0) qbar = 0.0;
    if (qbar > Z1)  qbar = Z1;

	result.meanQ  = qbar;

	// --- Schiwietz & Grande reduced-width model ---
    //
    // They define:
    //   d = [ Σ (q - q_m)^2 P_q ]^{1/2}   (RMS width in charge units)
    //   w = d * Z_p^{-0.27}
    //          * Z_t^{0.035 - 0.0009 Z_p}
    //          * f(q_m) * f(Z_p - q_m),
    //
    // with
    //   f(x) = sqrt((x + 0.37 Z_p^{0.6}) / x).
    //
    // Here we treat w as a user-supplied parameter and invert this relation
    // to obtain d = sigma_q.

	if (reducedWidthW <= 0.0 || qbar <= 0.0 || qbar >= Z1) {
        // At extreme charge fractions, the width should be small; we set σ=0.
        result.sigmaQ = 0.0;
        return result;
    }

	const double qp = Z1;
    const double qm = qbar;

	auto f_func = [qp](double x) -> double {
        if (x <= 0.0) return 1.0;  // avoid division by zero; width will be small anyway
        const double term = x + 0.37 * safePow(qp, 0.6);
        if (term <= 0.0) return 1.0;
        return std::sqrt(term / x);
    };

	const double f_qm      = f_func(qm);
    const double f_Zp_qm   = f_func(qp - qm);

	double S = safePow(qp, -0.27)
               * safePow(Z2, 0.035 - 0.0009 * qp)
               * f_qm * f_Zp_qm;

	if (S <= 0.0) {
        result.sigmaQ = 0.0;
        return result;
    }

	const double d = reducedWidthW / S;  // RMS width in charge units

    result.sigmaQ = (d > 0.0) ? d : 0.0;
    return result;


	// Width model: binomial-like variance scaled by an empirical factor.
    //
    //   var ≈ q̄ (1 - q̄ / Z_1)
    //   σ_q = scaleSigma * sqrt(var)
    //
    // scaleSigma can be tuned against experimental widths; here we keep
    // a moderate factor to avoid unrealistically narrow distributions.
    //double var = qbar * (1.0 - qbar / Z1);
    //if (var < 0.0) var = 0.0;
    //const double scaleSigma = 0.6;
    //double sigma = (var > 0.0) ? (scaleSigma * std::sqrt(var)) : 0.0;


    //result.sigmaQ = sigma;
    //return result;

}

ChargeStateResult ChargeStateCalculator::sample(
    int Zp, int Ap, double energyMeV, double targetZ, double reducedWidthW,
    std::mt19937_64& random) const
{
    ChargeStateResult result = computeMoments(Zp, Ap, energyMeV, targetZ, reducedWidthW);

    if (result.sigmaQ <= 0.0 || result.meanQ <= 0.0 || Zp <= 0) {
        result.sampledQ = 0;
        return result;
    }

    std::normal_distribution<double> gauss(result.meanQ, result.sigmaQ);
    double qSample = gauss(random);

    int qInt = static_cast<int>(std::floor(qSample + 0.5));
    qInt = std::max(0, std::min(Zp, qInt));

    result.sampledQ = qInt;
    return result;
}

} // namespace BeamMC