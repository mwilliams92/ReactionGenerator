#include "TwoBodyKinematics.h"

#include <cmath>
#include <stdexcept>

namespace BeamMC {
namespace {

const double PI = 3.14159265358979323846;
const double AMU_MEV = 931.49410242;

KinematicsVector unit(const KinematicsVector& v)
{
    const double magnitude = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (magnitude <= 0.0) return KinematicsVector(0.0, 0.0, 1.0);
    return KinematicsVector(v.x/magnitude, v.y/magnitude, v.z/magnitude);
}

KinematicsVector add(const KinematicsVector& a, const KinematicsVector& b)
{
    return KinematicsVector(a.x+b.x, a.y+b.y, a.z+b.z);
}

KinematicsVector scale(const KinematicsVector& v, double factor)
{
    return KinematicsVector(v.x*factor, v.y*factor, v.z*factor);
}

KinematicsVector cross(const KinematicsVector& a, const KinematicsVector& b)
{
    return KinematicsVector(a.y*b.z - a.z*b.y,
                            a.z*b.x - a.x*b.z,
                            a.x*b.y - a.y*b.x);
}

double magnitudeSquared(const KinematicsVector& v)
{
    return v.x*v.x + v.y*v.y + v.z*v.z;
}

// Interpret theta and phi as ejectile angles in the centre-of-mass frame,
// defined relative to the incident beam direction (not the laboratory z axis).
KinematicsVector rotateFromBeam(const KinematicsVector& beamDirection,
                                double theta, double phi)
{
    const KinematicsVector z = unit(beamDirection);
    const KinematicsVector reference = std::fabs(z.z) < 0.9
        ? KinematicsVector(0.0, 0.0, 1.0)
        : KinematicsVector(0.0, 1.0, 0.0);
    const KinematicsVector x = unit(cross(reference, z));
    const KinematicsVector y = cross(z, x);

    return unit(add(scale(z, std::cos(theta)),
                    add(scale(x, std::sin(theta)*std::cos(phi)),
                        scale(y, std::sin(theta)*std::sin(phi)))));
}

double angleX(const KinematicsVector& d) { return std::atan2(d.x, d.z); }
double angleY(const KinematicsVector& d) { return std::atan2(d.y, d.z); }

} // anonymous namespace

KinematicsResult::KinematicsResult()
    : ejectileEnergyMeV(0.0), recoilEnergyMeV(0.0),
      ejectileHorizontalAngleRad(0.0), ejectileVerticalAngleRad(0.0),
      recoilHorizontalAngleRad(0.0), recoilVerticalAngleRad(0.0),
      ejectileDirection(), recoilDirection() {}

TwoBodyKinematics::TwoBodyKinematics(int beamZ, int beamA,
                                     int targetZ, int targetA,
                                     int ejectileZ, int ejectileA,
                                     int recoilZ, int recoilA,
                                     double qValueMeV,
                                     double recoilExcitationMeV)
    : beamZ_(beamZ), beamA_(beamA), targetZ_(targetZ), targetA_(targetA),
      ejectileZ_(ejectileZ), ejectileA_(ejectileA),
      recoilZ_(recoilZ), recoilA_(recoilA),
      qValueMeV_(qValueMeV), recoilExcitationMeV_(recoilExcitationMeV),
      result_() {}

void TwoBodyKinematics::calculate(double beamEnergyMeV,
                                  double centreOfMassThetaRad,
                                  double centreOfMassPhiRad,
                                  const KinematicsVector& beamDirection)
{
    // Clear the previous event, including when the reaction is forbidden.
    result_ = KinematicsResult();

    if (!std::isfinite(beamEnergyMeV) || beamEnergyMeV < 0.0)
        throw std::invalid_argument("Invalid beam kinetic energy");
    if (!std::isfinite(centreOfMassThetaRad) ||
        !std::isfinite(centreOfMassPhiRad))
        throw std::invalid_argument("Invalid centre-of-mass angle");
    if (!std::isfinite(beamDirection.x) ||
        !std::isfinite(beamDirection.y) ||
        !std::isfinite(beamDirection.z) ||
        magnitudeSquared(beamDirection) <= 0.0)
        throw std::invalid_argument("Invalid beam direction");
    if (beamA_ <= 0 || targetA_ <= 0 || ejectileA_ <= 0 || recoilA_ <= 0)
        throw std::invalid_argument("Invalid reaction mass number");
    if (!std::isfinite(qValueMeV_) || !std::isfinite(recoilExcitationMeV_))
        throw std::invalid_argument("Invalid Q-value or excitation energy");

    // Non-relativistic masses in MeV/c^2; energy and momentum use c = 1.
    // This is intentionally the same mass approximation as the original BeamMC.
    const double mb = beamA_ * AMU_MEV;
    const double mt = targetA_ * AMU_MEV;
    const double me = ejectileA_ * AMU_MEV;
    const double mr = recoilA_ * AMU_MEV;

    // Q and excitation energy must be added AFTER transforming beam energy to CM.
    const double beamCmEnergy = beamEnergyMeV * mt / (mb + mt);
    const double availableCmEnergy =
        beamCmEnergy + qValueMeV_ - recoilExcitationMeV_;

    // Preserve original BeamMC behaviour for a closed reaction channel:
    // all result energies remain zero, so its transport is skipped.
    if (availableCmEnergy <= 0.0) return;

    const double reducedFinalMass = me * mr / (me + mr);
    const double finalCmMomentum =
        std::sqrt(2.0 * reducedFinalMass * availableCmEnergy);

    const double initialLabMomentum = std::sqrt(2.0 * mb * beamEnergyMeV);
    const double cmVelocity = initialLabMomentum / (mb + mt);

    const KinematicsVector zAxis = unit(beamDirection);
    const KinematicsVector ejectileCmDirection =
        rotateFromBeam(zAxis, centreOfMassThetaRad, centreOfMassPhiRad);

    // Galilean transformation from CM to lab: p_lab = p_CM + m*v_CM.
    // The two final CM momenta are equal and opposite.
    const KinematicsVector ejectileLabMomentum =
        add(scale(ejectileCmDirection, finalCmMomentum),
            scale(zAxis, me * cmVelocity));
    const KinematicsVector recoilLabMomentum =
        add(scale(ejectileCmDirection, -finalCmMomentum),
            scale(zAxis, mr * cmVelocity));

    result_.ejectileEnergyMeV = magnitudeSquared(ejectileLabMomentum) / (2.0 * me);
    result_.recoilEnergyMeV = magnitudeSquared(recoilLabMomentum) / (2.0 * mr);
    result_.ejectileDirection = unit(ejectileLabMomentum);
    result_.recoilDirection = unit(recoilLabMomentum);
    result_.ejectileHorizontalAngleRad = angleX(result_.ejectileDirection);
    result_.ejectileVerticalAngleRad = angleY(result_.ejectileDirection);
    result_.recoilHorizontalAngleRad = angleX(result_.recoilDirection);
    result_.recoilVerticalAngleRad = angleY(result_.recoilDirection);
}

double TwoBodyKinematics::getEjectileEnergyMeV() const { return result_.ejectileEnergyMeV; }
double TwoBodyKinematics::getRecoilEnergyMeV() const { return result_.recoilEnergyMeV; }
double TwoBodyKinematics::getEjectileHorizontalAngleRad() const { return result_.ejectileHorizontalAngleRad; }
double TwoBodyKinematics::getEjectileVerticalAngleRad() const { return result_.ejectileVerticalAngleRad; }
double TwoBodyKinematics::getRecoilHorizontalAngleRad() const { return result_.recoilHorizontalAngleRad; }
double TwoBodyKinematics::getRecoilVerticalAngleRad() const { return result_.recoilVerticalAngleRad; }
const KinematicsVector& TwoBodyKinematics::getEjectileDirection() const { return result_.ejectileDirection; }
const KinematicsVector& TwoBodyKinematics::getRecoilDirection() const { return result_.recoilDirection; }
const KinematicsResult& TwoBodyKinematics::getResult() const { return result_; }

} // namespace BeamMC
