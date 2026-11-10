#include "TwoBodyKinematics.h"

#include <algorithm>
#include <stdexcept>

namespace BeamMC {
namespace {

const double PI = 3.14159265358979323846;
const double AMU_MEV = 931.49410242;

KinematicsVector unit(const KinematicsVector& v)
{
    const double m = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (m <= 0.0) return KinematicsVector(0.0, 0.0, 1.0);
    return KinematicsVector(v.x/m, v.y/m, v.z/m);
}

KinematicsVector add(const KinematicsVector& a, const KinematicsVector& b)
{
    return KinematicsVector(a.x+b.x, a.y+b.y, a.z+b.z);
}

KinematicsVector scale(const KinematicsVector& a, double s)
{
    return KinematicsVector(a.x*s, a.y*s, a.z*s);
}

KinematicsVector cross(const KinematicsVector& a, const KinematicsVector& b)
{
    return KinematicsVector(a.y*b.z-a.z*b.y,
                            a.z*b.x-a.x*b.z,
                            a.x*b.y-a.y*b.x);
}

KinematicsVector rotateFromBeam(const KinematicsVector& beam,
                                double theta, double phi)
{
    const KinematicsVector z = unit(beam);
    KinematicsVector reference = std::fabs(z.z) < 0.9
        ? KinematicsVector(0.0, 0.0, 1.0)
        : KinematicsVector(0.0, 1.0, 0.0);
    const KinematicsVector x = unit(cross(reference, z));
    const KinematicsVector y = unit(cross(z, x));
    const KinematicsVector transverse = add(scale(x, std::cos(phi)),
                                            scale(y, std::sin(phi)));
    return unit(add(scale(z, std::cos(theta)),
                    scale(transverse, std::sin(theta))));
}

double angleX(const KinematicsVector& d) { return std::atan2(d.x, d.z); }
double angleY(const KinematicsVector& d) { return std::atan2(d.y, d.z); }

} // namespace

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
      ejectileZ_(ejectileZ), ejectileA_(ejectileA), recoilZ_(recoilZ),
      recoilA_(recoilA), qValueMeV_(qValueMeV),
      recoilExcitationMeV_(recoilExcitationMeV), result_() {}

void TwoBodyKinematics::calculate(double beamEnergyMeV,
                                  double theta,
                                  double phi,
                                  const KinematicsVector& beamDirection)
{
    if (beamEnergyMeV < 0.0) throw std::invalid_argument("Negative beam energy");
    if (beamA_ <= 0 || targetA_ <= 0 || ejectileA_ <= 0 || recoilA_ <= 0)
        throw std::invalid_argument("Invalid reaction mass number");

    const double mb = beamA_ * AMU_MEV;
    const double mt = targetA_ * AMU_MEV;
    const double me = ejectileA_ * AMU_MEV;
    const double mr = recoilA_ * AMU_MEV;
    const double available = std::max(0.0, beamEnergyMeV + qValueMeV_ - recoilExcitationMeV_);
    const double cmEnergy = available * mt / (mb + mt);
    const double p = std::sqrt(std::max(0.0, 2.0 * (me*mr/(me+mr)) * cmEnergy));
    const double ejectileCmE = p*p/(2.0*me);
    const double recoilCmE = p*p/(2.0*mr);
    const double vCm = std::sqrt(std::max(0.0, 2.0*beamEnergyMeV*mb)) / (mb+mt);
    const double ejectileLabE = ejectileCmE + 0.5*me*vCm*vCm
        + std::sqrt(2.0*ejectileCmE*me)*vCm*std::cos(theta);
    const double recoilTheta = PI - theta;
    const double recoilLabE = recoilCmE + 0.5*mr*vCm*vCm
        + std::sqrt(2.0*recoilCmE*mr)*vCm*std::cos(recoilTheta);

    result_.ejectileEnergyMeV = std::max(0.0, ejectileLabE);
    result_.recoilEnergyMeV = std::max(0.0, recoilLabE);
    result_.ejectileDirection = rotateFromBeam(beamDirection, theta, phi);
    result_.recoilDirection = rotateFromBeam(beamDirection, recoilTheta, phi + PI);
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
