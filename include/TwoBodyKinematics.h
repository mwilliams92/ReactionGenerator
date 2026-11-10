#ifndef BEAMMC_TWO_BODY_KINEMATICS_H
#define BEAMMC_TWO_BODY_KINEMATICS_H

#include <cmath>

namespace BeamMC {

struct KinematicsVector {
    double x;
    double y;
    double z;
    KinematicsVector() : x(0.0), y(0.0), z(1.0) {}
    KinematicsVector(double xValue, double yValue, double zValue)
        : x(xValue), y(yValue), z(zValue) {}
};

struct KinematicsResult {
    double ejectileEnergyMeV;
    double recoilEnergyMeV;
    double ejectileHorizontalAngleRad;
    double ejectileVerticalAngleRad;
    double recoilHorizontalAngleRad;
    double recoilVerticalAngleRad;
    KinematicsVector ejectileDirection;
    KinematicsVector recoilDirection;
    KinematicsResult();
};

class TwoBodyKinematics {
public:
    TwoBodyKinematics(int beamZ, int beamA,
                      int targetZ, int targetA,
                      int ejectileZ, int ejectileA,
                      int recoilZ, int recoilA,
                      double qValueMeV,
                      double recoilExcitationMeV = 0.0);

    void calculate(double beamEnergyMeV,
                   double centreOfMassThetaRad,
                   double centreOfMassPhiRad = 0.0,
                   const KinematicsVector& beamDirection = KinematicsVector());

    double getEjectileEnergyMeV() const;
    double getRecoilEnergyMeV() const;
    double getEjectileHorizontalAngleRad() const;
    double getEjectileVerticalAngleRad() const;
    double getRecoilHorizontalAngleRad() const;
    double getRecoilVerticalAngleRad() const;
    const KinematicsVector& getEjectileDirection() const;
    const KinematicsVector& getRecoilDirection() const;
    const KinematicsResult& getResult() const;

private:
    int beamZ_;
    int beamA_;
    int targetZ_;
    int targetA_;
    int ejectileZ_;
    int ejectileA_;
    int recoilZ_;
    int recoilA_;
    double qValueMeV_;
    double recoilExcitationMeV_;
    KinematicsResult result_;
};

} // namespace BeamMC

#endif
