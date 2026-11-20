// BeamMC.cpp
// -----------------------------------------------------------------------------
// Monte Carlo simulation of particle transport through a target and optional
// degrader. The program supports:
//   * SRIM stopping-power tables for compound materials
//   * Beam energy spread, transverse emittance, energy straggling, and scattering
//   * Beam transport mode
//   * Two-body reaction mode with a randomly sampled reaction depth
//   * Recoil transport through the remaining target and optional degrader
//
// Build:
//   rm -rf build
//   cmake -S . -B build -DBEAMMC_WITH_ROOT=ON
//   cmake --build build -j
//
// Run:
//   ./BeamMC beam
//   ./BeamMC reaction
// -----------------------------------------------------------------------------

#include "TwoBodyKinematics.h"
#include "ExcitationPopulation.h"
#include "ReactionCrossSection.h"
#include "ChargeState.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#ifdef BEAMMC_WITH_ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#endif
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double AMU_MEV = 931.49410242;

// -----------------------------------------------------------------------------
// Input data structures
// -----------------------------------------------------------------------------

struct Material {
    double thicknessMm = 0.0;
    double density = 0.0;              // g/cm^3
    double effectiveZ = 6.0;           // Used for approximate scattering only
    double effectiveA = 12.0;          // Used for approximate scattering only
    std::string srimFile;
    std::string energyUnit = "auto";   // "auto" reads SRIM header
    std::string stoppingUnit = "auto";
};

struct BeamConfig {
    int Z = 0;
    int A = 0;
    double energyValue = 0.0;
    std::string energyUnit = "MeV";

    double energySpreadValue = 0.0;
    std::string energySpreadMode = "sigma_abs";

    double energyMeV = 0.0;
    double energySpreadMeV = 0.0;

    double emittanceX = 0.0;           // mm mrad
    double emittanceY = 0.0;           // mm mrad
    double sigmaXmm = 1.0;
    double sigmaYmm = 1.0;

    double stepMm = 0.01;
    double stragglingCoefficient = 0.02;

    int particles = 10000;
    unsigned long long seed = 0;

    std::string output = "beam_output.csv";

    bool rootOutput = false;
    std::string rootFile = "beam_output.root";

	double reducedWidthW = 0.7;
};

struct ReactionConfig {
    int beamZ = 0;
    int beamA = 0;
    int targetZ = 0;
    int targetA = 0;
    int ejectileZ = 0;
    int ejectileA = 0;
    int recoilZ = 0;
    int recoilA = 0;

    double excitationMeV = 0.0;
    double qValueMeV = 0.0;
    double neutronSeparationMeV = -1.0;  // optional cutoff for Ex

    std::string excitationPopulationFile;
    std::string reactionCrossSectionFile;

    std::string beamTargetSRIM;
    std::string recoilTargetSRIM;
    std::string recoilDegraderSRIM;

    std::string output = "reaction_output.csv";
};

struct Point {
    double x;
    double y;
    double z;

    Point() :
        x(0.0), y(0.0), z(0.0) {}
    Point(double xValue, double yValue, double zValue)
        : x(xValue), y(yValue), z(zValue) {}
};

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() : x(0.0), y(0.0), z(1.0) {}
    Vec3(double xValue, double yValue, double zValue)
        : x(xValue), y(yValue), z(zValue) {}
};

struct Particle {
    double E;              // MeV
    Point position;        // mm
    Vec3 direction;        // unit vector
    double thetaX;         // mrad
    double thetaY;         // mrad
    int    chargeState;    // exit charge state (set at end)

    Particle()
        : E(0.0), position(), direction(),
          thetaX(0.0), thetaY(0.0), chargeState(0) {}
};

struct DegraderConfig : public Material {
    bool enabled = false;
};

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

std::string lower(const std::string& text)
{
    std::string result = text;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(result[i])));
    }
    return result;
}

std::string trim(const std::string& text)
{
    std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return std::string();
    std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

// -----------------------------------------------------------------------------
// Input parsing
// -----------------------------------------------------------------------------

std::map<std::string, std::string> readKeyValueFile(const std::string& filename)
{
    std::ifstream input(filename.c_str());
    if (!input) {
        throw std::runtime_error("Cannot open input file " + filename);
    }

    std::map<std::string, std::string> values;
    std::string line;

    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty()) continue;
        if (stripped[0] == '#') continue;

        std::istringstream stream(stripped);
        std::string key;
        std::string value;

        stream >> key;
        if (key.empty()) continue;

        if (stripped.find('=') != std::string::npos) {
            std::getline(stream, value);
            const std::size_t eqPos = value.find('=');
            if (eqPos != std::string::npos) {
                value = value.substr(eqPos + 1);
            }
        } else {
            std::getline(stream, value);
        }

        value = trim(value);
        values[lower(key)] = value;
    }

    return values;
}

int getInt(const std::map<std::string, std::string>& values,
           const std::string& key,
           int defaultValue = 0)
{
    const auto it = values.find(lower(key));
    if (it == values.end()) return defaultValue;
    return std::atoi(it->second.c_str());
}

double getDouble(const std::map<std::string, std::string>& values,
                 const std::string& key,
                 double defaultValue = 0.0)
{
    const auto it = values.find(lower(key));
    if (it == values.end()) return defaultValue;
    return std::atof(it->second.c_str());
}

std::string getString(const std::map<std::string, std::string>& values,
                      const std::string& key,
                      const std::string& defaultValue = std::string())
{
    const auto it = values.find(lower(key));
    if (it == values.end()) return defaultValue;
    return it->second;
}

BeamConfig readBeam(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    BeamConfig beam;

    beam.Z = getInt(values, "beam_Z");
    beam.A = getInt(values, "beam_A");
    beam.energyValue = getDouble(values, "beam_energy");
    beam.energyUnit = getString(values, "beam_energy_unit", beam.energyUnit);

    beam.energySpreadValue = getDouble(values, "energy_spread");
    beam.energySpreadMode = getString(values, "energy_spread_mode",
                                      beam.energySpreadMode);

    beam.emittanceX = getDouble(values, "emittance_x_mm_mrad");
    beam.emittanceY = getDouble(values, "emittance_y_mm_mrad");
    beam.sigmaXmm = getDouble(values, "beam_sigma_x_mm", beam.sigmaXmm);
    beam.sigmaYmm = getDouble(values, "beam_sigma_y_mm", beam.sigmaYmm);

    beam.particles = getInt(values, "particles", beam.particles);
    beam.stepMm = getDouble(values, "step_um", 0.1) * 0.001;
    beam.stragglingCoefficient = getDouble(
        values, "bohr_straggling_scale", beam.stragglingCoefficient);

    beam.seed = static_cast<unsigned long long>(
        getInt(values, "seed", static_cast<int>(beam.seed)));
    beam.output = getString(values, "output", beam.output);

    beam.rootOutput = (getInt(values, "root_output", 0) != 0);
    beam.rootFile = getString(values, "root_file", beam.rootFile);

	beam.reducedWidthW = getDouble(
        values, "charge_state_reduced_width", 0.7);

    if (beam.A <= 0 || beam.energyValue <= 0.0) {
        throw std::runtime_error("Beam input needs BEAM_A and BEAM_ENERGY > 0");
    }

    if (beam.energyUnit == "MeV" || beam.energyUnit == "mev") {
        beam.energyMeV = beam.energyValue;
        beam.energySpreadMeV = beam.energySpreadValue;
    } else {
        beam.energyMeV = beam.energyValue * beam.A;
        beam.energySpreadMeV = beam.energySpreadValue * beam.A;
    }

    if (beam.energySpreadMode == "sigma_percent" ||
        beam.energySpreadMode == "fwhm_percent") {
        const double fraction = beam.energySpreadValue / 100.0;
        if (beam.energySpreadMode == "sigma_percent") {
            beam.energySpreadMeV = fraction * beam.energyMeV;
        } else {
            beam.energySpreadMeV = fraction * beam.energyMeV /
                (2.0 * std::sqrt(2.0 * std::log(2.0)));
        }
    } else if (beam.energySpreadMode == "fwhm_abs") {
        beam.energySpreadMeV = beam.energySpreadMeV /
            (2.0 * std::sqrt(2.0 * std::log(2.0)));
    }

    return beam;
}

Material readMaterial(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    Material material;

    material.thicknessMm = getDouble(values, "thickness_mm");
    material.density = getDouble(values, "density_g_cm3");
    material.effectiveZ = getDouble(values, "effective_Z", material.effectiveZ);
    material.effectiveA = getDouble(values, "effective_A", material.effectiveA);
    material.srimFile = getString(values, "stopping_table");
    material.energyUnit = getString(values, "srim_energy_unit",
                                    material.energyUnit);
    material.stoppingUnit = getString(values, "srim_stopping_unit",
                                      material.stoppingUnit);

    if (material.thicknessMm <= 0.0 || material.density <= 0.0 ||
        material.srimFile.empty()) {
        throw std::runtime_error(
            "Target input needs THICKNESS_MM, DENSITY_G_CM3, and STOPPING_TABLE");
    }

    return material;
}

DegraderConfig readDegrader(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    DegraderConfig degrader;

    degrader.enabled = (getInt(values, "enabled", 0) != 0);
    degrader.thicknessMm = getDouble(values, "thickness_mm", degrader.thicknessMm);
    degrader.density = getDouble(values, "density_g_cm3", degrader.density);
    degrader.effectiveZ = getDouble(values, "effective_Z", degrader.effectiveZ);
    degrader.effectiveA = getDouble(values, "effective_A", degrader.effectiveA);
    degrader.srimFile = getString(values, "stopping_table", degrader.srimFile);
    degrader.energyUnit = getString(values, "srim_energy_unit", degrader.energyUnit);
    degrader.stoppingUnit = getString(values, "srim_stopping_unit",
                                      degrader.stoppingUnit);

    if (degrader.enabled) {
        if (degrader.thicknessMm <= 0.0 || degrader.density <= 0.0 ||
            degrader.srimFile.empty()) {
            throw std::runtime_error(
                "Degrader input needs THICKNESS_MM, DENSITY_G_CM3, and STOPPING_TABLE when ENABLED != 0");
        }
    }

    return degrader;
}

ReactionConfig readReaction(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    ReactionConfig reaction;

    reaction.beamZ = getInt(values, "beam_Z");
    reaction.beamA = getInt(values, "beam_A");
    reaction.targetZ = getInt(values, "target_Z");
    reaction.targetA = getInt(values, "target_A");
    reaction.ejectileZ = getInt(values, "ejectile_Z");
    reaction.ejectileA = getInt(values, "ejectile_A");
    reaction.recoilZ = getInt(values, "recoil_Z");
    reaction.recoilA = getInt(values, "recoil_A");

    reaction.excitationMeV = getDouble(values, "recoil_excitation_MeV");
    reaction.qValueMeV = getDouble(values, "q_value_MeV");
    reaction.neutronSeparationMeV = getDouble(
        values, "recoil_neutron_separation_mev", -1.0);

    reaction.beamTargetSRIM = getString(values, "beam_srim_target");
    reaction.recoilTargetSRIM = getString(values, "recoil_srim_target");
    reaction.recoilDegraderSRIM = getString(values, "recoil_srim_degrader");

    reaction.excitationPopulationFile = getString(
        values, "recoil_excitation_population_file");
    reaction.reactionCrossSectionFile = getString(
        values, "reaction_cross_section_file");

    reaction.output = getString(values, "output", reaction.output);

    if (reaction.beamA <= 0 || reaction.targetA <= 0 ||
        reaction.ejectileA <= 0 || reaction.recoilA <= 0 ||
        reaction.beamTargetSRIM.empty() || reaction.recoilTargetSRIM.empty()) {
        throw std::runtime_error(
            "Reaction input needs particle data, beam_srim_target, and "
            "recoil_srim_target");
    }

    return reaction;
}

// -----------------------------------------------------------------------------
// SRIM table handling
// -----------------------------------------------------------------------------

struct SRIMEntry {
    double energyMeV;
    double stoppingMeVPerMm;
};

struct SRIMTable {
    std::vector<SRIMEntry> entries;
};

SRIMTable readSRIMTable(const Material& material)
{
    std::ifstream input(material.srimFile.c_str());
    if (!input) {
        throw std::runtime_error("Cannot open SRIM stopping table " + material.srimFile);
    }

    std::string line;
    bool inTable = false;

    std::vector<double> energies;
    std::vector<double> totalStopping;

    while (std::getline(input, line)) {
        std::string stripped = trim(line);
        if (stripped.empty()) continue;

        if (!inTable) {
            // Accept both "Ion Energy ..." and "Ion dE/dx ..." header styles.
            if (stripped.find("Ion")    != std::string::npos &&
                (stripped.find("Energy") != std::string::npos ||
                 stripped.find("dE/dx")  != std::string::npos)) {
                inTable = true;
            }
            continue;
        }

        // Interpret data rows of the form:
        //   10.00 keV   4.323E+01  7.385E+02 ...
        std::istringstream stream(stripped);

        double energy = 0.0;
        if (!(stream >> energy)) {
            continue;
        }

        std::string token;
        if (!(stream >> token)) {
            continue;
        }

        double energyKeV = 0.0;
        double elec = 0.0;
        double nucl = 0.0;

        try {
            // Numeric token: actually electronic stopping; energy already keV.
            elec = std::stod(token);
            energyKeV = energy;
            if (!(stream >> nucl)) {
                continue;
            }
        } catch (const std::exception&) {
            const std::string unit = lower(token);
            double scaleToKeV = 1.0;

            if (unit.find("kev") != std::string::npos) {
                scaleToKeV = 1.0;
            } else if (unit.find("mev") != std::string::npos) {
                scaleToKeV = 1.0e3;
            } else if (unit.find("gev") != std::string::npos) {
                scaleToKeV = 1.0e6;
            } else {
                scaleToKeV = 1.0;
            }

            energyKeV = energy * scaleToKeV;

            if (!(stream >> elec >> nucl)) {
                continue;
            }
        }

        energies.push_back(energyKeV);
        totalStopping.push_back(elec + nucl);
    }

    if (energies.empty()) {
        throw std::runtime_error(
            "No data rows found in SRIM table " + material.srimFile);
    }

    double energyScale = 1.0;
    double stoppingScale = 1.0;

    const std::string energyUnit = lower(material.energyUnit);
    const std::string stoppingUnit = lower(material.stoppingUnit);

    if (energyUnit == "auto" || energyUnit.empty() ||
        energyUnit == "kev") {
        energyScale = 0.001;  // keV -> MeV
    } else if (energyUnit == "mev") {
        energyScale = 1.0;
    } else {
        throw std::runtime_error(
            "Unsupported SRIM energy unit: " + material.energyUnit);
    }

    if (stoppingUnit == "auto" || stoppingUnit.empty() ||
        stoppingUnit == "kev/micron") {
        stoppingScale = 0.001; // keV/micron -> MeV/micron
    } else if (stoppingUnit == "mev/micron") {
        stoppingScale = 1.0;
    } else if (stoppingUnit == "mev/mm") {
        stoppingScale = 1.0 / 1000.0;
    } else {
        throw std::runtime_error(
            "Unsupported SRIM stopping unit: " + material.stoppingUnit);
    }

    SRIMTable table;
    const double micronToMm = 0.001;

    for (std::size_t i = 0; i < energies.size(); ++i) {
        SRIMEntry entry;
        entry.energyMeV = energies[i] * energyScale;
        const double stoppingMeVPerMicron = totalStopping[i] * stoppingScale;
        entry.stoppingMeVPerMm = stoppingMeVPerMicron / micronToMm;
        table.entries.push_back(entry);
    }

    return table;
}

struct TransportModel {
    SRIMTable table;
    Material material;
};

TransportModel buildTransportModel(const Material& material)
{
    TransportModel model;
    model.table = readSRIMTable(material);
    model.material = material;
    return model;
}

double interpolateStoppingPower(const SRIMTable& table, double energyMeV)
{
    if (table.entries.empty()) return 0.0;

    if (energyMeV <= table.entries.front().energyMeV) {
        return table.entries.front().stoppingMeVPerMm;
    }

    if (energyMeV >= table.entries.back().energyMeV) {
        return table.entries.back().stoppingMeVPerMm;
    }

    std::size_t upperIndex = 1;
    while (upperIndex < table.entries.size() &&
           table.entries[upperIndex].energyMeV < energyMeV) {
        ++upperIndex;
    }

    if (upperIndex >= table.entries.size()) {
        return table.entries.back().stoppingMeVPerMm;
    }

    const std::size_t lowerIndex = upperIndex - 1;
    const double e1 = table.entries[lowerIndex].energyMeV;
    const double e2 = table.entries[upperIndex].energyMeV;
    const double s1 = table.entries[lowerIndex].stoppingMeVPerMm;
    const double s2 = table.entries[upperIndex].stoppingMeVPerMm;

    if (e2 <= e1) return s1;

    const double t = (energyMeV - e1) / (e2 - e1);
    return s1 + t * (s2 - s1);
}

double energyLossPerMillimetre(const TransportModel& model, double energyMeV)
{
    return interpolateStoppingPower(model.table, std::max(energyMeV, 1.0e-9));
}

// -----------------------------------------------------------------------------
// Reaction depth distribution using TALYS σ(E)
// -----------------------------------------------------------------------------

struct ReactionDepthDistribution {
    std::vector<double> depthsMm;
    std::vector<double> cdf;
    double totalWeight = 0.0;

    bool empty() const {
        return depthsMm.empty() || totalWeight <= 0.0;
    }

    double sample(std::mt19937_64& random,
                  const Material& target,
                  double stepMm) const
    {
        if (empty()) {
            double u = std::generate_canonical<double, 64>(random);
            double depth = target.thicknessMm * u;
            depth = std::max(
                std::min(depth, target.thicknessMm),
                std::min(stepMm, target.thicknessMm));
            return depth;
        }

        std::uniform_real_distribution<double> uni(0.0, totalWeight);
        double u = uni(random);

        auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
        std::size_t idx = static_cast<std::size_t>(it - cdf.begin());
        if (idx >= depthsMm.size()) {
            idx = depthsMm.size() - 1;
        }

        double depth = depthsMm[idx];
        depth = std::max(
            std::min(depth, target.thicknessMm),
            std::min(stepMm, target.thicknessMm));
        return depth;
    }
};

ReactionDepthDistribution buildReactionDepthDistribution(
    const BeamConfig& beam,
    const Material& target,
    const TransportModel& targetModel,
    const BeamMC::ReactionCrossSection& xs,
    const ReactionConfig& reaction)
{
    ReactionDepthDistribution dist;

    if (xs.empty() || target.thicknessMm <= 0.0) {
        return dist;
    }

    const double dz = beam.stepMm;
    double depthMm = 0.0;
    double energyMeV = beam.energyMeV;

    double totalWeight = 0.0;

    while (depthMm < target.thicknessMm && energyMeV > 0.0) {
        const double step = std::min(dz, target.thicknessMm - depthMm);

        // Map inverse-kinematics beam energy to equivalent projectile lab energy
        // for the direct TALYS calculation: E_proj,lab = E_beam,lab * (A_target/A_beam).
        const double projectileEnergyMeV =
            energyMeV * static_cast<double>(reaction.targetA) /
            static_cast<double>(reaction.beamA);
        const double sigmaMb = std::max(0.0, xs.evaluate(projectileEnergyMeV));
        const double weight = sigmaMb * step;

        if (weight > 0.0) {
            totalWeight += weight;
            dist.depthsMm.push_back(depthMm + 0.5 * step);
            dist.cdf.push_back(totalWeight);
        }

        const double dE =
            energyLossPerMillimetre(targetModel, energyMeV) * step;
        energyMeV = std::max(0.0, energyMeV - dE);
        depthMm += step;
    }

    dist.totalWeight = totalWeight;
    return dist;
}

// -----------------------------------------------------------------------------
// Vector and geometry helpers
// -----------------------------------------------------------------------------

Vec3 unitVector(const Vec3& vector)
{
    const double magnitude = std::sqrt(
        vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

    if (magnitude <= 0.0) {
        return Vec3(0.0, 0.0, 1.0);
    }

    return Vec3(
        vector.x / magnitude,
        vector.y / magnitude,
        vector.z / magnitude);
}

Point propagatePosition(const Point& position,
                        const Vec3& direction,
                        double distanceMm)
{
    Point newPosition;
    newPosition.x = position.x + direction.x * distanceMm;
    newPosition.y = position.y + direction.y * distanceMm;
    newPosition.z = position.z + direction.z * distanceMm;
    return newPosition;
}

double horizontalAngle(const Vec3& direction)
{
    return std::atan2(direction.x, direction.z) * 1000.0;
}

double verticalAngle(const Vec3& direction)
{
    return std::atan2(direction.y, direction.z) * 1000.0;
}

// -----------------------------------------------------------------------------
// Beam generation and transport
// -----------------------------------------------------------------------------

Particle generateInitialBeamParticle(const BeamConfig& beam,
                                     std::mt19937_64& random)
{
    Particle particle;

    if (beam.energySpreadMeV > 0.0) {
        std::normal_distribution<double> gauss(
            beam.energyMeV, beam.energySpreadMeV);
        particle.E = std::max(0.0, gauss(random));
    } else {
        particle.E = beam.energyMeV;
    }

    particle.position = Point(0.0, 0.0, 0.0);
    particle.direction = Vec3(0.0, 0.0, 1.0);

    particle.thetaX = 0.0;
    particle.thetaY = 0.0;

    return particle;
}

typedef void (*EnergyDepthRecorder)(double depthMm, double energyMeV);

void propagate(Particle& particle,
               const TransportModel& model,
               double distanceMm,
               double stepMm,
               int Z,
               int A,
               double stragglingCoefficient,
               std::mt19937_64& random,
               EnergyDepthRecorder recorder = nullptr)
{
    const Material& material = model.material;
    const double density = material.density;

    double remainingDistance = distanceMm;

    std::normal_distribution<double> gauss(0.0, 1.0);

    while (remainingDistance > 0.0 && particle.E > 0.0) {
        const double step = std::min(remainingDistance, stepMm);

        const double dEdxMeVPerMm = energyLossPerMillimetre(model, particle.E);
        const double meanLoss = dEdxMeVPerMm * step;

        double variance = 0.0;
        if (stragglingCoefficient > 0.0) {
            const double arealDensity = density * step * 0.1; // g/cm^2
            const double Zt = material.effectiveZ;
            const double At = material.effectiveA;

            const double K = 0.307075;
            const double z = static_cast<double>(Z);
            const double variancePerPath =
                4.0 * PI * K * z * z * arealDensity * Zt / At;

            variance = stragglingCoefficient * variancePerPath;
        }

        double loss = meanLoss;
        if (variance > 0.0) {
            const double sigma = std::sqrt(variance);
            loss += sigma * gauss(random);
        }

        loss = std::max(loss, 0.0);
        loss = std::min(loss, particle.E);

        particle.E -= loss;
        particle.position = propagatePosition(
            particle.position, particle.direction, step);

        if (recorder) {
            recorder(particle.position.z, particle.E);
        }

        const double radiationLengthMm = 1000.0;
        const double t = step / radiationLengthMm;
        const double theta0 =
            13.6 * std::sqrt(t) * (1.0 + 0.038 * std::log(t));

        const double thetaXrad = theta0 * gauss(random) / 1000.0;
        const double thetaYrad = theta0 * gauss(random) / 1000.0;

        const double cosX = std::cos(thetaXrad);
        const double sinX = std::sin(thetaXrad);
        const double cosY = std::cos(thetaYrad);
        const double sinY = std::sin(thetaYrad);

        Vec3 direction = particle.direction;

        double dx = direction.x * cosX + direction.z * sinX;
        double dz = -direction.x * sinX + direction.z * cosX;
        direction.x = dx;
        direction.z = dz;

        double dy = direction.y * cosY + direction.z * sinY;
        dz = -direction.y * sinY + direction.z * cosY;
        direction.y = dy;
        direction.z = dz;

        particle.direction = unitVector(direction);
        particle.thetaX = horizontalAngle(particle.direction);
        particle.thetaY = verticalAngle(particle.direction);

        remainingDistance -= step;
    }
}

// -----------------------------------------------------------------------------
// Reaction generation
// -----------------------------------------------------------------------------

Particle generateRecoil(const Particle& beamParticle,
                        const ReactionConfig& reaction,
                        std::mt19937_64& random,
                        const BeamMC::ExcitationPopulation* population)
{
    const int beamZ = reaction.beamZ;
    const int beamA = reaction.beamA;
    const int targetZ = reaction.targetZ;
    const int targetA = reaction.targetA;
    const int ejectileZ = reaction.ejectileZ;
    const int ejectileA = reaction.ejectileA;
    const int recoilZ = reaction.recoilZ;
    const int recoilA = reaction.recoilA;

    double excitationMeV = reaction.excitationMeV;
    if (population && !population->empty()) {
        excitationMeV = population->sample(random);
    }

    const double qValueMeV = reaction.qValueMeV;

    BeamMC::TwoBodyKinematics kinematics(
        beamZ, beamA,
        targetZ, targetA,
        ejectileZ, ejectileA,
        recoilZ, recoilA,
        qValueMeV,
        excitationMeV);

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double cosThetaCM = 2.0 * uniform01(random) - 1.0;
    const double thetaCM = std::acos(cosThetaCM);
    const double phiCM = 2.0 * PI * uniform01(random);

    const Vec3 beamDirection = beamParticle.direction;
    BeamMC::KinematicsVector beamUnit(
        beamDirection.x, beamDirection.y, beamDirection.z);

    kinematics.calculate(
        beamParticle.E,
        thetaCM,
        phiCM,
        beamUnit);

    const BeamMC::KinematicsResult& result = kinematics.getResult();

    Particle recoil;
    recoil.E = result.recoilEnergyMeV;
    recoil.position = beamParticle.position;
    recoil.direction = Vec3(
        result.recoilDirection.x,
        result.recoilDirection.y,
        result.recoilDirection.z);
    recoil.thetaX = result.recoilHorizontalAngleRad * 1000.0;
    recoil.thetaY = result.recoilVerticalAngleRad * 1000.0;

    return recoil;
}

// -----------------------------------------------------------------------------
// ROOT output helper
// -----------------------------------------------------------------------------

#ifdef BEAMMC_WITH_ROOT

struct RootWriter {
    TFile* file;
    TH1D* hEnergy;
    TH2D* hAngles;
    TH1D* hExitCharge;

    TH1D* hReactionDepth;        // 1D: reaction depth in target
    TH2D* hBeamEnergyVsDepth;    // 2D: beam energy vs depth
    TH2D* hRecoilEnergyVsDepth;  // 2D: recoil energy vs depth

    TTree* tree;

    double incidentEnergy;
    double incidentThetaX;
    double incidentThetaY;
    double recoilEnergy;
    double recoilThetaX;
    double recoilThetaY;
    double exitCharge;

    RootWriter(const std::string& filename,
               double targetThicknessMm,
               bool degraderEnabled,
               double degraderThicknessMm,
               double maxBeamEnergyMeV,
               double beamStepMm)
    {
        const double totalThicknessMm =
            targetThicknessMm +
            (degraderEnabled ? degraderThicknessMm : 0.0);

        int nDepthBins = std::max(
            1,
            static_cast<int>(std::ceil(totalThicknessMm / beamStepMm)));

        file = new TFile(filename.c_str(), "RECREATE");

        hEnergy = new TH1D(
            "hEnergy",
            "Recoil energy;E_{recoil} [MeV];Counts",
            400, 0.0, maxBeamEnergyMeV);

        hAngles = new TH2D(
            "hAngles",
            "Recoil angles;theta_{x} [mrad];theta_{y} [mrad]",
            400, -200.0, 200.0,
            400, -200.0, 200.0);

        hReactionDepth = new TH1D(
            "hReactionDepth",
            "Reaction depth in target;depth [mm];Counts",
            nDepthBins, 0.0, targetThicknessMm);

        hBeamEnergyVsDepth = new TH2D(
            "hBeamEnergyVsDepth",
            "Beam energy vs depth;depth [mm];E_{beam} [MeV]",
            nDepthBins, 0.0, totalThicknessMm,
            400, 0.0, maxBeamEnergyMeV * 1.1);

        hRecoilEnergyVsDepth = new TH2D(
            "hRecoilEnergyVsDepth",
            "Recoil energy vs depth;depth [mm];E_{recoil} [MeV]",
            nDepthBins, 0.0, totalThicknessMm,
            400, 0.0, maxBeamEnergyMeV * 1.1);

        hExitCharge = new TH1D(
            "hExitCharge",
            "Exit charge state;q;Counts",
            100, -0.5, 99.5);

        tree = new TTree("BeamMC", "BeamMC output");

        tree->Branch("incidentEnergy",  &incidentEnergy, "incidentEnergy/D");
        tree->Branch("incidentThetaX",  &incidentThetaX, "incidentThetaX/D");
        tree->Branch("incidentThetaY",  &incidentThetaY, "incidentThetaY/D");
        tree->Branch("recoilEnergy",    &recoilEnergy,   "recoilEnergy/D");
        tree->Branch("recoilThetaX",    &recoilThetaX,   "recoilThetaX/D");
        tree->Branch("recoilThetaY",    &recoilThetaY,   "recoilThetaY/D");
        tree->Branch("exitCharge",      &exitCharge,     "exitCharge/D");
    }

    ~RootWriter()
    {
        if (file) {
            file->cd();
            if (hEnergy)             hEnergy->Write();
            if (hAngles)            hAngles->Write();
            if (hReactionDepth)     hReactionDepth->Write();
            if (hBeamEnergyVsDepth) hBeamEnergyVsDepth->Write();
            if (hRecoilEnergyVsDepth) hRecoilEnergyVsDepth->Write();
            if (hExitCharge)        hExitCharge->Write();
            if (tree)               tree->Write();
            file->Close();
        }
    }

    void fill(const Particle& incident,
              const Particle& recoilOrBeam,
              int qExit)
    {
        incidentEnergy = incident.E;
        incidentThetaX = incident.thetaX;
        incidentThetaY = incident.thetaY;
        recoilEnergy   = recoilOrBeam.E;
        recoilThetaX   = recoilOrBeam.thetaX;
        recoilThetaY   = recoilOrBeam.thetaY;
        exitCharge     = static_cast<double>(qExit);

        if (hEnergy)      hEnergy->Fill(recoilEnergy);
        if (hAngles)      hAngles->Fill(recoilThetaX, recoilThetaY);
        if (hExitCharge)  hExitCharge->Fill(qExit);
        if (tree)         tree->Fill();
    }

    void fillReactionDepth(double depthMm, double beamEnergyMeV)
    {
        if (hReactionDepth)     hReactionDepth->Fill(depthMm);
        if (hBeamEnergyVsDepth) hBeamEnergyVsDepth->Fill(depthMm, beamEnergyMeV);
    }

    void fillBeamDepth(double depthMm, double beamEnergyMeV)
    {
        if (hBeamEnergyVsDepth) hBeamEnergyVsDepth->Fill(depthMm, beamEnergyMeV);
    }

    void fillRecoilDepth(double depthMm, double recoilEnergyMeV)
    {
        if (hRecoilEnergyVsDepth) hRecoilEnergyVsDepth->Fill(depthMm, recoilEnergyMeV);
    }
};

// Global pointer for recorder callbacks
RootWriter* gRootWriter = nullptr;

void recordBeamEnergyDepth(double depthMm, double energyMeV)
{
    if (gRootWriter) gRootWriter->fillBeamDepth(depthMm, energyMeV);
}

#endif // BEAMMC_WITH_ROOT

} // anonymous namespace

// -----------------------------------------------------------------------------
// Main program
// -----------------------------------------------------------------------------

int main(int argc, char** argv)
{
    if (argc != 2 ||
        (std::string(argv[1]) != "beam" &&
         std::string(argv[1]) != "reaction")) {
        std::cerr << "Usage: " << argv[0] << " beam|reaction\n";
        return 2;
    }

    try {
        const std::string mode = argv[1];

        const BeamConfig beam = readBeam("User_Input/beam_input.txt");
        Material target = readMaterial("User_Input/target_input.dat");
        DegraderConfig degrader = readDegrader("User_Input/degrader_input.txt");
        std::mt19937_64 random(beam.seed);

        TransportModel targetModel = buildTransportModel(target);
        TransportModel degraderModel;

        if (degrader.enabled) {
            degraderModel = buildTransportModel(degrader);
        }

        BeamMC::ChargeStateCalculator chargeCalculator;

        std::ofstream output;

#ifdef BEAMMC_WITH_ROOT
        RootWriter* rootWriter = nullptr;
        if (beam.rootOutput) {
            rootWriter = new RootWriter(
                beam.rootFile,
                target.thicknessMm,
                degrader.enabled,
                degrader.thicknessMm,
                beam.energyMeV,
                beam.stepMm);
            gRootWriter = rootWriter;
        }
#else
        if (beam.rootOutput) {
            throw std::runtime_error(
                "ROOT_OUTPUT is enabled, but BeamMC was compiled without ROOT support. "
                "Compile with -DBEAMMC_WITH_ROOT and root-config.");
        }
        RootWriter* rootWriter = nullptr;
#endif

        ReactionConfig reaction;

        if (mode == "reaction") {
            reaction = readReaction("User_Input/reaction_input.txt");

            BeamMC::ExcitationPopulation excitationPopulation;
            if (!reaction.excitationPopulationFile.empty()) {
                if (reaction.neutronSeparationMeV > 0.0) {
                    excitationPopulation.setMaxExcitationMeV(
                        reaction.neutronSeparationMeV);
                }
                excitationPopulation.load(reaction.excitationPopulationFile);
            }

            if (beam.Z != reaction.beamZ || beam.A != reaction.beamA) {
                std::cerr << "Warning: beam species differ between input files\n";
            }

            // Incoming beam loses energy according to the beam-in-target table.
            target.srimFile = reaction.beamTargetSRIM;
            targetModel = buildTransportModel(target);

            // Recoils use their own stopping table in the remaining target.
            Material recoilTarget = target;
            recoilTarget.srimFile = reaction.recoilTargetSRIM;
            const TransportModel recoilTargetModel =
                buildTransportModel(recoilTarget);

            BeamMC::ReactionCrossSection reactionCrossSection;
            ReactionDepthDistribution depthDistribution;

            if (!reaction.reactionCrossSectionFile.empty()) {
                reactionCrossSection.loadTALYSFile(
                    reaction.reactionCrossSectionFile);

                depthDistribution = buildReactionDepthDistribution(
                    beam, target, targetModel, reactionCrossSection, reaction);
            }

            // In reaction mode, the degrader table must describe the recoil.
            if (degrader.enabled) {
                if (reaction.recoilDegraderSRIM.empty()) {
                    throw std::runtime_error(
                        "Reaction mode with a degrader requires "
                        "recoil_srim_degrader");
                }
                degrader.srimFile = reaction.recoilDegraderSRIM;
                degraderModel = buildTransportModel(degrader);
            }

            output.open(reaction.output);
            output << "E_recoil_MeV,theta_x_mrad,theta_y_mrad,"
                      "x_mm,y_mm,z_mm,q_exit\n";

            long long kinematicallyRejectedEvents = 0;

            for (int event = 0; event < beam.particles; ++event) {
                Particle incoming = generateInitialBeamParticle(beam, random);
                const Particle incident = incoming;

                double depthMm = depthDistribution.sample(
                    random, target, beam.stepMm);

                // Transport the incoming beam to the reaction point.
                propagate(incoming, targetModel, depthMm, beam.stepMm,
                          beam.Z, beam.A, beam.stragglingCoefficient, random);

#ifdef BEAMMC_WITH_ROOT
                if (beam.rootOutput && rootWriter) {
                    rootWriter->fillReactionDepth(depthMm, incoming.E);
                }
#endif

                // Generate the recoil at the reaction point.
                const BeamMC::ExcitationPopulation* population =
                    excitationPopulation.empty() ? nullptr : &excitationPopulation;
                Particle recoil = generateRecoil(
                    incoming, reaction, random, population);

                if (recoil.E <= 0.0) {
                    ++kinematicallyRejectedEvents;
                    continue;
                }

                {
                    const double cosThetaZ = recoil.direction.z;
                    const double pathMm =
                        (cosThetaZ > 1.0e-6) ?
                        (target.thicknessMm - depthMm) / cosThetaZ : 0.0;

                    propagate(recoil, recoilTargetModel, pathMm, beam.stepMm,
                              reaction.recoilZ, reaction.recoilA,
                              beam.stragglingCoefficient, random);
                }

                if (degrader.enabled && recoil.direction.z > 0.0) {
                    const double cosThetaZ = recoil.direction.z;
                    const double pathMm =
                        (cosThetaZ > 1.0e-6) ?
                        degrader.thicknessMm / cosThetaZ : 0.0;

                    propagate(recoil, degraderModel, pathMm, beam.stepMm,
                              reaction.recoilZ, reaction.recoilA,
                              beam.stragglingCoefficient, random);
                }

                // Determine medium from which the recoil exits.
                const Material* exitMaterial = nullptr;
                if (degrader.enabled && recoil.direction.z > 0.0) {
                    exitMaterial = &degrader;
                } else {
                    exitMaterial = &target;
                }

                int qExit = 0;
                if (exitMaterial) {
                    BeamMC::ChargeStateResult cs =
                        chargeCalculator.sample(
                            reaction.recoilZ,
                            reaction.recoilA,
                            recoil.E,
                            exitMaterial->effectiveZ,
							beam.reducedWidthW,
                            random);
                    qExit = cs.sampledQ;
                    recoil.chargeState = qExit;
                }

                output << recoil.E << ','
                       << recoil.thetaX << ','
                       << recoil.thetaY << ','
                       << recoil.position.x << ','
                       << recoil.position.y << ','
                       << recoil.position.z << ','
                       << qExit << '\n';

#ifdef BEAMMC_WITH_ROOT
                if (beam.rootOutput && rootWriter) {
                    rootWriter->fill(incident, recoil, qExit);
                    rootWriter->fillRecoilDepth(recoil.position.z, recoil.E);
                }
#endif
            }

            std::cerr << "Reaction events rejected because the sampled "
                      << "excitation was kinematically inaccessible: "
                      << kinematicallyRejectedEvents << " / "
                      << beam.particles << '\n';

        } else { // beam mode
            output.open(beam.output);
            output << "E_MeV,theta_x_mrad,theta_y_mrad,x_mm,y_mm,z_mm,q_exit\n";

            for (int event = 0; event < beam.particles; ++event) {
                Particle particle = generateInitialBeamParticle(beam, random);

#ifdef BEAMMC_WITH_ROOT
                if (beam.rootOutput && rootWriter) {
                    propagate(particle, targetModel, target.thicknessMm,
                              beam.stepMm, beam.Z, beam.A,
                              beam.stragglingCoefficient, random,
                              &recordBeamEnergyDepth);
                } else {
                    propagate(particle, targetModel, target.thicknessMm,
                              beam.stepMm, beam.Z, beam.A,
                              beam.stragglingCoefficient, random);
                }
#else
                propagate(particle, targetModel, target.thicknessMm,
                          beam.stepMm, beam.Z, beam.A,
                          beam.stragglingCoefficient, random);
#endif

                if (degrader.enabled && particle.direction.z > 0.0) {
                    const double cosThetaZ = particle.direction.z;
                    const double pathMm =
                        (cosThetaZ > 1.0e-6) ?
                        degrader.thicknessMm / cosThetaZ : 0.0;

#ifdef BEAMMC_WITH_ROOT
                    if (beam.rootOutput && rootWriter) {
                        propagate(particle, degraderModel, pathMm, beam.stepMm,
                                  beam.Z, beam.A,
                                  beam.stragglingCoefficient, random,
                                  &recordBeamEnergyDepth);
                    } else {
                        propagate(particle, degraderModel, pathMm, beam.stepMm,
                                  beam.Z, beam.A,
                                  beam.stragglingCoefficient, random);
                    }
#else
                    propagate(particle, degraderModel, pathMm, beam.stepMm,
                              beam.Z, beam.A,
                              beam.stragglingCoefficient, random);
#endif
                }

                // Determine medium from which the beam exits.
                const Material* exitMaterial = nullptr;
                if (degrader.enabled && particle.direction.z > 0.0) {
                    exitMaterial = &degrader;
                } else {
                    exitMaterial = &target;
                }

                int qExit = 0;
                if (exitMaterial) {
                    BeamMC::ChargeStateResult cs =
                        chargeCalculator.sample(
                            beam.Z,
                            beam.A,
                            particle.E,
                            exitMaterial->effectiveZ,
							beam.reducedWidthW,
                            random);
                    qExit = cs.sampledQ;
                    particle.chargeState = qExit;
                }

                output << particle.E << ','
                       << particle.thetaX << ','
                       << particle.thetaY << ','
                       << particle.position.x << ','
                       << particle.position.y << ','
                       << particle.position.z << ','
                       << qExit << '\n';

#ifdef BEAMMC_WITH_ROOT
                if (beam.rootOutput && rootWriter) {
                    rootWriter->fill(particle, particle, qExit);
                }
#endif
            }
        }

#ifdef BEAMMC_WITH_ROOT
        if (rootWriter) {
            gRootWriter = nullptr;
            delete rootWriter;
            rootWriter = nullptr;
        }
#endif

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "BeamMC error: " << ex.what() << '\n';
        return 1;
    }
}