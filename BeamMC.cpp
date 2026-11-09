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
//   g++ -O2 -std=c++11 BeamMC.cpp -o BeamMC
//
// Run:
//   ./BeamMC beam
//   ./BeamMC reaction
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
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
    std::string energyUnit = "keV";
    std::string stoppingUnit = "MeV_cm2_g";
};

struct BeamConfig {
    int Z = 0;
    int A = 0;
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
};

struct DegraderConfig : Material {
    bool enabled = false;
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

    std::string beamTargetSRIM;
    std::string recoilTargetSRIM;
    std::string recoilDegraderSRIM;
    std::string output = "reaction_output.csv";
};

struct Point {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 1.0;
};

struct Particle {
    Point p;
    Vec3 dir;
    double E = 0.0;                    // MeV
};

// -----------------------------------------------------------------------------
// Small utilities for reading key=value or key:value input files
// -----------------------------------------------------------------------------

std::string trim(const std::string& s)
{
    const std::size_t first = s.find_first_not_of(" \t\r\n");
    const std::size_t last = s.find_last_not_of(" \t\r\n");

    if (first == std::string::npos) {
        return "";
    }
    return s.substr(first, last - first + 1);
}

std::string removeComments(std::string s)
{
    const std::size_t hash = s.find('#');
    if (hash != std::string::npos) {
        s.resize(hash);
    }

    const std::size_t semicolon = s.find(';');
    if (semicolon != std::string::npos) {
        s.resize(semicolon);
    }

    return trim(s);
}

std::map<std::string, std::string> readKeyValueFile(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Cannot open input file " + filename);
    }

    std::map<std::string, std::string> values;
    std::string line;

    while (std::getline(input, line)) {
        line = removeComments(line);
        if (line.empty()) {
            continue;
        }

        std::size_t separator = line.find('=');
        if (separator == std::string::npos) {
            separator = line.find(':');
        }

        if (separator != std::string::npos) {
            const std::string key = trim(line.substr(0, separator));
            const std::string value = trim(line.substr(separator + 1));
            values[key] = value;
        }
    }

    return values;
}

double getDouble(const std::map<std::string, std::string>& values,
                 const std::string& key,
                 double defaultValue = 0.0)
{
    const auto it = values.find(key);
    return it == values.end() ? defaultValue : std::stod(it->second);
}

int getInt(const std::map<std::string, std::string>& values,
           const std::string& key,
           int defaultValue = 0)
{
    return static_cast<int>(getDouble(values, key, defaultValue));
}

std::string getString(const std::map<std::string, std::string>& values,
                      const std::string& key,
                      const std::string& defaultValue = "")
{
    const auto it = values.find(key);
    return it == values.end() ? defaultValue : it->second;
}

// -----------------------------------------------------------------------------
// Input readers
// -----------------------------------------------------------------------------

Material readMaterial(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    Material material;

    material.thicknessMm = getDouble(values, "thickness_mm");
    material.density = getDouble(values, "density_g_cm3");
    material.effectiveZ = getDouble(
        values, "effective_Z", getDouble(values, "Z", 6.0));
    material.effectiveA = getDouble(
        values, "effective_A", getDouble(values, "A", 12.0));
    material.srimFile = getString(
        values, "stopping_table", getString(values, "srim_file"));
    material.energyUnit = getString(values, "srim_energy_unit", "keV");
    material.stoppingUnit = getString(
        values, "srim_stopping_unit", "MeV_cm2_g");

    if (material.thicknessMm <= 0.0 || material.density <= 0.0 ||
        material.srimFile.empty()) {
        throw std::runtime_error(
            "Material needs thickness, density, and stopping_table in " +
            filename);
    }

    return material;
}

BeamConfig readBeam(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    BeamConfig beam;

    beam.Z = getInt(values, "Z");
    beam.A = getInt(values, "A");
    beam.energyMeV = getDouble(values, "energy_MeV");
    beam.energySpreadMeV = getDouble(values, "energy_spread_MeV");

    beam.emittanceX = getDouble(
        values, "emittance_x_mm_mrad",
        getDouble(values, "emittance_mm_mrad"));
    beam.emittanceY = getDouble(values, "emittance_y_mm_mrad", beam.emittanceX);

    beam.sigmaXmm = getDouble(values, "beam_sigma_x_mm", 1.0);
    beam.sigmaYmm = getDouble(values, "beam_sigma_y_mm", beam.sigmaXmm);
    beam.particles = getInt(values, "particles", 10000);
    beam.stepMm = getDouble(values, "step_mm", 0.01);
    beam.seed = static_cast<unsigned long long>(getDouble(values, "seed", 0));
    beam.stragglingCoefficient = getDouble(
        values, "straggling_coefficient", 0.02);
    beam.output = getString(values, "output", beam.output);

    if (beam.Z <= 0 || beam.A <= 0 || beam.energyMeV <= 0.0) {
        throw std::runtime_error("Invalid beam input");
    }

    if (beam.sigmaXmm <= 0.0 || beam.sigmaYmm <= 0.0 ||
        beam.stepMm <= 0.0 || beam.particles <= 0) {
        throw std::runtime_error("Invalid beam width, step size, or particle count");
    }

    if (beam.seed == 0) {
        beam.seed = std::random_device()();
    }

    return beam;
}

DegraderConfig readDegrader(const std::string& filename)
{
    const auto values = readKeyValueFile(filename);
    DegraderConfig degrader;

    degrader.enabled = getInt(values, "enabled") != 0;
    if (!degrader.enabled) {
        return degrader;
    }

    degrader.thicknessMm = getDouble(values, "thickness_mm");
    degrader.density = getDouble(values, "density_g_cm3");
    degrader.effectiveZ = getDouble(values, "effective_Z", 6.0);
    degrader.effectiveA = getDouble(values, "effective_A", 12.0);
    degrader.srimFile = getString(values, "stopping_table");
    degrader.energyUnit = getString(values, "srim_energy_unit", "keV");
    degrader.stoppingUnit = getString(
        values, "srim_stopping_unit", "MeV_cm2_g");

    if (degrader.thicknessMm <= 0.0 || degrader.density <= 0.0 ||
        degrader.srimFile.empty()) {
        throw std::runtime_error(
            "Enabled degrader needs thickness, density, and stopping_table");
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

    reaction.excitationMeV = getDouble(
        values, "recoil_excitation_MeV");
    reaction.qValueMeV = getDouble(values, "q_value_MeV");

    reaction.beamTargetSRIM = getString(values, "beam_srim_target");
    reaction.recoilTargetSRIM = getString(values, "recoil_srim_target");
    reaction.recoilDegraderSRIM = getString(values, "recoil_srim_degrader");
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

struct StoppingTable {
    std::vector<double> energyMeV;
    std::vector<double> stoppingPower;
};

StoppingTable loadSRIMTable(const std::string& filename,
                            const std::string& energyUnit)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Cannot open SRIM file " + filename);
    }

    StoppingTable table;
    std::string line;

    while (std::getline(input, line)) {
        // Allow comma-, semicolon-, and tab-separated numeric tables.
        for (char& character : line) {
            if (character == ',' || character == ';' || character == '\t') {
                character = ' ';
            }
        }

        std::istringstream stream(line);
        std::vector<double> numbers;
        double value = 0.0;

        while (stream >> value) {
            numbers.push_back(value);
        }

        // Header and non-numeric lines are ignored.
        if (numbers.size() < 2 || numbers.front() <= 0.0 ||
            numbers.back() <= 0.0) {
            continue;
        }

        double energyMeV = numbers.front();
        if (energyUnit == "keV" || energyUnit == "kev") {
            energyMeV /= 1000.0;
        } else if (energyUnit == "eV" || energyUnit == "ev") {
            energyMeV /= 1.0e6;
        }

        table.energyMeV.push_back(energyMeV);
        table.stoppingPower.push_back(numbers.back());
    }

    if (table.energyMeV.size() < 2) {
        throw std::runtime_error(
            "No usable numeric SRIM table data in " + filename);
    }

    // Sort by energy and remove duplicate energy points.
    std::vector<std::size_t> order(table.energyMeV.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }

    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return table.energyMeV[a] < table.energyMeV[b];
    });

    StoppingTable sorted;
    for (const std::size_t index : order) {
        if (sorted.energyMeV.empty() ||
            table.energyMeV[index] > sorted.energyMeV.back()) {
            sorted.energyMeV.push_back(table.energyMeV[index]);
            sorted.stoppingPower.push_back(table.stoppingPower[index]);
        }
    }

    return sorted;
}

double interpolateStoppingPower(const StoppingTable& table, double energyMeV)
{
    if (energyMeV <= table.energyMeV.front()) {
        return table.stoppingPower.front();
    }
    if (energyMeV >= table.energyMeV.back()) {
        return table.stoppingPower.back();
    }

    const auto upper = std::upper_bound(
        table.energyMeV.begin(), table.energyMeV.end(), energyMeV);
    const std::size_t i = static_cast<std::size_t>(
        upper - table.energyMeV.begin() - 1);

    const double fraction =
        (energyMeV - table.energyMeV[i]) /
        (table.energyMeV[i + 1] - table.energyMeV[i]);

    return table.stoppingPower[i] + fraction *
        (table.stoppingPower[i + 1] - table.stoppingPower[i]);
}

struct TransportModel {
    Material material;
    StoppingTable table;
    double radiationLength = 1.0e30;  // g/cm^2, approximate
};

double approximateRadiationLength(double Z, double A)
{
    return 716.4 * A /
        (Z * (Z + 1.0) * std::log(287.0 / std::sqrt(Z)));
}

TransportModel buildTransportModel(const Material& material)
{
    TransportModel model;
    model.material = material;
    model.table = loadSRIMTable(material.srimFile, material.energyUnit);
    model.radiationLength = approximateRadiationLength(
        material.effectiveZ, material.effectiveA);
    return model;
}

// Convert the SRIM mass stopping power to MeV/mm.
double energyLossPerMillimetre(const TransportModel& model, double energyMeV)
{
    double stoppingPower = interpolateStoppingPower(
        model.table, std::max(energyMeV, 1.0e-9));

    // SRIM tables are commonly supplied in MeV cm^2/g or MeV cm^2/mg.
    if (model.material.stoppingUnit == "MeV_cm2_mg") {
        stoppingPower *= 1000.0;
    }

    // density [g/cm^3] × stopping power [MeV cm^2/g] gives MeV/cm.
    // Divide by 10 to convert to MeV/mm.
    return stoppingPower * model.material.density / 10.0;
}

// -----------------------------------------------------------------------------
// Vector and geometry helpers
// -----------------------------------------------------------------------------

Vec3 unitVector(const Vec3& vector)
{
    const double magnitude = std::sqrt(
        vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

    if (magnitude <= 0.0) {
        return {0.0, 0.0, 1.0};
    }

    return {
        vector.x / magnitude,
        vector.y / magnitude,
        vector.z / magnitude
    };
}

Vec3 add(const Vec3& a, const Vec3& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 scale(const Vec3& vector, double factor)
{
    return {factor * vector.x, factor * vector.y, factor * vector.z};
}

// Distance to the target surface reached by a particle emitted at depth.
double pathToTargetSurface(double depthMm,
                           double targetThicknessMm,
                           const Vec3& direction)
{
    if (direction.z > 1.0e-9) {
        return (targetThicknessMm - depthMm) / direction.z;
    }
    if (direction.z < -1.0e-9) {
        return depthMm / (-direction.z);
    }

    // A direction parallel to the target faces does not leave through a face.
    return 1.0e30;
}

// -----------------------------------------------------------------------------
// Particle generation and transport
// -----------------------------------------------------------------------------

Particle generateInitialBeamParticle(const BeamConfig& beam,
                                     std::mt19937_64& random)
{
    std::normal_distribution<double> gaussian(0.0, 1.0);
    Particle particle;

    particle.p = {
        beam.sigmaXmm * gaussian(random),
        beam.sigmaYmm * gaussian(random),
        0.0
    };

    // For a Gaussian beam, the approximate angular rms width is emittance / sigma.
    particle.dir = unitVector({
        beam.emittanceX / beam.sigmaXmm / 1000.0 * gaussian(random),
        beam.emittanceY / beam.sigmaYmm / 1000.0 * gaussian(random),
        1.0
    });

    particle.E = std::max(
        0.0, beam.energyMeV + beam.energySpreadMeV * gaussian(random));

    return particle;
}

void propagate(Particle& particle,
               const TransportModel& model,
               double lengthMm,
               double stepMm,
               double charge,
               double stragglingCoefficient,
               std::mt19937_64& random)
{
    if (lengthMm <= 0.0 || particle.E <= 0.0) {
        return;
    }

    std::normal_distribution<double> gaussian(0.0, 1.0);
    double remainingMm = lengthMm;
    const double massMeV = 12.0 * AMU_MEV;

    while (remainingMm > 1.0e-12 && particle.E > 1.0e-9) {
        const double distanceMm = std::min(stepMm, remainingMm);
        const double deterministicLoss = std::min(
            particle.E,
            std::max(0.0, energyLossPerMillimetre(model, particle.E) *
                distanceMm));

        // Simple Gaussian straggling model. The coefficient is user-configurable.
        const double sigmaLoss = std::sqrt(
            std::max(0.0, stragglingCoefficient * deterministicLoss *
                std::max(particle.E, 1.0e-6)));

        particle.E = std::max(
            0.0,
            particle.E - deterministicLoss + sigmaLoss * gaussian(random));

        particle.p.x += particle.dir.x * distanceMm;
        particle.p.y += particle.dir.y * distanceMm;
        particle.p.z += particle.dir.z * distanceMm;

        // Approximate multiple scattering using a Highland-like expression.
        const double beta = std::max(
            1.0e-5, std::sqrt(2.0 * particle.E / massMeV));
        const double momentumMeV = std::max(
            1.0e-5, std::sqrt(2.0 * massMeV * particle.E));
        const double arealDensity = model.material.density * distanceMm / 10.0;
        const double theta = 13.6 * charge / (beta * momentumMeV) *
            std::sqrt(std::max(0.0, arealDensity / model.radiationLength));

        particle.dir = unitVector({
            particle.dir.x + theta * gaussian(random) / std::sqrt(2.0),
            particle.dir.y + theta * gaussian(random) / std::sqrt(2.0),
            particle.dir.z
        });

        remainingMm -= distanceMm;
    }
}

// -----------------------------------------------------------------------------
// Two-body reaction kinematics
// -----------------------------------------------------------------------------

Particle generateRecoil(const Particle& incoming,
                        const ReactionConfig& reaction,
                        std::mt19937_64& random)
{
    Particle recoil;

    const double beamMass = reaction.beamA * AMU_MEV;
    const double targetMass = reaction.targetA * AMU_MEV;
    const double recoilMass = reaction.recoilA * AMU_MEV;
    const double ejectileMass = reaction.ejectileA * AMU_MEV;

    // Non-relativistic centre-of-mass energy and available final-state energy.
    const double centreOfMassEnergy =
        incoming.E * targetMass / (beamMass + targetMass);
    const double availableEnergy = centreOfMassEnergy + reaction.qValueMeV -
        reaction.excitationMeV;

    if (availableEnergy <= 0.0) {
        return recoil;
    }

    const double reducedFinalMass =
        ejectileMass * recoilMass / (ejectileMass + recoilMass);
    const double finalMomentum =
        std::sqrt(2.0 * reducedFinalMass * availableEnergy);
    const double initialMomentum =
        std::sqrt(2.0 * beamMass * incoming.E);

    // Isotropic emission in the centre-of-mass frame.
    const double cosTheta =
        2.0 * std::generate_canonical<double, 64>(random) - 1.0;
    const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    const double phi = 2.0 * PI *
        std::generate_canonical<double, 64>(random);

    // Construct an orthonormal basis around the incoming beam direction.
    const Vec3 zAxis = unitVector(incoming.dir);
    const Vec3 reference = std::fabs(zAxis.z) < 0.9
        ? Vec3{0.0, 0.0, 1.0}
        : Vec3{0.0, 1.0, 0.0};

    const Vec3 xAxis = unitVector({
        reference.y * zAxis.z - reference.z * zAxis.y,
        reference.z * zAxis.x - reference.x * zAxis.z,
        reference.x * zAxis.y - reference.y * zAxis.x
    });

    const Vec3 yAxis = {
        zAxis.y * xAxis.z - zAxis.z * xAxis.y,
        zAxis.z * xAxis.x - zAxis.x * xAxis.z,
        zAxis.x * xAxis.y - zAxis.y * xAxis.x
    };

    const Vec3 centreOfMassDirection = {
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta
    };

    const Vec3 ejectileMomentum = scale(
        add(add(scale(xAxis, centreOfMassDirection.x),
                scale(yAxis, centreOfMassDirection.y)),
            scale(zAxis, centreOfMassDirection.z)),
        finalMomentum);

    // Momentum conservation gives the recoil momentum in the lab frame.
    const Vec3 recoilMomentum = add(
        scale(zAxis, initialMomentum),
        scale(ejectileMomentum, -1.0));

    recoil.p = incoming.p;
    recoil.dir = unitVector(recoilMomentum);
    recoil.E = (recoilMomentum.x * recoilMomentum.x +
                recoilMomentum.y * recoilMomentum.y +
                recoilMomentum.z * recoilMomentum.z) /
               (2.0 * recoilMass);

    return recoil;
}

} // namespace

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

        const BeamConfig beam = readBeam("beam_input.txt");
        Material target = readMaterial("target_input.dat");
        DegraderConfig degrader = readDegrader("degrader_input.txt");
        std::mt19937_64 random(beam.seed);

        // Build the default target and degrader models. Reaction mode replaces
        // the relevant SRIM files with projectile-specific tables below.
        TransportModel targetModel = buildTransportModel(target);
        TransportModel degraderModel;

        if (degrader.enabled) {
            degraderModel = buildTransportModel(degrader);
        }

        std::ofstream output;
        ReactionConfig reaction;

        if (mode == "reaction") {
            reaction = readReaction("reaction_input.txt");

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
                      "x_mm,y_mm,z_mm\n";

            for (int event = 0; event < beam.particles; ++event) {
                Particle incoming = generateInitialBeamParticle(beam, random);

                // Uniform reaction depth between the front and back of target.
                double depthMm = target.thicknessMm *
                    std::generate_canonical<double, 64>(random);
                depthMm = std::max(
                    std::min(depthMm, target.thicknessMm),
                    std::min(beam.stepMm, target.thicknessMm));

                // Transport the incoming beam to the reaction point.
                propagate(incoming, targetModel, depthMm, beam.stepMm,
                          beam.Z, beam.stragglingCoefficient, random);

                // Generate the recoil at the reaction point.
                Particle recoil = generateRecoil(incoming, reaction, random);

                if (recoil.E > 0.0) {
                    // The recoil traverses the remainder of the target. The
                    // path is angle-corrected and supports forward/backward
                    // recoils exiting through either target surface.
                    const double remainingTargetPath = pathToTargetSurface(
                        depthMm, target.thicknessMm, recoil.dir);

                    propagate(recoil, recoilTargetModel, remainingTargetPath,
                              beam.stepMm, reaction.recoilZ,
                              beam.stragglingCoefficient, random);

                    // Only a forward-going recoil reaches a downstream degrader.
                    if (degrader.enabled && recoil.dir.z > 1.0e-9) {
                        propagate(recoil, degraderModel,
                                  degrader.thicknessMm / recoil.dir.z,
                                  beam.stepMm, reaction.recoilZ,
                                  beam.stragglingCoefficient, random);
                    }
                }

                output << std::setprecision(10)
                       << recoil.E << ','
                       << std::atan2(recoil.dir.x, recoil.dir.z) * 1000.0 << ','
                       << std::atan2(recoil.dir.y, recoil.dir.z) * 1000.0 << ','
                       << recoil.p.x << ','
                       << recoil.p.y << ','
                       << recoil.p.z << '\n';
            }
        } else {
            output.open(beam.output);
            output << "E_out_MeV,theta_x_mrad,theta_y_mrad,"
                      "x_mm,y_mm,z_mm\n";

            for (int event = 0; event < beam.particles; ++event) {
                Particle particle = generateInitialBeamParticle(beam, random);

                // Transport the beam through the complete target.
                propagate(particle, targetModel, target.thicknessMm,
                          beam.stepMm, beam.Z,
                          beam.stragglingCoefficient, random);

                // In beam mode, a downstream degrader uses the beam stopping
                // table specified in degrader_input.txt.
                if (degrader.enabled && particle.dir.z > 1.0e-9) {
                    propagate(particle, degraderModel,
                              degrader.thicknessMm / particle.dir.z,
                              beam.stepMm, beam.Z,
                              beam.stragglingCoefficient, random);
                }

                output << std::setprecision(10)
                       << particle.E << ','
                       << std::atan2(particle.dir.x, particle.dir.z) * 1000.0 << ','
                       << std::atan2(particle.dir.y, particle.dir.z) * 1000.0 << ','
                       << particle.p.x << ','
                       << particle.p.y << ','
                       << particle.p.z << '\n';
            }
        }

        if (!output) {
            throw std::runtime_error("Could not open output file");
        }

        std::cerr << "Completed " << beam.particles
                  << " particles using seed " << beam.seed << "\n";
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << "\n";
        return 1;
    }

    return 0;
}