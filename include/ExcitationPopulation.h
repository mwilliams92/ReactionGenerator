#ifndef BEAMMC_EXCITATION_POPULATION_H
#define BEAMMC_EXCITATION_POPULATION_H

#include <cstddef>
#include <random>
#include <string>
#include <vector>

namespace BeamMC {

struct ExcitationLevel {
    int bin;
    std::string type;
    double excitationMeV;
    double populationCrossSectionMb;
    double probability;
    double cumulativeProbability;
};

class ExcitationPopulation {
public:
    ExcitationPopulation();
    explicit ExcitationPopulation(const std::string& filename);

    void load(const std::string& filename);
    bool empty() const;
    std::size_t size() const;

    // Samples rows according to the probability column. Probabilities are
    // renormalised internally, so they need not sum exactly to one.
    double sample(std::mt19937_64& random) const;

    const std::vector<ExcitationLevel>& levels() const;

private:
    std::vector<ExcitationLevel> levels_;
    std::vector<double> cumulativeWeights_;
    double totalWeight_;
};

} // namespace BeamMC

#endif
