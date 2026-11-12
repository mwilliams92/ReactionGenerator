#ifndef BEAMMC_REACTION_CROSS_SECTION_H
#define BEAMMC_REACTION_CROSS_SECTION_H

#include <vector>
#include <string>

namespace BeamMC {

class ReactionCrossSection {
public:
    struct Point {
        double energyMeV;
        double crossSectionMb;
    };

    ReactionCrossSection() {}

    // Load a TALYS YANDF-style file (e.g. *.tot).
    // The parser ignores all comment lines starting with '#'
    // and reads lines of the form:  E  xs
    void loadTALYSFile(const std::string& filename);

    bool empty() const { return data_.empty(); }

    // Piecewise-linear interpolation of sigma(E) in mb.
    // Returns 0.0 if there is no data.
    double evaluate(double energyMeV) const;

private:
    std::vector<Point> data_;
};

} // namespace BeamMC

#endif // BEAMMC_REACTION_CROSS_SECTION_H