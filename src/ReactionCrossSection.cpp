#include "ReactionCrossSection.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace BeamMC {

namespace {

std::string trim(const std::string& s)
{
    std::size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return std::string();
    std::size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

} // namespace

void ReactionCrossSection::loadTALYSFile(const std::string& filename)
{
    data_.clear();

    std::ifstream input(filename.c_str());
    if (!input) {
        throw std::runtime_error(
            "Cannot open reaction cross-section file " + filename);
    }

    std::string line;
    while (std::getline(input, line)) {
        std::string stripped = trim(line);
        if (stripped.empty()) continue;
        if (stripped[0] == '#') continue;

        std::istringstream stream(stripped);
        double energyMeV = 0.0;
        double xsMb = 0.0;
        if (!(stream >> energyMeV >> xsMb)) {
            continue; // not a data line
        }

        Point p;
        p.energyMeV = energyMeV;
        p.crossSectionMb = xsMb;
        data_.push_back(p);
    }

    if (data_.empty()) {
        throw std::runtime_error(
            "No numeric data found in reaction cross-section file " + filename);
    }

    std::sort(
        data_.begin(), data_.end(),
        [](const Point& a, const Point& b) {
            return a.energyMeV < b.energyMeV;
        });
}

double ReactionCrossSection::evaluate(double energyMeV) const
{
    if (data_.empty()) {
        return 0.0;
    }

    if (energyMeV <= data_.front().energyMeV) {
        return data_.front().crossSectionMb;
    }
    if (energyMeV >= data_.back().energyMeV) {
        return data_.back().crossSectionMb;
    }

    // upper_bound: first point with energy > energyMeV
    auto it = std::upper_bound(
        data_.begin(), data_.end(), energyMeV,
        [](double E, const Point& p) { return E < p.energyMeV; });

    std::size_t idx = static_cast<std::size_t>(it - data_.begin());
    if (idx == 0) {
        return data_.front().crossSectionMb;
    }

    const Point& p1 = data_[idx - 1];
    const Point& p2 = data_[idx];

    const double dE = p2.energyMeV - p1.energyMeV;
    if (dE <= 0.0) {
        return p1.crossSectionMb;
    }

    const double t = (energyMeV - p1.energyMeV) / dE;
    return p1.crossSectionMb + t * (p2.crossSectionMb - p1.crossSectionMb);
}

} // namespace BeamMC