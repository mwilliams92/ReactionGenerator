#include "ExcitationPopulation.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace BeamMC {
namespace {

std::string trim(const std::string& text)
{
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string normaliseHeader(std::string text)
{
    text = trim(text);

    // Strip UTF-8 BOM if present.
    if (text.size() >= 3 &&
        static_cast<unsigned char>(text[0]) == 0xef &&
        static_cast<unsigned char>(text[1]) == 0xbb &&
        static_cast<unsigned char>(text[2]) == 0xbf) {
        text.erase(0, 3);
    }

    for (char& c : text) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return text;
}

std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;

    while (std::getline(stream, field, ',')) {
        fields.push_back(trim(field));
    }
    return fields;
}

std::size_t findColumn(const std::vector<std::string>& header,
                       const std::vector<std::string>& names)
{
    for (std::size_t i = 0; i < header.size(); ++i) {
        const std::string name = normaliseHeader(header[i]);
        for (const std::string& candidateRaw : names) {
            const std::string candidate = normaliseHeader(candidateRaw);
            if (name == candidate) {
                return i;
            }
        }
    }
    return header.size();
}

double readNumber(const std::vector<std::string>& fields,
                  std::size_t column,
                  const std::string& description,
                  std::size_t lineNumber)
{
    if (column >= fields.size() || fields[column].empty()) {
        throw std::runtime_error(
            "Missing " + description +
            " at CSV line " + std::to_string(lineNumber));
    }

    try {
        return std::stod(fields[column]);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "Invalid " + description +
            " at CSV line " + std::to_string(lineNumber) +
            ": " + fields[column]);
    }
}

} // namespace

ExcitationPopulation::ExcitationPopulation()
    : levels_(), cumulativeWeights_(), totalWeight_(0.0), maxExcitationMeV_(-1.0)
{
}

ExcitationPopulation::ExcitationPopulation(const std::string& filename)
    : levels_(), cumulativeWeights_(), totalWeight_(0.0), maxExcitationMeV_(-1.0)
{
    load(filename);
}

void ExcitationPopulation::load(const std::string& filename)
{

    std::ifstream input(filename.c_str());
    if (!input) {
        throw std::runtime_error(
            "Cannot open excitation population file: " + filename);
    }

    levels_.clear();
    cumulativeWeights_.clear();
    totalWeight_ = 0.0;

    std::string line;
    std::vector<std::string> header;
    std::size_t lineNumber = 0;

    // Read header (skip blank/comment lines).
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string cleaned = trim(line);
        if (cleaned.empty() || cleaned[0] == '#') continue;
        header = splitCsv(line);
        break;
    }

    if (header.empty()) {
        throw std::runtime_error(
            "Excitation population file has no CSV header: " + filename);
    }

    const std::size_t binColumn =
        findColumn(header, {"bin"});

    const std::size_t typeColumn =
        findColumn(header, {"type"});

    const std::size_t excitationColumn =
        findColumn(header, {
            "excitation in mev",
            "excitation (mev)",
            "excitation_mev",
            "excitation",
            "ex_mev"
        });

    const std::size_t crossSectionColumn =
        findColumn(header, {
            "population cross-section (in mb)",
            "population cross-section",
            "population_cross_section_mb",
            "population_mb"
        });

    const std::size_t probabilityColumn =
        findColumn(header, {
            "probability",
            "population probability"
        });

    const std::size_t cumulativeColumn =
        findColumn(header, {
            "cumulative probability",
            "cumulative_probability"
        });

    if (excitationColumn == header.size() ||
        probabilityColumn == header.size()) {
        throw std::runtime_error(
            "Excitation population CSV must contain excitation and "
            "probability columns: " + filename);
    }

    // Read data rows.
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string cleaned = trim(line);
        if (cleaned.empty() || cleaned[0] == '#') continue;

        const std::vector<std::string> fields = splitCsv(line);

        const double excitation = readNumber(
            fields, excitationColumn, "excitation", lineNumber);
        const double probability = readNumber(
            fields, probabilityColumn, "probability", lineNumber);

        if (!std::isfinite(excitation) || excitation < 0.0) {
            throw std::runtime_error(
                "Invalid excitation at CSV line " +
                std::to_string(lineNumber));
        }
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::runtime_error(
                "Invalid probability at CSV line " +
                std::to_string(lineNumber));
        }

		if (maxExcitationMeV_ > 0.0 && excitation > maxExcitationMeV_) {
        	continue;  // do not include this level
    	}

        // Skip zero-probability rows.
        if (probability == 0.0) continue;

        ExcitationLevel level;

        // Optional integer bin.
        if (binColumn < fields.size() && !fields[binColumn].empty()) {
            try {
                level.bin = static_cast<int>(std::stod(fields[binColumn]));
            } catch (...) {
                level.bin = 0;
            }
        } else {
            level.bin = 0;
        }

        // Optional type string.
        level.type = (typeColumn < fields.size())
                         ? fields[typeColumn]
                         : std::string();

        level.excitationMeV = excitation;

        // Optional population cross section.
        if (crossSectionColumn < fields.size() &&
            !fields[crossSectionColumn].empty()) {
            try {
                level.populationCrossSectionMb =
                    std::stod(fields[crossSectionColumn]);
            } catch (...) {
                level.populationCrossSectionMb = 0.0;
            }
        } else {
            level.populationCrossSectionMb = 0.0;
        }

        level.probability = probability;

        // Optional cumulative probability (not used in sampling).
        if (cumulativeColumn < fields.size() &&
            !fields[cumulativeColumn].empty()) {
            try {
                level.cumulativeProbability =
                    std::stod(fields[cumulativeColumn]);
            } catch (...) {
                level.cumulativeProbability = 0.0;
            }
        } else {
            level.cumulativeProbability = 0.0;
        }

        totalWeight_ += probability;
        levels_.push_back(level);
        cumulativeWeights_.push_back(totalWeight_);
    }

    if (levels_.empty() || totalWeight_ <= 0.0) {
        throw std::runtime_error(
            "Excitation population file contains no positive probabilities: " +
            filename);
    }
}

bool ExcitationPopulation::empty() const
{
    return levels_.empty();
}

std::size_t ExcitationPopulation::size() const
{
    return levels_.size();
}

const std::vector<ExcitationLevel>& ExcitationPopulation::levels() const
{
    return levels_;
}

double ExcitationPopulation::sample(std::mt19937_64& random) const
{
    if (levels_.empty() || totalWeight_ <= 0.0) {
        throw std::runtime_error(
            "Cannot sample an empty excitation population");
    }

    std::uniform_real_distribution<double> distribution(0.0, totalWeight_);
    const double draw = distribution(random);

    const auto selected = std::lower_bound(
        cumulativeWeights_.begin(), cumulativeWeights_.end(), draw);

    const std::size_t index = static_cast<std::size_t>(
        selected - cumulativeWeights_.begin());

    return levels_[std::min(index, levels_.size() - 1)].excitationMeV;
}

} // namespace BeamMC