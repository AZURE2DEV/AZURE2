#include "ParameterLabel.h"

#include <cmath>
#include <sstream>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "JGroup.h"
#include "PPair.h"

namespace AZURELabel {
namespace {

// Light nuclei only; R-matrix work does not go far up the chart, and an unknown
// charge degrades to "Z=<n>" rather than guessing.
const char* kElements[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr"};
const int kNumElements = (int)(sizeof(kElements) / sizeof(kElements[0]));

}  // namespace

std::string Spin(double j) {
    std::ostringstream out;
    const double twice = 2.0 * j;
    const long twiceRounded = std::lround(twice);

    if (std::fabs(twice - double(twiceRounded)) > 1.e-6) {
        // Not a multiple of 1/2; print it as-is rather than inventing a fraction.
        out << j;
        return out.str();
    }

    if (twiceRounded % 2 == 0) {
        out << (twiceRounded / 2);
    } else {
        out << twiceRounded << "/2";
    }
    return out.str();
}

std::string JPi(double j, int parity) {
    return Spin(j) + (parity < 0 ? "-" : "+");
}

std::string Nuclide(int z, double mass) {
    const int a = (int)std::lround(mass);

    // Conventional shorthands for the light projectiles. The alpha keeps its
    // "4He" form: a bare "a" is easy to misread in a pair like "n+a".
    if (z == 0 && a == 1) return "n";
    if (z == 1 && a == 1) return "p";
    if (z == 1 && a == 2) return "d";
    if (z == 1 && a == 3) return "t";

    std::ostringstream out;
    out << a;
    if (z >= 0 && z < kNumElements) {
        out << kElements[z];
    } else {
        out << "(Z=" << z << ")";
    }
    return out.str();
}

std::string Pair(CNuc* compound, int pairNum) {
    if (!compound || pairNum < 1 || pairNum > compound->NumPairs()) {
        std::ostringstream out;
        out << "pair " << pairNum;
        return out.str();
    }

    PPair* pair = compound->GetPair(pairNum);

    // PType 10 marks a photon pair; naming its two "particles" is meaningless.
    if (pair->GetPType() == 10) return "gamma";

    return Nuclide(pair->GetZ(1), pair->GetM(1)) + "+" +
           Nuclide(pair->GetZ(2), pair->GetM(2));
}

std::string Level(JGroup* jgroup, ALevel* level, int jgroupIndex, int levelIndex) {
    std::ostringstream out;
    if (!jgroup || !level) {
        out << "level " << levelIndex << " of J-group " << jgroupIndex;
        return out.str();
    }

    out.precision(4);
    out << "J^pi=" << JPi(jgroup->GetJ(), jgroup->GetPi()) << ", E=" << std::fixed
        << level->GetE() << " MeV (J-group " << jgroupIndex << ", level "
        << levelIndex << ")";
    return out.str();
}

std::string Channel(CNuc* compound, JGroup* jgroup, int channelIndex) {
    std::ostringstream out;
    out << "channel " << channelIndex;

    if (!jgroup || channelIndex < 1 || channelIndex > jgroup->NumChannels()) {
        return out.str();
    }

    AChannel* channel = jgroup->GetChannel(channelIndex);
    const char radType = channel->GetRadType();

    out << ": ";
    if (radType == 'E' || radType == 'M') {
        // Radiative channel: the multipolarity is the informative part.
        out << radType << channel->GetL() << " gamma to "
            << Pair(compound, channel->GetPairNum());
    } else {
        out << Pair(compound, channel->GetPairNum()) << ", L=" << channel->GetL()
            << ", S=" << Spin(channel->GetS());
        if (radType != 'P') out << " (type " << radType << ")";
    }
    return out.str();
}

std::string LevelAndChannel(CNuc* compound, JGroup* jgroup, ALevel* level,
                            int jgroupIndex, int levelIndex, int channelIndex) {
    return Level(jgroup, level, jgroupIndex, levelIndex) + ", " +
           Channel(compound, jgroup, channelIndex);
}

int NumParameters(CNuc* compound, EData* data) {
    int n = 0;
    if (compound) {
        for (int j = 1; j <= compound->NumJGroups(); j++) {
            JGroup* jgroup = compound->GetJGroup(j);
            n += jgroup->NumLevels() * (1 + jgroup->NumChannels());
        }
    }
    if (data) {
        std::vector<ESegment>& segments = data->GetSegments();
        for (size_t s = 0; s < segments.size(); s++)
            if (segments[s].IsVaryNorm()) n++;
        n += (int)segments.size();
    }
    return n;
}

std::string Parameter(CNuc* compound, EData* data, int minuitIndex) {
    if (minuitIndex < 0) return std::string();

    // Same order CNuc::FillMnParams uses: per J-group, per level, one energy
    // then one width per channel.
    int index = 0;
    if (compound) {
        for (int j = 1; j <= compound->NumJGroups(); j++) {
            JGroup* jgroup = compound->GetJGroup(j);
            for (int la = 1; la <= jgroup->NumLevels(); la++) {
                ALevel* level = jgroup->GetLevel(la);
                if (index == minuitIndex)
                    return "level energy: " + Level(jgroup, level, j, la);
                index++;
                for (int ch = 1; ch <= jgroup->NumChannels(); ch++) {
                    if (index == minuitIndex)
                        return "width: " +
                               LevelAndChannel(compound, jgroup, level, j, la, ch);
                    index++;
                }
            }
        }
    }

    if (!data) return std::string();

    // Then EData::FillMnParams: norms for segments that vary one, then an energy
    // shift for every segment.
    std::vector<ESegment>& segments = data->GetSegments();
    for (size_t s = 0; s < segments.size(); s++) {
        if (!segments[s].IsVaryNorm()) continue;
        if (index == minuitIndex) {
            std::ostringstream out;
            out << "normalization of segment " << segments[s].GetSegmentKey()
                << " (" << segments[s].GetDataFile() << ")";
            return out.str();
        }
        index++;
    }
    for (size_t s = 0; s < segments.size(); s++) {
        if (index == minuitIndex) {
            std::ostringstream out;
            out << "energy shift of segment " << segments[s].GetSegmentKey()
                << " (" << segments[s].GetDataFile() << ")";
            return out.str();
        }
        index++;
    }

    return std::string();
}

}  // namespace AZURELabel
