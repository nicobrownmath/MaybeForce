#include "LifeAPI.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <thread>
#include <mutex>
#include <stack>
#include <queue>
#include <chrono>
#include <map>
#include <set>
#include <numeric>

#define NUM_TRANSFORMS 16
#define NUM_SYMMETRIES 22

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}
//from catforce
std::string GetRLE(const std::vector<std::vector<bool>> &life2d) {
  if (life2d.empty())
    return "";

  if (life2d[0].empty())
    return "";

  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < life2d[0].size(); j++) {
    bool last_val = life2d[0][j];
    unsigned run_count = 0;

    for (const auto & i : life2d) {
      bool val = i[j];

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == !last_val) {
        if (run_count > 1)
          result << run_count;

        if (last_val == 1)
          result << "o";
        else
          result << "b";

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val) {
      if (run_count > 1)
        result << run_count;

      result << "o";

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}
//from catforce
std::string GetRLE(const LifeState &s) {
  std::vector<std::vector<bool>> vec(N, std::vector<bool>(N));

  for (unsigned j = 0; j < N; j++)
    for (unsigned i = 0; i < N; i++)
      vec[i][j] = s.GetCell(i - 32, j - 32) == 1;

  return GetRLE(vec);
}
std::vector<SymmetryTransform> SymmetryGroupFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {Identity};
  case StaticSymmetry::D2AcrossX:
    return {Identity, ReflectAcrossX};
  case StaticSymmetry::D2AcrossXEven:
    return {Identity, ReflectAcrossXEven};
  case StaticSymmetry::D2AcrossY:
    return {Identity, ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {Identity, ReflectAcrossYEven};
  case StaticSymmetry::D2diagodd:
    return {Identity, ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {Identity, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Identity, Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Identity, Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Identity, Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Identity, Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Identity, Rotate90, Rotate270, Rotate180OddBoth};
  case StaticSymmetry::C4even:
    return {Identity, Rotate90Even, Rotate270Even, Rotate180EvenBoth};
  case StaticSymmetry::D4:
    return {Identity, ReflectAcrossX, ReflectAcrossY, Rotate180OddBoth};
  case StaticSymmetry::D4even:
    return {Identity, ReflectAcrossXEven, ReflectAcrossYEven, Rotate180EvenBoth};
  case StaticSymmetry::D4horizontaleven:
    return {Identity, ReflectAcrossX, ReflectAcrossYEven, Rotate180EvenHorizontal};
  case StaticSymmetry::D4verticaleven:
    return {Identity, ReflectAcrossXEven, ReflectAcrossY, Rotate180EvenVertical};
  case StaticSymmetry::D4diag:
    return {Identity, Rotate180OddBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {Identity, Rotate180EvenBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Identity, ReflectAcrossX, ReflectAcrossY, Rotate90, Rotate270, Rotate180OddBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D8even:
    return {Identity, ReflectAcrossXEven, ReflectAcrossYEven, Rotate90Even, Rotate270Even, Rotate180EvenBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  default: //D2negdiagevenSubgroupsOnly
    return {Identity, ReflectAcrossYeqNegX};
  }
}
std::vector<SymmetryTransform> CharToTransforms(char ch) {
  switch (ch) {
  case '.':
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  case '|':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossY);
  case '-':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossX);
  case '\\':
    return SymmetryGroupFromEnum(StaticSymmetry::D2diagodd);
  case '/':
    return SymmetryGroupFromEnum(StaticSymmetry::D2negdiagodd);
  case '+':
  case '@':
    return SymmetryGroupFromEnum(StaticSymmetry::C4);
  case 'x':
    return {Identity, Rotate90, ReflectAcrossX, ReflectAcrossYeqX};
  case '*':
    return SymmetryGroupFromEnum(StaticSymmetry::D8);
  default:
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  }
}
StaticSymmetry SymmetryFromString(const std::string &name) {
  std::string start = name.substr(0, 2);
  std::string rest = name.substr(2);
  if (start == "D2") {
    if (rest == "-" or rest == "vertical") {
      return StaticSymmetry::D2AcrossX;
    } else if (rest == "-even" or rest == "verticaleven") {
      return StaticSymmetry::D2AcrossXEven;
    } else if (rest == "|" or rest == "horizontal") {
      return StaticSymmetry::D2AcrossY;
    } else if (rest == "|even" or rest == "horizontaleven") {
      return StaticSymmetry::D2AcrossYEven;
    } else if (rest == "/" or rest == "/odd") {
      return StaticSymmetry::D2negdiagodd;
    } else if (rest == "\\" or rest == "\\odd") {
      return StaticSymmetry::D2diagodd;
    }
  } else if (start == "C2") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C2;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C2even;
    } else if (rest == "horizontaleven" or rest == "|even") {
      return StaticSymmetry::C2horizontaleven;
    } else if (rest == "verticaleven" or rest == "-even" or rest == "_2") {
      return StaticSymmetry::C2verticaleven;
    }
  } else if (start == "C4") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C4;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C4even;
    }
  } else if (start == "D4") {
    std::string evenOddInfo = rest.substr(1);
    if (rest[0] == '+' or (rest.size() > 1 and rest[1] == '+')) {
      if (evenOddInfo == "" or rest == "_+1") {
        return StaticSymmetry::D4;
      } else if (evenOddInfo == "even" or rest == "_+4") {
        return StaticSymmetry::D4even;
      } else if (evenOddInfo == "verticaleven" or evenOddInfo == "-even" or
                 rest == "_+2") {
        return StaticSymmetry::D4verticaleven;
      } else if (evenOddInfo == "horizontaleven" or evenOddInfo == "|even") {
        return StaticSymmetry::D4horizontaleven;
      }
    } else if (rest[0] == 'x' or (rest.size() > 1 and rest[1] == 'x')) {
      if (evenOddInfo == "" or rest == "_x1") {
        return StaticSymmetry::D4diag;
      } else if (evenOddInfo == "even" or rest == "_x4") {
        return StaticSymmetry::D4diageven;
      }
    }
  } else if (start == "D8") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::D8;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::D8even;
    }
  }
  return StaticSymmetry::C1;
}
std::vector<SymmetryTransform> SymmetryChainFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {};
  case StaticSymmetry::D2AcrossX:
    return {ReflectAcrossX};
  case StaticSymmetry::D2AcrossXEven:
    return {ReflectAcrossXEven};
  case StaticSymmetry::D2AcrossY:
    return {ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {ReflectAcrossYEven};
  case StaticSymmetry::D2diagodd:
    return {ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Rotate90, Rotate180OddBoth};
  case StaticSymmetry::C4even:
    return {Rotate90Even, Rotate180EvenBoth};
  case StaticSymmetry::D4:
    return {ReflectAcrossX, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {ReflectAcrossXEven, ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {ReflectAcrossYEven, ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {ReflectAcrossXEven, ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {ReflectAcrossX, ReflectAcrossY, Rotate90};
  case StaticSymmetry::D8even:
    return {ReflectAcrossXEven, ReflectAcrossYEven, Rotate90Even};
  default: //D2negdiagevenSubgroupsOnly
    return {ReflectAcrossYeqNegX};
  }
}

class SearchParams {
    public:
    unsigned threads = 1;
    int minSaveInterval = 60;
    std::string allOutputFile = "output-all.rle";
    std::string partialsOutputFile = "output-part.rle";
    std::string shuttlesOutputFile = "output-shuttle.rle";
    std::string oscsOutputFile = "output-osc.rle";
    bool outputAll = true;
    bool findShuttles = true;
    unsigned matchToGen = 1;
    unsigned maxOutputsPerRow = ~0U;

    std::string pattern = "";
    int patternX = 0;
    int patternY = 0;
    bool branchUseUniqueRepresentative = false;
    StaticSymmetry symmetry;
    
    unsigned maxGen = 0;
    unsigned maxCatGen = 0;

    unsigned maxCats = 0;
    unsigned maxActiveCats = ~0U;
};

class LifeEnvelope {
    public:
    LifeState oneNeighbor;
    LifeState twoNeighbors;
    LifeState threeNeighbors;
    LifeState anyNeighbors;

    bool operator==(const LifeEnvelope &b) const {
        return oneNeighbor == b.oneNeighbor && twoNeighbors == b.twoNeighbors && threeNeighbors == b.threeNeighbors && anyNeighbors == b.anyNeighbors;
    }
    friend bool operator <(const LifeEnvelope& a, const LifeEnvelope& b) {
        if (a.oneNeighbor < b.oneNeighbor) return true;
            else if (a.oneNeighbor == b.oneNeighbor) {
                if (a.twoNeighbors < b.twoNeighbors) return true;
                else if (a.twoNeighbors == b.twoNeighbors) {
                    if (a.threeNeighbors < b.threeNeighbors) return true;
                    else if (a.threeNeighbors == b.threeNeighbors) {
                        return a.anyNeighbors < b.anyNeighbors;
                    }
                    return false;
                }
                return false;
            }
        return false;
    }

    void Print() const {
        for (int j = 0; j < 64; j++) {
            for (int i = 0; i < N; i++) {
                if (anyNeighbors.GetCell(i - 32, j - 32) == 0) {
                    printf(".");
                } else if (oneNeighbor.GetCell(i - 32, j - 32) == 1) {
                    printf("1");
                } else if (twoNeighbors.GetCell(i - 32, j - 32) == 1) {
                    printf("2");
                } else if (threeNeighbors.GetCell(i - 32, j - 32) == 1) {
                    printf("3");
                } else
                    printf("O");
            }
            printf("\n");
        }

        printf("\n\n\n\n\n\n");
        }

    static inline LifeEnvelope Envelope(const LifeState &lifeState) {
        LifeEnvelope output;
        output.oneNeighbor = lifeState.OneNeighbor();
        output.twoNeighbors = lifeState.TwoNeighbors();
        output.threeNeighbors = lifeState.ThreeNeighbors();
        output.anyNeighbors = lifeState.ZOI();
        return output;
    }

    void Join(const LifeEnvelope &newEnvelope) {
        oneNeighbor.Join(newEnvelope.oneNeighbor);
        twoNeighbors.Join(newEnvelope.twoNeighbors);
        threeNeighbors.Join(newEnvelope.threeNeighbors);
        anyNeighbors.Join(newEnvelope.anyNeighbors);
    }

    void Add(const LifeEnvelope &newEnvelope) {
        threeNeighbors = (threeNeighbors & ~newEnvelope.anyNeighbors) | (twoNeighbors & newEnvelope.oneNeighbor) | (oneNeighbor & newEnvelope.twoNeighbors) | (~anyNeighbors & newEnvelope.threeNeighbors);
        twoNeighbors = (twoNeighbors & ~newEnvelope.anyNeighbors) | (oneNeighbor & newEnvelope.oneNeighbor) | (~anyNeighbors & newEnvelope.twoNeighbors);
        oneNeighbor = (oneNeighbor & ~newEnvelope.anyNeighbors) | (~anyNeighbors & newEnvelope.oneNeighbor);
        anyNeighbors = anyNeighbors | newEnvelope.anyNeighbors;
    }

    void AddJoin(const LifeEnvelope &newEnvelope) {
        Add(newEnvelope);
        Join(newEnvelope);
    }

    void Inverse() {
        oneNeighbor.Inverse();
        twoNeighbors.Inverse();
        threeNeighbors.Inverse();
        anyNeighbors.Inverse();
    }
    LifeEnvelope& operator&=(const LifeEnvelope &other) {
        oneNeighbor &= other.oneNeighbor;
        twoNeighbors &= other.twoNeighbors;
        threeNeighbors &= other.threeNeighbors;
        anyNeighbors &= other.anyNeighbors;
        return *this;
    }

    LifeEnvelope Moved(int x, int y) const {
        LifeEnvelope output = *this;
        output.oneNeighbor.Move(x,y);
        output.twoNeighbors.Move(x,y);
        output.threeNeighbors.Move(x,y);
        output.anyNeighbors.Move(x,y);
        return output;
    }

    //Note: these only work for interaction detection when patterns are *not* touching
    bool Intersects(const LifeEnvelope &otherEnvelope) const {
        return oneNeighbor.Intersects(otherEnvelope.twoNeighbors | otherEnvelope.threeNeighbors) || twoNeighbors.Intersects(otherEnvelope.oneNeighbor | otherEnvelope.threeNeighbors) || threeNeighbors.Intersects(otherEnvelope.anyNeighbors) || anyNeighbors.Intersects(otherEnvelope.threeNeighbors);
    }

    bool Intersects(const LifeEnvelope &otherEnvelope, int dX, int dY) const {
        return oneNeighbor.Intersects(otherEnvelope.twoNeighbors | otherEnvelope.threeNeighbors, dX, dY) || twoNeighbors.Intersects(otherEnvelope.oneNeighbor | otherEnvelope.threeNeighbors, dX, dY) || threeNeighbors.Intersects(otherEnvelope.anyNeighbors, dX, dY) || anyNeighbors.Intersects(otherEnvelope.threeNeighbors, dX, dY);
    }
};

template <>
struct std::hash<LifeState>
{
    std::size_t operator()(const LifeState& k) const
    {
        std::size_t output = 0;
        for (unsigned i = 0; i < N; i++) {
            output ^= hash<uint64_t>{}(k.state[i]) + 0x9e3779b9 + (output<<6) + (output>>2);
        }
        return output;
    }
};

template <class T>
class CategoryContainer {
    public:
    std::map<T, std::vector<LifeState>> categories;
    std::vector<T> keys;

    bool hasBeenUpdated = false;
    unsigned size = 0;

    unsigned maxOutputsPerRow = ~0U;

    bool Add(const LifeState &startState, const LifeState &finalState, const T &key) {
        if (categories.find(key) == categories.end()) keys.push_back(key);
        if (categories[key].size() < maxOutputsPerRow * 2)
        {
            categories[key].push_back(startState);
            categories[key].push_back(finalState);
            hasBeenUpdated = true;
            size++;
            return true;
        }
        return false;
    }

    std::string CategoriesRLE() {
        std::stringstream output;
        for (const T &key : keys) {
            std::vector<LifeState> results = categories[key];

            const unsigned Dist = 36 + 64;

            unsigned howmany = results.size();

            unsigned width = Dist * howmany;
            unsigned height = Dist;

            std::vector<std::vector<bool>> vec(width, std::vector<bool>(height));

            for (unsigned l = 0; l < howmany; l++)
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                    vec[Dist * l + i][j] = results[l].GetCell(i - 32, j - 32) == 1;

            output << GetRLE(vec);
            output << "36$";
        }

        return output.str();
    }
};

class Catalyst {
    public:
    LifeState state;
    LifeState locus;
    
    bool firstUnderTransform[NUM_TRANSFORMS];
    LifeState fundamentalDomains[NUM_SYMMETRIES];
    LifeState symmetryExclusions[NUM_SYMMETRIES];

    LifeState locusOneNeighborR180;
    LifeState locusTwoNeighborsR180;
    LifeState locusZOIR180;
    LifeEnvelope envelope;

    std::vector<std::pair<LifeState, LifeState>> forbidden;
    std::vector<LifeState> required;
    std::vector<LifeState> antirequired;
    std::vector<std::pair<LifeState, LifeState>> anyrequired; //consist of the required and antirequired parts

    //TODO: Could probably cache some ZOIs here
    bool SeparatedFrom(const int &catX, const int &catY, const LifeState &nextState) const {
        if (!nextState.Contains(state, catX, catY)) return false;
        if ((nextState ^ state.Moved(catX, catY)).ZOI().Intersects(state, catX, catY)) return false;
        return !LifeEnvelope::Envelope(nextState ^ state.Moved(catX, catY)).Intersects(envelope, catX, catY);
    }
    
    bool InteractsWith(const int &catX, const int &catY, const LifeState &nextState, const LifeEnvelope &nextEnvelope) const {
        return nextState.ZOI().Intersects(state, catX, catY) || nextEnvelope.Intersects(envelope, catX, catY);
    }

    bool Invalid(const LifeState &testState, const int &x, const int &y) const {
        for (auto &requiredState : required) {
            LifeState missingCells = requiredState;
            missingCells.Move(x, y);
            missingCells.Copy(testState, ANDNOT);
            if (!missingCells.IsEmpty()) return true;
        }
        for (auto &antirequiredState : antirequired) {
            LifeState presentState = antirequiredState;
            presentState.Move(x, y);
            presentState.Copy(testState, AND);
            if (!presentState.IsEmpty()) return true;
        }
        if (!anyrequired.empty())
        {
            bool anyRequirementMet = false;
            for (auto [requiredState, antirequiredState] : anyrequired) {
                LifeState errorState = requiredState;
                errorState.Move(x, y);
                errorState.Copy(testState, ANDNOT);
                LifeState antirequiredError = antirequiredState;
                antirequiredError.Move(x, y);
                antirequiredError.Copy(testState, AND);
                errorState.Join(std::move(antirequiredError));

                if (errorState.IsEmpty()) {
                    anyRequirementMet = true;
                    break;
                }
            }
            if (!anyRequirementMet) return true;
        }
        return false;
    }

    bool Forbidden(const LifeState &testState, const int &x, const int &y) const {
        for (auto [forbiddenState, forbiddenStateBorder] : forbidden) {
            LifeState forbiddenError = forbiddenState;
            forbiddenError.Move(x, y);
            forbiddenError.Copy(testState, ANDNOT);
            LifeState forbiddenBorderError = forbiddenStateBorder;
            forbiddenBorderError.Move(x, y);
            forbiddenBorderError.Copy(testState, AND);
            if (forbiddenError.IsEmpty() && forbiddenBorderError.IsEmpty()) {
                return true;
            }
        }
        return false;
    }

    void Transform(SymmetryTransform transform) {
        state.Transform(transform);
        locus.Transform(transform);
        
        for (int i = 0; i < (int)forbidden.size(); i++) {
            forbidden[i].first.Transform(transform);
            forbidden[i].second.Transform(transform);
        }
        for (int i = 0; i < (int)required.size(); i++) {
            required[i].Transform(transform);
        }
        for (int i = 0; i < (int)antirequired.size(); i++) {
            antirequired[i].Transform(transform);
        }
        for (int i = 0; i < (int)anyrequired.size(); i++) {
            anyrequired[i].first.Transform(transform);
            anyrequired[i].second.Transform(transform);
        }
    }
    void FillOutData() {
        envelope = LifeEnvelope::Envelope(state);

        locus.Copy(state, AND);

        locusZOIR180 = locus.ZOI();
        locusZOIR180.Transform(Rotate180OddBoth);
        locusOneNeighborR180 = locus.OneNeighbor();
        locusOneNeighborR180.Transform(Rotate180OddBoth);
        locusTwoNeighborsR180 = locus.TwoNeighbors();
        locusTwoNeighborsR180.Transform(Rotate180OddBoth);

        for (int i = 0; i < (int)required.size(); i++) {
            required[i].Copy(state, AND);
        }
        for (int i = 0; i < (int)antirequired.size(); i++) {
            antirequired[i].Copy(state, ANDNOT);
        }
        for (int i = 0; i < (int)anyrequired.size(); i++) {
            anyrequired[i].first.Copy(state, AND);
            anyrequired[i].second.Copy(state, ANDNOT);
        }

        //TODO: These are very slow

        //fill out fundamental domains
        //registering these isn't super optimized but we only need to do it once
        //TODO: Need to make sure these work properly with periodic catalysts
        for (unsigned i = 0; i < NUM_SYMMETRIES; i++) {
            fundamentalDomains[i] = LifeState();
            std::vector<SymmetryTransform> validTransforms;
            for (auto &testTransform : SymmetryGroupFromEnum(static_cast<StaticSymmetry>(i))) {
                if (testTransform == Identity) continue;
                LifeState transformedCat = state;
                transformedCat.Transform(testTransform);
                if (!transformedCat.FindMatches(state).IsEmpty()) {
                    validTransforms.push_back(testTransform);
                }
            }
            if (!validTransforms.empty()) {
                std::set<LifeState> accountedForStates;
                for (unsigned x = 0; x < 64; x++) {
                    for (unsigned y = 0; y < 64; y++) {
                        LifeState testPlacedState = state;
                        testPlacedState.Move(x, y);
                        if (accountedForStates.find(testPlacedState) == accountedForStates.end()) {
                            fundamentalDomains[i].Set(x, y);
                            for (auto &transform : validTransforms) {
                                LifeState addState = testPlacedState;
                                addState.Transform(transform);
                                accountedForStates.insert(addState);
                            }
                        }
                    }
                }
            }
            else {
                fundamentalDomains[i].Inverse();
            }
        }
        
        //fill out symmetry exclusions
        //TODO: These need to be the same for each step of the pattern
        for (unsigned i = 0; i < NUM_SYMMETRIES; i++) {
            symmetryExclusions[i] = LifeState();
            for (int x = 0; x < 64; x++) {
                for (int y = 0; y < 64; y++) {
                    LifeState testPlacedState = state;
                    testPlacedState.Move(x, y);
                    LifeState testPlacedStateWithSymmetry = testPlacedState.GetSymChain(SymmetryChainFromEnum(static_cast<StaticSymmetry>(i)));
                    testPlacedStateWithSymmetry.Copy(testPlacedState, ANDNOT);
                    //check if the zoi from testPlacedState intersects the zoi from testPlacedStateWithSymmetry
                    //Note: This is very slightly overzealous: It doesn't allow catalysts which share a ZOI but are still stable
                    //The reason for this is that they would interfere with each other's neighbor counts
                    //TODO: However, it should be possible to allow (stable) configurations like this provided neither's locus is too close to the other
                    if (!testPlacedState.ZOI().AreDisjoint(testPlacedStateWithSymmetry.ZOI())) {
                        //invalid position
                        symmetryExclusions[i].Set(x, y);
                    }
                }
            }
        }
    }

    static void AddCatalyst(std::vector<Catalyst> *catalysts, const std::vector<std::string> &elems) {
        if (elems.size() <= 3) {
            std::cout << "Bad catalyst: Missing parameters" << std::endl;
            exit(0);
        }

        Catalyst catalyst;
        catalyst.state = LifeState::Parse(elems[1].c_str(), stoi(elems[3]), stoi(elems[4]));
        catalyst.locus = catalyst.state;

        //TODO: catalyst.recoveryTime = stoi(elems[2]);
        
        std::vector<SymmetryTransform> transforms = CharToTransforms(elems[5].at(0));

        int index = 6;
        
        while (index < (int)elems.size()) {
            if (elems[index] == "forbidden") {
                LifeState forbiddenState = LifeState::Parse(elems[index + 1].c_str(), stoi(elems[index + 2]), stoi(elems[index + 3]));
                LifeState forbiddenStateBorder = forbiddenState.ZOI();
                forbiddenStateBorder.Copy(forbiddenState, ANDNOT);
                catalyst.forbidden.push_back({forbiddenState, forbiddenStateBorder});
                index += 4;
            }
            else if (elems[index] == "required") {
                catalyst.required.push_back(LifeState::Parse(elems[index + 1].c_str(), stoi(elems[index + 2]), stoi(elems[index + 3])));
                index += 4;
            }
            else if (elems[index] == "antirequired") {
                catalyst.antirequired.push_back(LifeState::Parse(elems[index + 1].c_str(), stoi(elems[index + 2]), stoi(elems[index + 3])));
                index += 4;
            }
            else if (elems[index] == "anyrequired") {
                catalyst.anyrequired.push_back({LifeState::Parse(elems[index + 1].c_str(), stoi(elems[index + 2]), stoi(elems[index + 3])), LifeState::Parse(elems[index + 4].c_str(), stoi(elems[index + 5]), stoi(elems[index + 6]))});
                index += 7;
            }
            else if (elems[index] == "locus") {
                catalyst.locus = LifeState::Parse(elems[index + 1].c_str(), stoi(elems[index + 2]), stoi(elems[index + 3]));
                index += 4;
            }
            else if (elems[index] == "slots") {
//TODO:                catalyst.slots = stoi(elems[index + 1]);
                index += 2;
            }
            else if (elems[index] == "check-recovery") {
//TODO:                catalyst.checkRecovery = true;
                index += 1;
            }
            else if (elems[index] == "check-recovery-always") {
//TODO:                catalyst.checkRecoveryAlways = true;
                index += 1;
            }
            else {
                printf("Bad catalyst: Invalid parameter %s\n", elems[index].c_str());
                exit(0);
            }
        }

        int startingIndex = catalysts->size();
        for (auto &transform : transforms) {
            Catalyst transformedCatalyst = catalyst;
            transformedCatalyst.Transform(transform);
            transformedCatalyst.FillOutData();

            //fill out transforms
            for (int i = 0; i < NUM_TRANSFORMS; i++) {
                transformedCatalyst.firstUnderTransform[i] = true;
                for (int testCatIndex = startingIndex; testCatIndex < (int)catalysts->size(); testCatIndex++) {
                    LifeState earlierCatState = (*catalysts)[testCatIndex].state;
                    earlierCatState.Transform(static_cast<SymmetryTransform>(i));
                    if (!earlierCatState.FindMatches(transformedCatalyst.state).IsEmpty()) {
                        transformedCatalyst.firstUnderTransform[i] = false;
                        break;
                    }
                }
            }
            
            printf("Loaded catalyst %d\n", (int)catalysts->size());
            catalysts->push_back(std::move(transformedCatalyst));
        }
    }
    
    bool FirstUnderSymmetry(const std::vector<SymmetryTransform> transforms) {
        for (auto &transform : transforms) {
            if (!firstUnderTransform[static_cast<int>(transform)]) return false;
        }
        return true;
    }

    LifeState FindNewInteractionPositions(const LifeEnvelope &interactWithEnvelope, const LifeEnvelope &excludeEnvelope, StaticSymmetry patternSymmetry) {
        LifeState catalystPositions = locusOneNeighborR180.Convolve(interactWithEnvelope.twoNeighbors);
        catalystPositions |= locusTwoNeighborsR180.Convolve(interactWithEnvelope.oneNeighbor);
        catalystPositions |= locusZOIR180.Convolve(interactWithEnvelope.threeNeighbors);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        //TODO: Removals need to not be locus-based
        catalystPositions &= ~locusOneNeighborR180.Convolve(excludeEnvelope.twoNeighbors);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        catalystPositions &= ~locusTwoNeighborsR180.Convolve(excludeEnvelope.oneNeighbor);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        catalystPositions &= ~locusZOIR180.Convolve(excludeEnvelope.threeNeighbors);
        
        //TODO: make this currentSymmetry?
        if (patternSymmetry != C1)
            catalystPositions &= fundamentalDomains[static_cast<int>(patternSymmetry)];

        if (patternSymmetry != C1)
            catalystPositions &= ~symmetryExclusions[static_cast<int>(patternSymmetry)];

        return catalystPositions;
    }
};

class SearchData {
    public:
    unsigned generation;
    LifeState startState;
    std::map<
        uint64_t,
        std::set<std::vector<std::tuple<unsigned, unsigned, unsigned, bool>>>
    > statePossibilities;
        //takes state-catalystsState
        //to sets of:
        //  lists of catalyst data (ID, x, y, isActive) + extraEnvelope for this positioning
        //   A catalyst is active if:
        //     any of its cells are missing or
        //     removing it would impact the pattern's next step
        //   We do not remove active catalysts from the pattern when determining the key
};

class Searcher {
    public:
    std::chrono::time_point<std::chrono::steady_clock> begin;
    std::chrono::time_point<std::chrono::steady_clock> lastSaveTime;
    
    uint64_t iterations = 0;
    uint64_t nextDataSize = 0;
    unsigned mainDataProgress;

    SearchParams params;
    SearchData mainData;

    std::vector<Catalyst> catalysts;
    
    CategoryContainer<LifeState> allOutputsCategoryContainer;
    CategoryContainer<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> partialsCategoryContainer;
    CategoryContainer<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> shuttlesCategoryContainer;
    CategoryContainer<unsigned> oscillatorsCategoryContainer;

    std::set<uint64_t> avoidOscs;

    std::mutex reportLock;
    std::mutex categoryContainerLock;

    //initialize search
    void Init(const std::string& fname) {
        std::ifstream infile;
        infile.open(fname.c_str(), std::ifstream::in);
        if (!infile.good()) {
            std::cout << "Could not open file!" << std::endl;
            exit(1);
        }

        std::string line;
        while (std::getline(infile, line)) {
            std::vector<std::string> elems;
            split(line, ' ', elems);

            if (elems.size() < 2)
                continue;
                
            if (elems[0] == "threads") {
                params.threads = stoi(elems[1]);
            }
            if (elems[0] == "outputFile") {
                params.allOutputFile = elems[1] + "-all.rle";
                params.partialsOutputFile = elems[1] + "-part.rle";
                params.shuttlesOutputFile = elems[1] + "-shuttle.rle";
                params.oscsOutputFile = elems[1] + "-osc.rle";
            }
            if (elems[0] == "outputAll") {
                params.outputAll = (elems[1] == "true");
            }
            if (elems[0] == "findShuttles") {
                params.findShuttles = (elems[1] == "true");
            }
            if (elems[0] == "minSaveInterval") {
                params.minSaveInterval = stoi(elems[1]);
            }
            if (elems[0] == "matchToGen") {
                params.matchToGen = stoi(elems[1]);
            }
            if (elems[0] == "maxOutputsPerRow") {
                params.maxOutputsPerRow = stoi(elems[1]);
            }

            if (elems[0] == "pattern") {
                params.pattern = elems[1];

                if (elems.size() > 3) {
                    params.patternX = stoi(elems[2]);
                    params.patternY = stoi(elems[3]);
                }
            }
            if (elems[0] == "symmetry") {
                params.symmetry = SymmetryFromString(elems[1]);
            }
            if (elems[0] == "maxGen") {
                params.maxGen = stoull(elems[1]);
            }
            if (elems[0] == "maxCatGen") {
                params.maxCatGen = stoull(elems[1]);
            }
            if (elems[0] == "maxCats") {
                params.maxCats = stoull(elems[1]);
            }
            if (elems[0] == "maxActiveCats") {
                params.maxActiveCats = stoull(elems[1]);
            }
            if (elems[0] == "branchUseUniqueRepresentative") {
                params.branchUseUniqueRepresentative = (elems[1] == "true");
            }

            if (elems[0] == "cat")
            {
                Catalyst::AddCatalyst(&catalysts, elems);
            }
        }

        allOutputsCategoryContainer.maxOutputsPerRow = params.maxOutputsPerRow;
        partialsCategoryContainer.maxOutputsPerRow = params.maxOutputsPerRow;
        shuttlesCategoryContainer.maxOutputsPerRow = params.maxOutputsPerRow;
        oscillatorsCategoryContainer.maxOutputsPerRow = params.maxOutputsPerRow;
    }

    void Run() {
        mainData = GetStartData();

        for (unsigned gen = 0; gen < params.maxGen; gen++) {
            mainData = GetNextStep(std::move(mainData));
            mainDataProgress = 0;
            Report();
        }
    }

    SearchData GetStartData() {
        SearchData output;
        output.generation = 0;
        output.startState = LifeState::Parse(params.pattern.c_str(), params.patternX, params.patternY);
        output.startState.SymChain(SymmetryChainFromEnum(params.symmetry));

        output.statePossibilities[std::hash<LifeState>{}(output.startState)] = {{}};

        return output;
    }

    std::vector<std::thread> threads;
    std::mutex getNextStepOutputLock;
    std::mutex getNextStepIndexLock;
    std::mutex getNextStepSearchDataLock;

    //gets the possibilities for the next step
    SearchData GetNextStep(const SearchData &searchData) {
        SearchData output;
        output.generation = searchData.generation + 1;
        output.startState = searchData.startState;
        
        mainDataProgress = 0;

        std::map<uint64_t, std::set<std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool>>>>::const_iterator indexInSearchData = searchData.statePossibilities.begin();
        
        //initialize threads
        threads = std::vector<std::thread>();
        for (unsigned i = 0; i < params.threads; i++) {
            threads.push_back(std::thread(&Searcher::SearchThroughMap, this, &output, &searchData, &indexInSearchData));
        }

        //wait for threads to finish
        for (unsigned i = 0; i < params.threads; i++)
        {
            threads[i].join();
        }

        return output;
    }

    void SearchThroughMap(SearchData *output, const SearchData *searchData, std::map<uint64_t, std::set<std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool>>>>::const_iterator *indexInSearchData) {
        getNextStepSearchDataLock.lock();
        unsigned searchDataGeneration = searchData->generation;
        LifeState searchDataStartState = searchData->startState;
        auto searchDataEnd = searchData->statePossibilities.end();
        getNextStepSearchDataLock.unlock();

        std::vector<SymmetryTransform> symmetryChain = SymmetryChainFromEnum(params.symmetry);

        while (true) {
            getNextStepIndexLock.lock();

            if (*indexInSearchData == searchDataEnd) {
                getNextStepIndexLock.unlock();
                break;
            }

            std::set<std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool>>> catalystsDatas = (*indexInSearchData)->second;

            (*indexInSearchData)++;
            getNextStepIndexLock.unlock();

            //find next-gen info and add to output.statePossibilities
            //reconstruct current state
            LifeState currentState = searchDataStartState;
            for (auto [catIndex, catX, catY, catActive] : *(catalystsDatas.begin())) {
                currentState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
            }
            for (unsigned gen = 0; gen < searchDataGeneration; gen++) {
                currentState.Step();
            }
            for (auto [catIndex, catX, catY, catActive] : *(catalystsDatas.begin())) {
                if(!catActive) currentState &= ~catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
            }
            
            //find the next state of the pattern
            LifeState nextState = currentState; nextState.Step();
            
            //find matches for startState and descendants of startState
            //exclude partials where the pattern intersects an active catalyst
            LifeState activeCatalystsState;
            for (auto [catIndex, catX, catY, catActive] : *(catalystsDatas.begin())) {
                if (catActive) activeCatalystsState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
            }
            std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> matchPoints = FindPartialMatches(nextState, searchDataStartState, activeCatalystsState, std::min(params.matchToGen, searchDataGeneration+1), searchDataGeneration+1);
            //exclude matches also present in stateData
            std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> newMatchPoints;
            for (auto [genOffset, x, y, transform] : matchPoints) {
                if (genOffset == searchDataGeneration + 1) {
                    newMatchPoints.push_back({genOffset, x, y, transform});
                } else {
                    LifeState tryMatch = searchDataStartState;
                    for (unsigned gen = 0; gen < searchDataGeneration - genOffset; gen++) {
                        tryMatch.Step();
                    }
                    tryMatch.Transform(transform);
                    LifeState tryMatchBorder = tryMatch.ZOI() & ~tryMatch;
                    if (!(currentState.Contains(tryMatch, x, y) && currentState.AreDisjoint(tryMatchBorder, x, y)))
                        newMatchPoints.push_back({genOffset, x, y, transform});
                }
            }

            LifeEnvelope nextEnvelope = LifeEnvelope::Envelope(nextState);
            LifeEnvelope envelopeIntersection; envelopeIntersection.Inverse();

            //add new catalysts in all possible configurations
            std::vector<
                std::pair<
                    LifeState,
                    std::pair<std::vector<std::tuple<unsigned, unsigned, unsigned, bool>>, LifeEnvelope>
                >
            > newNextInfos;

            bool couldBeANewCatalyst = false;
            
            //adjust given all possible catalystsData
            for (auto &catalystsData : catalystsDatas) {
                bool isValid = true;

                LifeState newDefaultNextState = nextState;
                std::vector<std::tuple<unsigned, unsigned, unsigned, bool>> newDefaultCatalystsDataList = catalystsData;

                bool allCatsFormerlyInactive = true;
                bool allCatsCurrentlyInactive = true;
                unsigned activeCats = 0;

                //remove (deactivate) all active catalysts no longer interacting with the rest of the pattern, and add (activate) all catalysts which are now interacting with the rest of the pattern
                unsigned indexInData = 0;
                for (auto [catIndex, catX, catY, catActive] : catalystsData) {
                    if (catActive) {
                        //validity checking
                        //TODO: Optional required recovery time
                        if (catalysts[catIndex].Invalid(nextState, catX, catY)) {
                            isValid = false;
                            break;
                        }

                        allCatsFormerlyInactive = false;

                        //cat is active previously
                        if (catalysts[catIndex].SeparatedFrom(catX, catY, nextState)) {
                            //need to deactivate this
                            newDefaultNextState &= ~catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                            std::get<3>(newDefaultCatalystsDataList[indexInData]) = false;
                        }
                        else {
                            allCatsCurrentlyInactive = false;
                            activeCats++;
                        }
                    }
                    else {
                        //cat is inactive previously
                        if (catalysts[catIndex].InteractsWith(catX, catY, nextState, nextEnvelope)) {
                            //need to activate this
                            newDefaultNextState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                            std::get<3>(newDefaultCatalystsDataList[indexInData]) = true;

                            allCatsCurrentlyInactive = false;
                            activeCats++;
                        }
                    }
                    indexInData++;
                }

                if (!isValid) continue;

                //get envelope
                LifeState thisCurrentState = searchDataStartState;
                for (auto [catIndex, catX, catY, catActive] : catalystsData) {
                    thisCurrentState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                }
                LifeEnvelope currentEnvelope;
                for (unsigned gen = 0; gen < searchDataGeneration; gen++) {
                    currentEnvelope.Join(LifeEnvelope::Envelope(thisCurrentState));
                    thisCurrentState.Step();
                }
                currentEnvelope.Join(LifeEnvelope::Envelope(thisCurrentState));
                if (newDefaultCatalystsDataList.size() < params.maxCats && activeCats < params.maxActiveCats) {
                    couldBeANewCatalyst = true;
                    envelopeIntersection &= currentEnvelope;
                }

                //detect oscillation
                LifeState detectOscillation = searchDataStartState;
                for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                    detectOscillation |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                }
                LifeState nextStateWithCats = newDefaultNextState;
                for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                    if(!catActive) nextStateWithCats |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                }
                //detectOscillation.Print();
                //nextStateWithCats.Print();
                unsigned oscPeriod = 0;
                for (unsigned gen = 0; gen < searchDataGeneration+1; gen++) {
                    if (detectOscillation == nextStateWithCats) {
                        //this is an oscillator
                        oscPeriod = searchDataGeneration+1 - gen;
                        break;
                    }
                    detectOscillation.Step();
                }
                //detectOscillation.Print();

                //if all catalysts are newly inactive, this is a solution
                categoryContainerLock.lock();
                if (params.outputAll && allCatsCurrentlyInactive && !allCatsFormerlyInactive) {
                    LifeState exampleStartState = searchDataStartState;
                    LifeState exampleEndState = newDefaultNextState;
                    for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                        exampleStartState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                        if (!catActive)
                            exampleEndState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                    }

                    allOutputsCategoryContainer.Add(exampleStartState, exampleEndState, newDefaultNextState);
                }
                if (!newMatchPoints.empty()) {
                    LifeState exampleStartState = searchDataStartState;
                    LifeState exampleEndState = newDefaultNextState;
                    for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                        exampleStartState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                        if (!catActive)
                            exampleEndState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                    }

                    for (auto &matchPointData : newMatchPoints)
                    {
                        if (std::get<3>(matchPointData) == Identity)
                            partialsCategoryContainer.Add(exampleStartState, exampleEndState, matchPointData);
                        else
                            shuttlesCategoryContainer.Add(exampleStartState, exampleEndState, matchPointData);
                    }
                }
                if (oscPeriod > 2) {
                    if (AddOscillatorsToCollection(nextStateWithCats, oscPeriod)) {
                        LifeState exampleStartState = searchDataStartState;
                        LifeState exampleEndState = newDefaultNextState;
                        for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                            exampleStartState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                            if (!catActive)
                                exampleEndState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                        
                        }
                        oscillatorsCategoryContainer.Add(exampleStartState, exampleEndState, oscPeriod);
                    }
                }
                categoryContainerLock.unlock();

                if (oscPeriod == 0) newNextInfos.push_back({newDefaultNextState, {newDefaultCatalystsDataList, currentEnvelope}});
            }

            //find all possible positions for new catalysts
            if (couldBeANewCatalyst) {
                std::vector<std::tuple<unsigned, unsigned, unsigned>> possibleNewCatalysts;
                if (searchDataGeneration < params.maxCatGen) {
                    std::vector<SymmetryTransform> currentSymGroup = SymmetryGroupFromEnum(params.symmetry);
                    for (unsigned catIndex = 0; catIndex < catalysts.size(); catIndex++) {
                        if (catalysts[catIndex].FirstUnderSymmetry(currentSymGroup)) {
                            //find new interaction positions
                            //instead of using current envelope, can use the intersection of all possible envelopes
                            LifeState possibleCatPositions = catalysts[catIndex].FindNewInteractionPositions(nextEnvelope, envelopeIntersection, params.symmetry);

                            if (!possibleCatPositions.IsEmpty()) {
                            for (unsigned x = 0; x < 64; x++) {
                                if (possibleCatPositions.state[x] != 0) {
                                    for (unsigned y = 0; y < 64; y++) {
                                        if (possibleCatPositions.GetCell(x,y) == 1) {
                                            //possible new catalyst position

                                            //avoid forbidden positions
                                            LifeState testState = nextState;
                                            testState |= catalysts[catIndex].state.Moved(x, y).GetSymChain(symmetryChain);
                                            if (catalysts[catIndex].Forbidden(testState, x, y)) continue;
                                            
                                            //checkRecovery
                                            //TODO: Make optional
                                            bool recovers = false;
                                            for (unsigned testGen = searchDataGeneration + 2; testGen < params.maxGen; testGen++) {
                                                testState.Step();
                                                if (catalysts[catIndex].Invalid(testState, x, y)) {
                                                    break;
                                                }
                                                if (catalysts[catIndex].SeparatedFrom(x, y, testState)) {
                                                    recovers = true;
                                                    break;
                                                }
                                            }
                                            if (!recovers) continue;

                                            possibleNewCatalysts.push_back({catIndex, x, y});
                                        }
                                    }
                                }
                            }
                        }
                        }
                    }
                }

                for (auto [catIndex, catX, catY] : possibleNewCatalysts) {
                    unsigned newNextInfosLength = newNextInfos.size();
                    for (unsigned infosIndex = 0; infosIndex < newNextInfosLength; infosIndex++) {
                        LifeState buildFromNextState = newNextInfos[infosIndex].first;
                        std::vector<std::tuple<unsigned, unsigned, unsigned, bool>> buildFromCatalystsDataList = newNextInfos[infosIndex].second.first;
                        LifeEnvelope buildFromNextEnvelope = newNextInfos[infosIndex].second.second;

                        unsigned activeCats = 0;
                        for (auto [catIndex, catX, catY, catActive] : buildFromCatalystsDataList) {
                            if (catActive) activeCats++;
                        }

                        if (buildFromCatalystsDataList.size() < params.maxCats && activeCats < params.maxActiveCats) {
                            //add the new catalyst and add the configuration to the infos list if possible
                            /*printf("Thing\n");
                            buildFromNextState.Print();
                            catalysts[catIndex].envelope.Print();
                            catalysts[catIndex].envelope.Moved(catX, catY).Print();
                            buildFromNextEnvelope.Print();*/

                            if (!catalysts[catIndex].InteractsWith(catX, catY, buildFromNextState, buildFromNextEnvelope)) {
                                buildFromNextState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                                buildFromCatalystsDataList.push_back({catIndex, catX, catY, true});
                                buildFromNextEnvelope.AddJoin(catalysts[catIndex].envelope.Moved(catX, catY));

                                /*buildFromNextState.Print();*/

                                newNextInfos.push_back({buildFromNextState, {buildFromCatalystsDataList, buildFromNextEnvelope}});
                            }
                        }
                    }
                }
            }

            for (auto &newInfo : newNextInfos) {
                //TODO: checkRecoveryAlways
                //TODO: should check this doesn't prune excessively
                LifeState testState = newInfo.first;
                for (auto [catIndex, catX, catY, catActive] : newInfo.second.first) {
                    if (!catActive) testState |= catalysts[catIndex].state.Moved(catX, catY).GetSymChain(symmetryChain);
                }
                bool allCatsRecover = false;
                for (unsigned testGen = searchDataGeneration + 2; testGen < params.maxGen; testGen++) {
                    testState.Step();
                    allCatsRecover = true;
                    bool aCatFails = false;
                    for (auto [catIndex, catX, catY, catActive] : newInfo.second.first) {
                        if (catalysts[catIndex].Invalid(testState, catX, catY)) {
                            aCatFails = true;
                            allCatsRecover = false;
                            break;
                        }
                        if (!catalysts[catIndex].SeparatedFrom(catX, catY, testState)) {
                            allCatsRecover = false;
                            break;
                        }
                    }
                    if (aCatFails) break;
                    if (allCatsRecover) break;
                }
                if (!allCatsRecover) continue;

                //add our new data to the corresponding key
                uint64_t hashValue = std::hash<LifeState>{}(newInfo.first);

                getNextStepOutputLock.lock();
                if (output->statePossibilities.find(hashValue) == output->statePossibilities.end()) {
                    output->statePossibilities[hashValue] = {};
                }
                if (output->statePossibilities[hashValue].empty() || !params.branchUseUniqueRepresentative)
                    output->statePossibilities[hashValue].insert(newInfo.second.first);
                nextDataSize = output->statePossibilities.size();
                getNextStepOutputLock.unlock();
            }

            //try reporting
            reportLock.lock();
            iterations++;
            mainDataProgress++;
            std::chrono::time_point<std::chrono::steady_clock> timePoint = std::chrono::steady_clock::now();
            std::chrono::duration<double> currentTime = timePoint - lastSaveTime;
            if (currentTime.count() > params.minSaveInterval) { //only save reports after intervals because it's expensive and doesn't thread well
                lastSaveTime = timePoint;
                Report();
            }
            reportLock.unlock();
        }
    }

    std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> FindPartialMatches(const LifeState &lifeState, const LifeState &startState, const LifeState &avoidState, const unsigned &maxMatchingGen, const unsigned &stateGen) {
        std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> output;
        LifeState matchState = startState;
        for (unsigned matchGen = 0; matchGen < maxMatchingGen; matchGen++) {
            LifeState matchStateBorder = matchState.ZOI() & ~matchState;
            unsigned generationOffset = stateGen - matchGen;
            CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, Identity);
            if (params.findShuttles) {
                switch (params.symmetry) {
                    case C1:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 64, ReflectAcrossX);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossY);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossYeqX, -1);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossYeqNegXP1, 1);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 64, Rotate90);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 64, Rotate180OddBoth);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 64, Rotate270);
                        break;
                    case D2AcrossX:
                    case D2AcrossXEven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossY);
                        break;
                    case D2AcrossY:
                    case D2AcrossYEven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 64, ReflectAcrossX);
                        break;
                    case D2diagodd:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossYeqNegXP1, 1);
                        break;
                    case D2negdiagodd:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 64, 1, ReflectAcrossYeqX, -1);
                        break;
                    case C2:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, Rotate90);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossX);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossYeqX);
                        break;
                    case C2even:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, Rotate90Even);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossXEven);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossYeqX);
                        break;
                    case C2verticaleven:
                    case C4even:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossXEven);
                        break;
                    case C2horizontaleven:
                    case C4:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, ReflectAcrossX);
                        break;
                    case D4:
                    case D4diag:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, Rotate90);
                        break;
                    case D4even:
                    case D4diageven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, avoidState, generationOffset, 1, 1, Rotate90Even);
                        break;
                    default: break;
                }
            }
            matchState.Step();
        }
        return output;
    }
    
    void CheckMatches(std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> *newMatchPoints, const LifeState &testState, const LifeState &matchState, const LifeState &antiMatchState, const LifeState &avoidState, const unsigned &generation, const unsigned &maxXOffset, const unsigned &maxYOffset, const SymmetryTransform &transform, const int &diagonality = 0) {
        LifeState transformedMatchState = matchState;
        transformedMatchState.Transform(transform);
        LifeState transformedAntiMatchState = antiMatchState;
        transformedAntiMatchState.Transform(transform);

        if (maxXOffset == 64 && maxYOffset == 64) {
            //pattern-match everywhere at once if we're testing everywhere
            LifeState matches = transformedMatchState.FindMatches(testState);
            for (unsigned x = 0; x < maxXOffset; x++) {
                if (matches.state[x] != 0)
                    for (unsigned y = 0; y < maxYOffset; y++) {
                        if (matches.GetCell(x, y) == 1) {
                            if (testState.AreDisjoint(transformedAntiMatchState, x, y) && avoidState.AreDisjoint(transformedMatchState, x, y))
                                newMatchPoints->push_back({generation, x, y, transform});
                        }
                    }
            }
        } else {
            //otherwise just checking each cell individually is faster
            for (unsigned x = 0; x < maxXOffset; x++) {
                for (unsigned yIndex = 0; yIndex < maxYOffset; yIndex++) {
                    unsigned y = (yIndex + (64 + diagonality) * x) % 64;
                    if (testState.Contains(transformedMatchState, x, y) && testState.AreDisjoint(transformedAntiMatchState, x, y) && avoidState.AreDisjoint(transformedMatchState, x, y))
                        newMatchPoints->push_back({generation, x, y, transform});
                }
            }
        }
    }

    bool AddOscillatorsToCollection(const LifeState &state, const unsigned &period) {
        //add all oscillators of period > 2
        LifeState stateCopy = state;
        bool foundNewOsc = false;
        for (unsigned x = 0; x < 64; x++) {
            if (stateCopy.state[x] != 0ULL) {
                for (unsigned y = 0; y < 64; y++) {
                    if (stateCopy.GetCell(x, y) == 1) {
                        //need to flood every cell based on evolution impact
                        LifeState oscMask;
                        oscMask.SetCell(x, y, 1);
                        unsigned gen = 0;
                        unsigned lastGenUpated = 0;
                        LifeState currentState = stateCopy;
                        while (lastGenUpated + period * 2 >= gen) {
                            gen++;
                            LifeState testStateNotOsc = currentState;
                            testStateNotOsc.Copy(oscMask, ANDNOT);
                            LifeState testStateOsc = currentState;
                            testStateOsc.Copy(oscMask, AND);
                            currentState.Step();
                            testStateNotOsc.Step();
                            testStateOsc.Step();
                            LifeState testState = testStateOsc;
                            testState.Join(testStateNotOsc);
                            
                            //patterns affect each other
                            //expand oscMask
                            //this must expand to contain:
                            //  any active cell or forcing cell in the ZOI of the current oscMask
                            //  any cell which is on in either osc or not osc but not in currentState
                            //TODO: This may not quite work perfectly if multi-unices are any indication
                            LifeState newOscMaskCells = oscMask.ZOI();
                            newOscMaskCells.Copy(currentState, AND);
                            LifeState newOscMaskCellsFromErrors = testState;
                            newOscMaskCellsFromErrors.Copy(currentState, ANDNOT);
                            newOscMaskCells.Join(newOscMaskCellsFromErrors);
                            if (!oscMask.Contains(newOscMaskCells)) lastGenUpated = gen;
                            oscMask.Join(newOscMaskCells);
                        }
                        LifeState isolatedOsc = stateCopy;
                        foundNewOsc |= AddOscillatorToCollection(std::move(isolatedOsc));

                        stateCopy.Copy(std::move(oscMask), ANDNOT);
                    }
                }
            }
        }
        return foundNewOsc;
    }

    bool AddOscillatorToCollection(const LifeState &oscillator) {
        unsigned period = 0;
        LifeState copy = oscillator;
        do {
            period++;
            copy.Step();
        } while (!(copy == oscillator));
        if (period > 2) {
            LifeState standardizedState = std::get<0>(copy.StandardizedWithTransforms(C1));
            for (unsigned i = 1; i < period; i++) {
                copy.Step();

                standardizedState = std::min(standardizedState, std::get<0>(copy.StandardizedWithTransforms(C1)));
            }

            if (avoidOscs.find(std::hash<LifeState>{}(standardizedState)) == avoidOscs.end()) {
                avoidOscs.insert(std::hash<LifeState>{}(standardizedState));
                return true;
            }
            else {
                return false;
            }
        }
        else return false;
    }

    void Report() {
        std::chrono::duration<double> timePassed = std::chrono::steady_clock::now() - begin;
        printf("\nSaving output at %f seconds...\n", timePassed.count());
        
        //rule string
        std::string ruleString = "x = 0, y = 0, rule = B3/S23\n";

        categoryContainerLock.lock();
        if (allOutputsCategoryContainer.hasBeenUpdated) {
            std::ofstream allResultsFile(params.allOutputFile.c_str());
            allResultsFile << ruleString;
            allResultsFile << allOutputsCategoryContainer.CategoriesRLE();
            allResultsFile.close();
            allOutputsCategoryContainer.hasBeenUpdated = false;
        }
        if (partialsCategoryContainer.hasBeenUpdated) {
            std::ofstream allResultsFile(params.partialsOutputFile.c_str());
            allResultsFile << ruleString;
            allResultsFile << partialsCategoryContainer.CategoriesRLE();
            allResultsFile.close();
            partialsCategoryContainer.hasBeenUpdated = false;
        }
        if (shuttlesCategoryContainer.hasBeenUpdated) {
            std::ofstream allResultsFile(params.shuttlesOutputFile.c_str());
            allResultsFile << ruleString;
            allResultsFile << shuttlesCategoryContainer.CategoriesRLE();
            allResultsFile.close();
            shuttlesCategoryContainer.hasBeenUpdated = false;
        }
        if (oscillatorsCategoryContainer.hasBeenUpdated) {
            std::ofstream allResultsFile(params.oscsOutputFile.c_str());
            allResultsFile << ruleString;
            allResultsFile << oscillatorsCategoryContainer.CategoriesRLE();
            allResultsFile.close();
            oscillatorsCategoryContainer.hasBeenUpdated = false;
        }

        std::chrono::time_point<std::chrono::steady_clock> currentTimePoint = std::chrono::steady_clock::now();
        std::chrono::duration<double> currentTime = currentTimePoint - begin;
        printf("Output saved, took %f seconds:\n  Total Iterations: %llu\n  Iterations per second: %f\n  Generation: %u\n  Main data progress: %u/%llu\n  Next data size: %llu\n  Solution types found: %d/%u\n  Partial types found: %d/%u\n  Shuttle types found: %d/%u\n  Oscillator types found: %d/%u\n",
            currentTime.count() - timePassed.count(),
            iterations,
            iterations / currentTime.count(),
            mainData.generation,
            mainDataProgress, 
            mainData.statePossibilities.size(),
            nextDataSize,
            (int)allOutputsCategoryContainer.keys.size(), 
            allOutputsCategoryContainer.size,
            (int)partialsCategoryContainer.keys.size(), 
            partialsCategoryContainer.size,
            (int)shuttlesCategoryContainer.keys.size(), 
            shuttlesCategoryContainer.size,
            (int)oscillatorsCategoryContainer.keys.size(), 
            oscillatorsCategoryContainer.size
            );

        categoryContainerLock.unlock();
    }
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./LocalForce <in file>" << std::endl;
        exit(0);
    }

    Searcher searcher;
    searcher.begin = std::chrono::steady_clock::now();
    
    searcher.Init(argv[1]);

    std::chrono::duration<double> timePassed = std::chrono::steady_clock::now() - searcher.begin;
    printf("Initialization time: %f seconds\n", timePassed.count());

    searcher.Run();

    timePassed = std::chrono::steady_clock::now() - searcher.begin;
    printf("\nSearch Ended\n");
    printf("Total elapsed time: %f seconds\n", timePassed.count());
}