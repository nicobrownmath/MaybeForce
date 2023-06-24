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

class SearchParams {
    public:
    int minSaveInterval;
    std::string allOutputFile = "output-all.rle";
    std::string partialsOutputFile = "output-part.rle";
    bool outputAll = true;
    bool findShuttles = true;
    unsigned matchToGen = 1;

    std::string pattern = "";
    int patternX = 0;
    int patternY = 0;
    
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
    }

    static void AddCatalyst(std::vector<Catalyst> *catalysts, const std::vector<std::string> &elems) {
        if (elems.size() <= 3) {
            std::cout << "Bad catalyst: Missing parameters" << std::endl;
            exit(0);
        }

        Catalyst catalyst;
        catalyst.state = LifeState::Parse(elems[1].c_str(), stoi(elems[2]), stoi(elems[3]));
        catalyst.locus = catalyst.state;
        
        std::vector<SymmetryTransform> transforms = CharToTransforms(elems[4].at(0));

        int index = 5;
        
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
            else {
                printf("Bad catalyst: Invalid parameter %s\n", elems[index].c_str());
                exit(0);
            }
        }

        for (auto &transform : transforms) {
            Catalyst transformedCatalyst = catalyst;
            transformedCatalyst.Transform(transform);
            transformedCatalyst.FillOutData();
            
            printf("Loaded catalyst %d\n", (int)catalysts->size());
            catalysts->push_back(std::move(transformedCatalyst));
        }
    }

    LifeState FindNewInteractionPositions(const LifeEnvelope &interactWithEnvelope, const LifeEnvelope &excludeEnvelope) {
        LifeState catalystPositions = locusOneNeighborR180.Convolve(interactWithEnvelope.twoNeighbors);
        catalystPositions |= locusTwoNeighborsR180.Convolve(interactWithEnvelope.oneNeighbor);
        catalystPositions |= locusZOIR180.Convolve(interactWithEnvelope.threeNeighbors);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        catalystPositions &= ~locusOneNeighborR180.Convolve(excludeEnvelope.twoNeighbors);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        catalystPositions &= ~locusTwoNeighborsR180.Convolve(excludeEnvelope.oneNeighbor);
        if (catalystPositions.IsEmpty()) return catalystPositions;
        catalystPositions &= ~locusZOIR180.Convolve(excludeEnvelope.threeNeighbors);
        return catalystPositions;
    }
};

class SearchData {
    public:
    unsigned generation;
    LifeState startState;
    std::map<
        LifeState,
        std::set<std::pair<std::vector<std::tuple<unsigned, unsigned, unsigned, bool>>, LifeEnvelope>>
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
                
            if (elems[0] == "outputFile") {
                params.allOutputFile = elems[1] + "-all.rle";
                params.partialsOutputFile = elems[1] + "-part.rle";
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

            if (elems[0] == "pattern") {
                params.pattern = elems[1];

                if (elems.size() > 3) {
                    params.patternX = stoi(elems[2]);
                    params.patternY = stoi(elems[3]);
                }
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

            if (elems[0] == "cat")
            {
                Catalyst::AddCatalyst(&catalysts, elems);
            }
        }
    }

    void Run() {
        mainData = GetStartData();

        for (unsigned gen = 0; gen < params.maxGen; gen++) {
            /*for (auto [stateData, catalystsDatas] : mainData.statePossibilities) {
                bool stillInteracting = false;
                LifeState exampleStartState = mainData.startState;
                LifeState exampleEndState = stateData;
                for (auto [catIndex, catX, catY, catActive] : catalystsDatas.begin()->first) {
                    exampleStartState |= catalysts[catIndex].state.Moved(catX, catY);
                    if (!catActive) {
                        exampleEndState |= catalysts[catIndex].state.Moved(catX, catY);
                    }
                    else {
                        stillInteracting = true;
                        break;
                    }
                }
                if (stillInteracting) continue; 
                exampleStartState.Print();
                exampleEndState.Print();
            }*/

            mainData = GetNextStep(mainData);
        }
        Report();

        /*for (auto [stateData, catalystsDatas] : mainData.statePossibilities) {
            LifeState exampleStartState = mainData.startState;
            LifeState exampleEndState = stateData;
            for (auto [catIndex, catX, catY, catActive] : catalystsDatas.begin()->first) {
                exampleStartState |= catalysts[catIndex].state.Moved(catX, catY);
                if (!catActive) {
                    exampleEndState |= catalysts[catIndex].state.Moved(catX, catY);
                }
            }
            exampleStartState.Print();
            exampleEndState.Print();
        }*/
    }

    SearchData GetStartData() {
        SearchData output;
        output.generation = 0;
        output.startState = LifeState::Parse(params.pattern.c_str(), params.patternX, params.patternY);

        output.statePossibilities[output.startState] = {{{}, LifeEnvelope::Envelope(output.startState)}};

        return output;
    }

    //gets the possibilities for the next step
    //TODO: Detect stabilization
    SearchData GetNextStep(const SearchData &searchData) {
        SearchData output;
        output.generation = searchData.generation + 1;
        output.startState = searchData.startState;

        mainDataProgress = 0;
        for (auto [stateData, catalystsDatas] : searchData.statePossibilities) {
            //find next-gen info and add to output.statePossibilities
            
            //find the next state of the pattern
            LifeState nextState = stateData; nextState.Step();

            //find matches for startState and descendants of startState
            std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> matchPoints = FindPartialMatches(nextState, searchData.startState, std::min(params.matchToGen, output.generation), output.generation);
            //exclude matches also present in stateData
            std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> newMatchPoints;
            for (auto [genOffset, x, y, transform] : matchPoints) {
                if (genOffset == output.generation) {
                    newMatchPoints.push_back({genOffset, x, y, transform});
                } else {
                    LifeState tryMatch = searchData.startState;
                    for (unsigned gen = 0; gen < output.generation - genOffset - 1; gen++) {
                        tryMatch.Step();
                    }
                    tryMatch.Transform(transform);
                    LifeState tryMatchBorder = tryMatch.ZOI() & ~tryMatch;
                    if (!(stateData.Contains(tryMatch, x, y) && stateData.AreDisjoint(tryMatchBorder, x, y)))
                        newMatchPoints.push_back({genOffset, x, y, transform});
                }
            }

            LifeEnvelope currentEnvelope = LifeEnvelope::Envelope(stateData);
            LifeEnvelope nextEnvelope = LifeEnvelope::Envelope(nextState);

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
                std::vector<std::tuple<unsigned, unsigned, unsigned, bool>> newDefaultCatalystsDataList = catalystsData.first;

                bool allCatsFormerlyInactive = true;
                bool allCatsCurrentlyInactive = true;
                unsigned activeCats = 0;

                //remove (deactivate) all active catalysts no longer interacting with the rest of the pattern, and add (activate) all catalysts which are now interacting with the rest of the pattern
                unsigned indexInData = 0;
                for (auto [catIndex, catX, catY, catActive] : catalystsData.first) {
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
                            newDefaultNextState &= ~catalysts[catIndex].state.Moved(catX, catY);
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
                            newDefaultNextState |= catalysts[catIndex].state.Moved(catX, catY);
                            std::get<3>(newDefaultCatalystsDataList[indexInData]) = true;

                            allCatsCurrentlyInactive = false;
                            activeCats++;
                        }
                    }
                    indexInData++;
                }

                if (!isValid) continue;

                if (newDefaultCatalystsDataList.size() < params.maxCats && activeCats < params.maxActiveCats) couldBeANewCatalyst = true;

                //if all catalysts are newly inactive, this is a solution
                if (params.outputAll && allCatsCurrentlyInactive && !allCatsFormerlyInactive) {
                    LifeState exampleStartState = searchData.startState;
                    LifeState exampleEndState = newDefaultNextState;
                    for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                        exampleStartState |= catalysts[catIndex].state.Moved(catX, catY);
                        if (!catActive)
                            exampleEndState |= catalysts[catIndex].state.Moved(catX, catY);
                    }

                    allOutputsCategoryContainer.Add(exampleStartState, exampleEndState, newDefaultNextState);
                }
                if (!newMatchPoints.empty()) {
                    LifeState exampleStartState = searchData.startState;
                    LifeState exampleEndState = newDefaultNextState;
                    for (auto [catIndex, catX, catY, catActive] : newDefaultCatalystsDataList) {
                        exampleStartState |= catalysts[catIndex].state.Moved(catX, catY);
                        if (!catActive)
                            exampleEndState |= catalysts[catIndex].state.Moved(catX, catY);
                    }

                    //TODO: Add for all members of newMatchPoints
                    partialsCategoryContainer.Add(exampleStartState, exampleEndState, newMatchPoints[0]);
                }

                newNextInfos.push_back({newDefaultNextState, {newDefaultCatalystsDataList, catalystsData.second}});
            }

            //find all possible positions for new catalysts
            if (couldBeANewCatalyst) {
                std::vector<std::tuple<unsigned, unsigned, unsigned>> possibleNewCatalysts;
                if (output.generation <= params.maxCatGen) {
                    
                    for (unsigned catIndex = 0; catIndex < catalysts.size(); catIndex++) {

                        //find new interaction positions
                        //TODO: Instead of using current envelope, can use the intersection of all possible envelopes
                        LifeState possibleCatPositions = catalysts[catIndex].FindNewInteractionPositions(nextEnvelope, currentEnvelope);

                        if (!possibleCatPositions.IsEmpty()) {
                            for (unsigned x = 0; x < 64; x++) {
                                if (possibleCatPositions.state[x] != 0) {
                                    for (unsigned y = 0; y < 64; y++) {
                                        if (possibleCatPositions.GetCell(x,y) == 1) {
                                            //possible new catalyst position
                                            
                                            //checkRecovery
                                            //TODO: Make optional
                                            bool recovers = false;
                                            LifeState testState = nextState;
                                            testState |= catalysts[catIndex].state.Moved(x, y);
                                            for (unsigned testGen = searchData.generation + 2; testGen < params.maxGen; testGen++) {
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
                                buildFromNextState |= catalysts[catIndex].state.Moved(catX, catY);
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
                //adjust envelope
                newInfo.second.second.Join(LifeEnvelope::Envelope(newInfo.first));

                //add our new data to the corresponding key
                if (output.statePossibilities.find(newInfo.first) == output.statePossibilities.end()) {
                    output.statePossibilities[newInfo.first] = {};
                }
                output.statePossibilities[newInfo.first].insert(newInfo.second);
                nextDataSize = output.statePossibilities.size();
            }

            //try reporting
            iterations++;
            mainDataProgress++;
            std::chrono::time_point<std::chrono::steady_clock> timePoint = std::chrono::steady_clock::now();
            std::chrono::duration<double> currentTime = timePoint - lastSaveTime;
            if (currentTime.count() > params.minSaveInterval) { //only save reports after intervals because it's expensive and doesn't thread well
                lastSaveTime = timePoint;
                Report();
            }
        }
        return output;
    }

    std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> FindPartialMatches(const LifeState &lifeState, const LifeState &startState, const unsigned &maxMatchingGen, const unsigned &stateGen) {
        std::vector<std::tuple<unsigned, unsigned, unsigned, SymmetryTransform>> output;
        LifeState matchState = startState;
        for (unsigned matchGen = 0; matchGen < maxMatchingGen; matchGen++) {
            LifeState matchStateBorder = matchState.ZOI() & ~matchState;
            LifeState catalystsState;
            unsigned generationOffset = stateGen - matchGen;
            CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, Identity);
            if (params.findShuttles) {
                switch (C1) {//TODO: searchData.symmetry) {
                    case C1:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 64, ReflectAcrossX);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossY);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossYeqX, -1);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossYeqNegXP1, 1);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 64, Rotate90);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 64, Rotate180OddBoth);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 64, Rotate270);
                        break;
                    case D2AcrossX:
                    case D2AcrossXEven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossY);
                        break;
                    case D2AcrossY:
                    case D2AcrossYEven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 64, ReflectAcrossX);
                        break;
                    case D2diagodd:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossYeqNegXP1, 1);
                        break;
                    case D2negdiagodd:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 64, 1, ReflectAcrossYeqX, -1);
                        break;
                    case C2:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, Rotate90);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossX);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossYeqX);
                        break;
                    case C2even:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, Rotate90Even);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossXEven);
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossYeqX);
                        break;
                    case C2verticaleven:
                    case C4even:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossXEven);
                        break;
                    case C2horizontaleven:
                    case C4:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, ReflectAcrossX);
                        break;
                    case D4:
                    case D4diag:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, Rotate90);
                        break;
                    case D4even:
                    case D4diageven:
                        CheckMatches(&output, lifeState, matchState, matchStateBorder, catalystsState, generationOffset, 1, 1, Rotate90Even);
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

    void Report() {
        std::chrono::duration<double> timePassed = std::chrono::steady_clock::now() - begin;
        printf("\nSaving output at %f seconds...\n", timePassed.count());
        
        //rule string
        std::string ruleString = "x = 0, y = 0, rule = B3/S23\n";

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

        std::chrono::time_point<std::chrono::steady_clock> currentTimePoint = std::chrono::steady_clock::now();
        std::chrono::duration<double> currentTime = currentTimePoint - begin;
        printf("Output saved, took %f seconds:\n  Total Iterations: %llu\n  Iterations per second: %f\n  Generation: %u\n  Main data progress: %u/%llu\n  Next data size: %llu\n  Solution types found: %d/%u\n  Partial types found: %d/%u\n",
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
            partialsCategoryContainer.size
            );
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