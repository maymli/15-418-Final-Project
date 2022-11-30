#include <algorithm>
#include <unordered_set>
#include "graph.h"

class OpenMPColorGraph : public ColorGraph {
public:
  void buildGraph(std::vector<graphNode> &nodes, std::vector<std::pair<int, int>> &pairs,
                  std::unordered_map<graphNode, std::vector<graphNode>> &graph) {
    // note: I don't think you can actually parallelize this part?
    // #pragma omp parallel for shared(nodes, graph)
    for (auto &node : nodes) {
      graph[node] = {};
    }
  
    size_t numPairs = pairs.size();

    // #pragma omp parallel for schedule(static, 1) shared(pairs, graph)
    for (size_t i = 0; i < numPairs; i++) {
      int first = pairs[i].first;
      int second = pairs[i].second;

      graph[first].push_back(second);

      graph[second].push_back(first);
    }
  }

  int firstAvailableColor(int node, std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                          std::unordered_map<graphNode, color> &colors) {
    std::unordered_set<int> usedColors;
    for (const auto &nbor : graph[node]) {
      if (nbor < node && colors.count(nbor) > 0) {
        usedColors.insert(colors[nbor]);
      }
    }

    int minColor = 0;
    while(true) {
      if (usedColors.find(minColor) == usedColors.end()) {
        return minColor;
      }
      minColor++;
    }
  }

  void colorGraph(std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                  std::unordered_map<graphNode, color> &colors) {
    // TODO: in order to parallelize this, I think we just need to make sure that
    // all the variables are shared, but with the current method, we'd probably 
    // run into issues with the used color set if the colors are being updated in parallel
    // and it's not recorded in the set
    // #pragma omp parallel for shared(graph, colors)
    // also we'll need to change this into a vector to iterat through with pragma
   
    int numNodes = (int) graph.size(); 
    for (int i = 0; i < numNodes; i++) {
      colors[i] = -1;
    }

    #pragma omp parallel for schedule(dynamic, 12)
    for (int i = 0; i < numNodes; i++) { // take advantage of nodes that are [0, numNodes)
      int color = firstAvailableColor(i, graph, colors);
      colors[i] = color;
    }

    int numColors = 0;
    for (int i = 0; i < numNodes; i++) {
      numColors = std::max(numColors, colors[i] + 1);
    }

    #pragma omp parallel for shared(graph, colors, numColors)
    for (int i = 0 ; i < numNodes; i++) {
      int color = colors[i];
      for (auto &nbor : graph[i]) {
        if (color == colors[nbor]) {
          #pragma omp atomic capture
          colors[i] = numColors++;
          break;
        }
      }
    }
  }
};

std::unique_ptr<ColorGraph> createOpenMPColorGraph() {
  return std::make_unique<OpenMPColorGraph>();
}
