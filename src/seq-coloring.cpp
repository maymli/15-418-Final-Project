#include <unordered_set>
#include "graph.h"


class SeqColorGraph : public ColorGraph {
public:
  void buildGraph(std::vector<graphNode> &nodes, std::vector<std::pair<int, int>> &pairs,
                  std::unordered_map<graphNode, std::vector<graphNode>> &graph) {
    for (auto &node : nodes) {
      graph[node] = {};
    }
  
    for (auto &pair : pairs) {
      graph[pair.first].push_back(pair.second);
      graph[pair.second].push_back(pair.first);
    }
  }

  int firstAvailableColor(int node, std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                          std::unordered_map<graphNode, color> &colors) {
    std::unordered_set<int> usedColors;
    for (const auto &nbor : graph[node]) {
      if (colors.count(nbor) > 0) {
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
    int numNodes = (int) graph.size();
    for (int i = 0; i < numNodes; i++) {
      int color = firstAvailableColor(i, graph, colors);
      colors[i] = color;
    }
  }
};

std::unique_ptr<ColorGraph> createSeqColorGraph() {
  return std::make_unique<SeqColorGraph>();
}
