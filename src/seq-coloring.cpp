#include <unordered_set>
#include "graph.h"

class SeqColorGraph : public ColorGraph {
public:
  void buildGraph(std::vector<graphNode> &nodes, std::vector<std::pair<int, int>> &pairs) {
    for (auto &node : nodes) {
      graph[node] = {};
    }
  
    for (auto &pair : pairs) {
      graph[pair.first].push_back(pair.second);
      graph[pair.second].push_back(pair.first);
    }
  }

  int firstAvailableColor(int node, std::unordered_map<graphNode, color> &colors) {
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

  void colorGraph(std::unordered_map<graphNode, color> &colors) {
    int graphSize = (int) graph.size();
    for (const auto &node : graph) {
      int color = firstAvailableColor(node.first, colors);
      colors[node.first] = color;
    }
  }
};

int main() {
}
