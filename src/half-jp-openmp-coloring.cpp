#include <algorithm>
#include <iostream>
#include <unordered_set>
#include "graph.h"

class HalfJPOpenMPColorGraph : public ColorGraph {
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
      if (colors[nbor] != -1) {
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
      colors[i] = -1;
    }

    #pragma omp parallel for schedule(dynamic, 12) shared(graph, colors)
    for (int i = 0; i < numNodes; i++) {
      int color = firstAvailableColor(i, graph, colors);
      colors[i] = color;
    }

    int numWrong = 0;
    #pragma omp parallel for schedule(dynamic, 2) shared(graph, colors, numWrong)
    for (int i = 0; i < numNodes; i++) {
      color c = colors[i];
      for (const auto &nbor : graph[i]) {
        if (colors[nbor] == c && i < nbor) {
          colors[i] = -1;
	      #pragma omp atomic
	      numWrong++;
          break;
	    }
      }
    }

    while (numWrong > 0) {
      #pragma omp parallel for schedule(dynamic, 2) shared(graph, colors, numWrong)
      for (int i = 0; i < numNodes; i++) {
        if (colors[i] == -1) {
          bool colorNow = true;
          for (const auto &nbor : graph[i]) {
            if (colors[nbor] == -1 && i < nbor) {
              colorNow = false;
              break;
            }
          }
          if (colorNow) {
            colors[i] = firstAvailableColor(i, graph, colors);
            #pragma omp atomic
            numWrong--;
          }
        }
      }
    }
  }
};

std::unique_ptr<ColorGraph> createHalfJPOpenMPColorGraph() {
  return std::make_unique<HalfJPOpenMPColorGraph>();
}
