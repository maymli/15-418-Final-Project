#include "graph.h"
#include <cstring>
#include <string>


// can add more Sequential Types
enum class ColoringType { Sequential };

struct StartupOptions {
  std::string inputFile = "";
  ColoringType coloringType = ColoringType::Sequential;
};

StartupOptions parseOptions(int argc, const char **argv) {
  StartupOptions so;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-f") == 0) {
      so.inputFile = argv[i+1];
    } else if (strcmp(argv[i], "-seq") == 0) {
      so.coloringType = ColoringType::Sequential;
    }
  }
  return so;
}

bool checkCorrectness(std::vector<graphNode> &nodes,
                      std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                      std::unordered_map<graphNode, color> &colors) {
  for (auto &node : nodes) {
    if (colors.count(node) == 0)
      return false;

    color curr = colors[node];

    for (auto &nbor : graph[node]) {
      if (colors.count(nbor) == 0) {
        return false;
      }
      if (colors[nbor] == curr) {
        return false;
      }
    }
  }
  return true;
}

int main(int argc, const char **argv) {
  StartupOptions options = parseOptions(argc, argv);

  // TODO: add a read nodes + pairs from file option here
  
  return 0;
}
