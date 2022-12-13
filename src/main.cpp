#include "graph.h"
#include "timing.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


// can add more Sequential Types
enum class ColoringType { Sequential, OpenMP, JPOpenMP };

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
    } else if (strcmp(argv[i], "-openmp") == 0) {
      so.coloringType = ColoringType::OpenMP;
    } else if (strcmp(argv[i], "-jpop") == 0) {
      so.coloringType = ColoringType::JPOpenMP;
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
    if (curr == -1) std::cout << "Negative color\n";

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

bool readGraphFromFile(std::string fileName, std::vector<graphNode> &nodes,
                            std::vector<std::pair<graphNode, graphNode>> &pairs) {
  std::ifstream inFile;

  inFile.open(fileName);
  if (!inFile) {
    return false;
  }

  std::string line;

  std::getline(inFile, line);
  std::stringstream sstream(line);
  std::string str;
  std::getline(sstream, str, '\n');
  int numVertices = (int) atoi(str.c_str());

  for (int i = 0; i < numVertices; i++) {
    nodes.push_back(i);
  }

  while(std::getline(inFile, line)) {
    std::stringstream sstream2(line);
    std::getline(sstream2, str, ' ');
    int v1 = (int) atoi(str.c_str());
    std::getline(sstream2, str, '\n');
    int v2 = (int) atoi(str.c_str());

    pairs.push_back(std::make_pair(v1, v2));
  }

  return true;
}

void createCompleteTest(std::vector<graphNode> &nodes,
                        std::vector<std::pair<graphNode, graphNode>> &pairs) {
  int numNodes = 5000;
  nodes.resize(numNodes);
  for (int i = 0; i < numNodes; i++) {
    nodes[i] = i;
  }
  pairs.clear();
  for (int i = 0; i < numNodes; i++) {
    for (int j = i + 1; j < numNodes; j++) {
      pairs.push_back(std::make_pair(i, j));
    }
  }
}

int main(int argc, const char **argv) {
  StartupOptions options = parseOptions(argc, argv);

  // TODO: add a read nodes + pairs from file option here
  std::vector<graphNode> nodes;
  std::vector<std::pair<graphNode, graphNode>> pairs;
  if (!readGraphFromFile(options.inputFile, nodes, pairs)) {
    createCompleteTest(nodes, pairs);
    // std::cerr << "Failed to read graph from input file\n";
  }

  std::unique_ptr<ColorGraph> cg;

  switch (options.coloringType) {
    case ColoringType::Sequential:
      cg = createSeqColorGraph();
      break;
    case ColoringType::OpenMP:
      cg = createOpenMPColorGraph();
      break;
    case ColoringType::JPOpenMP:
      cg = createJPOpenMPColorGraph();

  }

  Timer t;

  std::unordered_map<graphNode, std::vector<graphNode>> graph;
  std::unordered_map<graphNode, color> colors;
  cg->buildGraph(nodes, pairs, graph);
  t.reset();
  cg->colorGraph(graph, colors);

  double time_spent = t.elapsed();
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(5);
  std::cout << "Time spent: " << time_spent << std::endl;
  std::cout << "Colored with ";
  int max = 0;
  for (auto &color : colors) {
    max = std::max(max, color.second);
  }
  std::cout << max + 1 << " colors\n"; 

  if (!checkCorrectness(nodes, graph, colors)) {
    std::cout << "Failed to color graph correctly\n";
    return -1;
  }

  return 0;
}
