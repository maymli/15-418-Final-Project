#include "src/graph.h"
#include "mpi.h"
#include "src/timing.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>

struct StartupOptions {
  std::string inputFile = "";
  int max_iters;
};

StartupOptions parseOptions(int argc, char **argv) {
  StartupOptions so;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-f") == 0) {
      so.inputFile = argv[i+1];
    }
    else if (strcmp (argv[i], "-iters") == 0) {
      so.max_iters = std::stoi(argv[i+1]);
    }
  }
  return so;
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

  inFile.close();
  return true;
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

std::vector<graphNode> findIncorrectNodes(std::vector<graphNode> &nodes,
                      std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                      std::unordered_map<graphNode, color> &colors) {
  std::vector<graphNode> nodesToRecolor;
  for (auto &node : nodes) {
    color curr = colors[node];
    for (auto &nbor : graph[node]) {
      if (colors[nbor] == curr) {
        if (node + 23 % 5 < nbor + 23 % 5) {
          // so that only one node has to re-color
          nodesToRecolor.push_back(node);
          break;
        }
      }
    }
  }
  return nodesToRecolor;
}

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

void colorGraph(std::vector<graphNode> &nodes, std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                     std::unordered_map<graphNode, color> &colors, std::vector<int> &newColors, int diff) {

    // std::cerr << "coloring graph! " << diff << "\n";
    for (auto node : nodes) {
      int color = firstAvailableColor(node, graph, colors);
      // TODO: this is assuming the node == i of node in nodes
      // it's not as generalizable, but it should be easy enough to re-map node values to indices
      colors[node] = color;
      newColors[node - diff] = color;
    }
    // std::cerr << "colored graph! " << diff << "\n";
}

void headCorrection(std::unordered_map<graphNode, std::vector<graphNode>> graph, std::vector<graphNode> nodes, std::unordered_map<graphNode, color> &colors) {
  //int wrong = 0;
  for (size_t i = 0 ; i < nodes.size(); i++) {
    int node = nodes[i];
    int color = colors[node];

    for (auto &nbor : graph[node]) {
        if (color == colors[nbor]) {
          //wrong += 1;
          colors[node] = firstAvailableColor(node, graph, colors);//numColors++;
          break;
        }
    }
  }
  //std::cerr << "wrongs: " << wrong << std::endl;
}

int main(int argc, char *argv[]) {
  
  #pragma region init
  // MPI initialization:
  int pid;
  int nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Variable initialization
  StartupOptions options = parseOptions(argc, argv);
  std::vector<graphNode> nodes;
  std::vector<std::pair<graphNode, graphNode>> pairs;
  std::unordered_map<graphNode, color> colors;

  // Head reads in file values
  if (pid == 0) {
    readGraphFromFile(options.inputFile, nodes, pairs);
  }
  
  // Broadcast number of nodes and pairs of nodes
  int numNodes = (int) nodes.size();
  MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  int numPairs = (int) pairs.size();
  MPI_Bcast(&numPairs, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  
  // Resize task processes' variables
  if (pid != 0) {
      nodes.resize(numNodes);
      pairs.resize(numPairs);
  }

  MPI_Datatype MPI_PAIR;
  MPI_Type_contiguous(2, MPI_INT, &MPI_PAIR);
  MPI_Type_commit(&MPI_PAIR);

  // Broadcast nodes and pairs from head to task processes
  MPI_Bcast(&nodes[0], numNodes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pairs[0], numPairs, MPI_PAIR, 0, MPI_COMM_WORLD);

  // Now all processors know the nodes and the pairs/edges

  // Each processor builds the full graph
  std::unordered_map<graphNode, std::vector<graphNode>> graph;
  buildGraph(nodes, pairs, graph);
  #pragma endregion

  MPI_Barrier(MPI_COMM_WORLD);
  Timer t;
  t.reset();

  std::vector<color> colors_v;
  colors_v.resize(numNodes);

  // TODO: assume num nodes divisible by num procs for now
  int nodeChunkSize = numNodes / nproc;
  int start = nodeChunkSize * pid;
  int end = nodeChunkSize * (pid + 1);
  
  std::vector<graphNode> nodesToColor;
  for (int i = start; i < end; i ++) {
    nodesToColor.push_back(nodes[i]);
  }

  std::vector<color> newColors;
  newColors.resize(nodeChunkSize);

  int max_iters = options.max_iters;
  int num_iters = 0;
  while (num_iters < max_iters) {

    // each processor colors its respective chunk of nodes w/ knowledge of full graph
    colorGraph(nodesToColor, graph, colors, newColors, start);

    // update all other processors
    // experiment w/ using AllGather vs Gather + Scatter
    MPI_Allgather(&newColors[0], nodeChunkSize, MPI_INT, &colors_v[0], nodeChunkSize, MPI_INT, MPI_COMM_WORLD);
    
    // update all nodes in colors map with updated colors vector
    for (int i = 0; i < numNodes; i ++) {
      colors[nodes[i]] = colors_v[i];
    }

    // check if its assignment is valid 
    std::vector<graphNode> nodesToRecolor = findIncorrectNodes(nodesToColor, graph, colors);
    nodesToColor.swap(nodesToRecolor);
    num_iters += 1;
  }

  #pragma region finish
  // Return simulation time + correction
  if (pid == 0) {

    // Correct any mis-colored nodes
    headCorrection(graph, nodes, colors);

    // Finish timing and print information
    double time_spent = t.elapsed();
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "Time spent: " << time_spent << std::endl;
    if (!checkCorrectness(nodes, graph, colors)) {
      std::cout << "Failed to color graph correctly\n";
      return -1;
    } else {
      std::cout << "Colored with ";
      int max = 0;
      for (auto &node : nodes) {
        max = std::max(max, colors[node]);
      }
      std::cout << max + 1 << " colors\n"; 
    }
  }

  MPI_Finalize();
  return 0;
  #pragma endregion
}
