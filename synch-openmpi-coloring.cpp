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
};

StartupOptions parseOptions(int argc, char **argv) {
  StartupOptions so;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-f") == 0) {
      so.inputFile = argv[i+1];
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

int nodeToProc(graphNode node, int totalNodes, int nproc) {
  return (int) (node / (((float) totalNodes) / nproc));
}

std::unordered_map<graphNode, color> colorNodes(std::unordered_map<graphNode, std::vector<graphNode>> &graph, 
                int pid, int nproc, int totalNodes) {
  
  std::unordered_map<graphNode, color> colors;
  int startNode = pid * (((float) totalNodes) / nproc);
  int endNode = (pid + 1) * (((float) totalNodes) / nproc);
  std::cout << startNode << ", " << endNode << std::endl;

  for (int node = startNode; node < endNode; node++) {
    std::vector<graphNode> smallerNbors;
    std::vector<graphNode> largerNbors;
    for (const auto &nbor : graph[node]) {
      if (nbor < startNode || nbor > endNode) {
        if (node < nbor) {
          largerNbors.emplace_back(nbor);
        } else {
          smallerNbors.emplace_back(nbor);
        }
      }
    }

    for (const auto &nbor : largerNbors) {
      MPI_Status status;
      int color = 0;
      MPI_Recv(&color, 1, MPI_INT, nodeToProc(nbor, totalNodes, nproc), nbor, MPI_COMM_WORLD, &status);
      colors[nbor] = color;
    }

    // wait for colors to come in, then assign colors 
    int color = firstAvailableColor(node, graph, colors);
    colors[node] = color;

    // send color to nbor that needs it and is smaller
    for (const auto &nbor : smallerNbors) {
      MPI_Request send_req;
      MPI_Isend(&color, 1, MPI_INT, nodeToProc(nbor, totalNodes, nproc), node, MPI_COMM_WORLD, &send_req);
    }
  }
  return colors;
}
  
int main(int argc, char *argv[]) {
  
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

  // Head reads in file values
  if (pid == 0) {
    readGraphFromFile(options.inputFile, nodes, pairs);
  }
  
  // Start timer
  MPI_Barrier(MPI_COMM_WORLD);
  Timer t;
  t.reset();

  // Broadcast number of nodes and pairs of nodes
  int numNodes = (int) nodes.size();
  std::cout << nodes.size() << std::endl;
  MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  int numPairs = (int) pairs.size();
  MPI_Bcast(&numPairs, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  // Resize task processes' variables
  if (pid != 0) {
      nodes.resize(numNodes);
      pairs.resize(numPairs);
  }
  // std::cerr << "hello i am " << pid << "! i have " << numNodes << " nodes and " << numPairs << " pairs\n";

  // Broadcast nodes and pairs from head to task processes
  MPI_Bcast(&nodes[0], numNodes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pairs[0], numPairs, MPI_INT, 0, MPI_COMM_WORLD);
  // Now all processors know the nodes and the pairs/edges

  // Each processor builds the full graph
  std::unordered_map<graphNode, std::vector<graphNode>> graph;
  buildGraph(nodes, pairs, graph);

  std::cout << nproc << ", " << numNodes << std::endl;

  // Start head and task processes
  // Head allocates a chunk of nodes that each task process colors
  auto colors = colorNodes(graph, pid, nproc, numNodes);
  
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << pid << ": " << colors.size() << std::endl;
  
  int startNode = pid * (((float) numNodes) / nproc);
  int endNode = (pid + 1) * (((float) numNodes) / nproc);
  int colorArray[(endNode - startNode) * 2];

  for (int i = startNode; i < endNode; i++) {
    int index = 2 * (i - startNode);
    colorArray[index] = i;
    colorArray[index++] = colors[i];
  }

  int chunkSizes[nproc];
  int startPos[nproc];
  startPos[0] = 0;
  chunkSizes[0] = ((float) numNodes) / nproc;
  for (int i = 1; i < nproc; i++) {
    int start = i * (((float) numNodes) / nproc);
    int end = (i + 1) * (((float) numNodes) / nproc);
    startPos[i] = startPos[i - 1] + chunkSizes[i - 1];
    chunkSizes[i] = (end - start) * 2;
  }

  // does this work
  int allColors[2 * numNodes];

  MPI_Gatherv(colorArray, chunkSizes[pid], MPI_INT, allColors, chunkSizes, startPos, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  
  // Return simulation time + check correctness
  if (pid == 0) {
    // Finish timing and print information
    double time_spent = t.elapsed();
    colors.clear();
    for (int i = 0; i < numNodes; i++) {
      colors[allColors[2 * i]] = allColors[2 * i + 1];
    }
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
}
