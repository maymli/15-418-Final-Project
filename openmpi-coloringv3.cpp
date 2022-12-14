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
  int s;
};
StartupOptions parseOptions(int argc, char **argv) {
  StartupOptions so;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-f") == 0) {
      so.inputFile = argv[i+1];
    }
    else if (strcmp (argv[i], "-s") == 0) {
      so.s = std::stoi(argv[i+1]);
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
                      std::unordered_map<graphNode, color> &colors, int pid) {                 
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
                  std::unordered_map<graphNode, std::vector<graphNode>> &graph, int pid) {
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

  // Start timer
  MPI_Barrier(MPI_COMM_WORLD);
  Timer t;
  t.reset();

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
  buildGraph(nodes, pairs, graph, pid);
  #pragma endregion
  
  // Each processor has nodes + graph now
  // TODO: assume num nodes divisible by num procs for now, fix remainder later
  
  #pragma region partition
  // using "natural" ordering (has implicit mapping of nodes to processors)
  // also assuming node = i in nodes
  int nodeChunkSize = numNodes / nproc;
  int start = nodeChunkSize * pid;
  int end = nodeChunkSize * (pid + 1);
  
  std::vector<graphNode> nodesToColor;
  for (int i = start; i < end; i ++) {
    nodesToColor.push_back(nodes[i]);
  }
  
  // create mapping from boundary nodes to processors that will need it
  std::unordered_map<graphNode, std::vector<graphNode>> boundaryMap;

  int proc;
  for (auto node : nodesToColor) {
    for (auto nbor : graph[node]) {
      proc = nbor / nodeChunkSize;
      if (proc != pid) {
        // then node is a boundary node, whose nbor is in proc
        if (boundaryMap.find(node) == boundaryMap.end()) {
          boundaryMap[node] = {};
          boundaryMap[node].resize(nproc);
        }
        boundaryMap[node][proc] = 1;
      }
    }
  }
  
  #pragma endregion

  int s = options.s; // let the user control this?
  int num_supersteps = numNodes / s;

  std::vector<int> unfinishedProcessors;
  unfinishedProcessors.resize(nproc);
  int unfinished = 1;
  
  std::unordered_set<int> otherPids;
  for (int proc = 0; proc < nproc; proc ++) {
    if (proc != pid) {
      otherPids.insert(proc);
    }
  }

  while (unfinished != 0) {
    if (nodesToColor.size() > 0) {
    std::vector<graphNode> boundaryNodes;

    // for each superstep
    for (int step = 0; step < num_supersteps; step++) 
    {
      std::vector<std::vector<color>> colorsToSend(nproc);
      
      int limit = std::min((step + 1) * s, (int) nodesToColor.size());
      for (int i = step * s; i < limit; i++) {
        int node = nodesToColor[i];

        // assign color to node
        int color = firstAvailableColor(node, graph, colors);
        colors[node] = color;

        if (boundaryMap.find(node) != boundaryMap.end()) {
          boundaryNodes.push_back(node);

          for (int proc = 0; proc < nproc; proc ++) {
            if (boundaryMap[node][proc]) {
              colorsToSend[proc].push_back(node);
              colorsToSend[proc].push_back(color);
            }
          }
        }
      }

      #pragma region communication
      // send over number of boundary nodes that will be sent
      for (auto proc : otherPids) {
        int numNew = colorsToSend[proc].size();
        MPI_Request numColors;
        MPI_Isend(&numNew, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &numColors); // tag is 0, doesn't matter
      }
      // receive number of boundary nodes that will be received
      std::vector<int> numsNew;
      numsNew.resize(nproc);
      for (auto proc : otherPids) {
        MPI_Status status;
        MPI_Recv(&numsNew[proc], 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &status); // tag is 0, doesn't matter
      }
      // TODO: try doing IRecv and Waitall?
      
      // actually send over the boundary nodes
      for (auto proc : otherPids) {
        MPI_Request newColors;
        MPI_Isend(&colorsToSend[proc][0], colorsToSend[proc].size(), MPI_INT, proc, 0, MPI_COMM_WORLD, &newColors); // tag is 0, doesn't matter
      }

      std::vector<std::vector<int>> newColors;
      newColors.resize(nproc);
      for (auto proc : otherPids) {
        newColors[proc].resize(numsNew[proc]);
      }
      // receive the boundary nodes
      for (auto proc : otherPids) {
        MPI_Status status;
        MPI_Recv(&newColors[proc][0], numsNew[proc], MPI_INT, proc, 0, MPI_COMM_WORLD, &status); // tag is 0, doesn't matter
      }
      #pragma endregion

      // update boundary node neighbors w/ new colors
      for (auto proc : otherPids) {
        int numColors = newColors[proc].size();
        for (int i = 0; i < numColors; i += 2) {
          int node = newColors[proc][i];
          int color = newColors[proc][i + 1];
          colors[node] = color;
        }
      }
    }

    // at this point all colors must have been received
    std::unordered_set<int> added;
    std::vector<graphNode> nodesToRecolor;

    //std::cerr << "\nP" << pid << ": boundary nodes = " << boundaryNodes.size() << "\n";
    for (auto node : boundaryNodes) {
      for (auto nbor : graph[node]) {
        if (colors[node] == colors[nbor]) {
          if (node + 23 % 5 < nbor + 23 % 5) { // TODO: use a different pseudorandom function
            if (added.find(node) == added.end()) {
              added.insert(node);
              nodesToRecolor.push_back(node);
            }
          }
        }
      }
    }
    nodesToColor.swap(nodesToRecolor);
    }

    #pragma region check
    // Check to see if all processors are finished!
    unfinished = nodesToColor.size(); // 0 if finished, > 0 if unfinished!

    //MPI_Allgather(&unfinished, 1, MPI_INT, &unfinishedProcessors[0], 1, MPI_INT, MPI_COMM_WORLD);
    // send over number of boundary nodes that will be sent
    for (auto proc : otherPids) {
      MPI_Request finished;
      MPI_Isend(&unfinished, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &finished); // tag is 0, doesn't matter
    }
    // receive number of boundary nodes that will be received
    for (auto proc : otherPids) {
      MPI_Status status;
      MPI_Recv(&unfinishedProcessors[proc], 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &status); // tag is 0, doesn't matter
    }
    
    if (!unfinished) { // finished, aka nodesToColor.size() == 0
      otherPids.clear();
    }
    else { // unfinished
      for (int proc = 0; proc < nproc; proc ++) {
        if (proc != pid && !unfinishedProcessors[proc]) { // if another proc finished, remove it
          otherPids.erase(proc);
        }
      }
    }
    #pragma endregion
  }

  // gather all node colors (insignificant time cost)
  std::vector<color> finalColors;
  finalColors.resize(nodeChunkSize);
  for (int i = start; i < end; i ++) {
    finalColors[i - start] = colors[nodes[i]];
  }

  std::vector<color> rootColors;
  rootColors.resize(numNodes);
  MPI_Gather(&finalColors[0], nodeChunkSize, MPI_INT, &rootColors[0], nodeChunkSize, MPI_INT, 0, MPI_COMM_WORLD);

  if (pid == 0) {
    for (int i = 0; i < numNodes; i ++) {
      colors[i] = rootColors[i];
    }
  }

  #pragma region finish
  // Return simulation time + correction
  if (pid == 0) {
    double time_spent = t.elapsed();
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "Time spent: " << time_spent << std::endl;
    if (!checkCorrectness(nodes, graph, colors, pid)) {
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
