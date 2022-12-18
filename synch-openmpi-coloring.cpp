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
    if (colors.count(node) == 0 || colors[node] == -1) {
      std::cout << "missing color for " << node << std::endl;
      return false;
    }

    color curr = colors[node];

    for (auto &nbor : graph[node]) {
      if (colors.count(nbor) == 0) {
        std::cout << "missing color for " << nbor << std::endl;
        return false;
      }
      if (colors[nbor] == curr) {
        std::cout << "same color " << curr << " for " << node << " and " << nbor << std::endl;
        return false;
      }
    }
  }
  return true;
}

void buildGraph(std::vector<graphNode> &nodes, std::vector<std::pair<int, int>> &pairs,
                  std::unordered_map<graphNode, std::vector<graphNode>> &graph) {
  std::unordered_map<graphNode, std::unordered_set<graphNode>> edges;
  for (auto &node : nodes) {
    graph[node] = {};
    edges[node] = {};
  }

  for (auto &pair : pairs) {
    if (edges[pair.first].count(pair.second) == 0) {
      graph[pair.first].push_back(pair.second);
      graph[pair.second].push_back(pair.first);

      edges[pair.first].insert(pair.second);
      edges[pair.second].insert(pair.first);
    }
  }
}

void reorder(int startNode, int endNode, std::unordered_map<graphNode, std::vector<graphNode>> &graph,
      std::vector<graphNode> &order) {
  std::unordered_map<graphNode, graphNode> sums;
  for (int i = startNode; i < endNode; i++) {
    sums[i] = 0;
    for (const auto &nbor : graph[i]) {
      if (nbor >= endNode) {
        sums[i] += 1;
      }
    }
  }
  std::vector<std::pair<graphNode, int>> pairs(sums.begin(), sums.end());
  std::sort(pairs.begin(), pairs.end(), 
                  [] (const std::pair<graphNode, int> &lhs, const std::pair<graphNode, int> &rhs) {
                    return lhs.second < rhs.second;
                  });
  for (const auto &pair : pairs) {
    order.emplace_back(pair.first);
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

int nodeToProc(graphNode node, int totalNodes, int nproc) {
  int pid = (node / (((float) totalNodes) / nproc));
  int endNode = (pid + 1) * (((float) totalNodes) / nproc);
  if (pid != nproc - 1 && node >= endNode) pid++;

  return pid;
}

/* void colorNodes(std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                std::unordered_map<graphNode, color> &colors,
                int pid, int nproc, int totalNodes) {
  
  int startNode = pid * (((float) totalNodes) / nproc);
  int endNode = (pid + 1) * (((float) totalNodes) / nproc);
  std::vector<MPI_Request> requests;
  // requests.resize(totalNodes);
  int numReqs = 0;
  int numWait = 0;

  // map pids to sets they need
  
  std::vector<graphNode> order;
  reorder(startNode, endNode, graph, order);

 for (int node = startNode; node < endNode; node++) {
 // for (auto &node : order) {
    std::vector<graphNode> smallerNbors;
    std::vector<graphNode> largerNbors;
    for (const auto &nbor : graph[node]) {
      if (nbor < startNode || nbor >= endNode) {
        if (colors[nbor] == -1 && node < nbor) {
          largerNbors.emplace_back(nbor);
        } else if (node > nbor) {
          smallerNbors.emplace_back(nbor);
        }
      }
    }

    for (const auto &nbor : largerNbors) {
      // MPI_Status status;
      int color2 = 0;
      MPI_Recv(&color2, 1, MPI_INT, nodeToProc(nbor, totalNodes, nproc), (int) nbor, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      colors[nbor] = color2;
      numWait++;
    }

    // wait for colors to come in, then assign colors 
    int color = firstAvailableColor(node, graph, colors);
    colors[node] = color;
    
    std::vector<bool> sent_pids;
    sent_pids.resize(nproc, false);

    // send color to nbor that needs it and is smaller
    for (const auto &nbor : smallerNbors) {
      int send_pid = nodeToProc(nbor, totalNodes, nproc);
      if (!sent_pids[send_pid]) {
	requests.emplace_back();
        MPI_Isend(&colors[node], 1, MPI_INT, send_pid, node, MPI_COMM_WORLD, &requests[numReqs++]);
        sent_pids[send_pid] = true;
      }
    }
  }
  MPI_Waitall(numReqs, &requests[0], MPI_STATUS_IGNORE);
} */

void colorNodes(std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                std::unordered_map<graphNode, color> &colors,
                int pid, int nproc, int totalNodes) {
  
  int startNode = pid * (((float) totalNodes) / nproc);
  int endNode = (pid + 1) * (((float) totalNodes) / nproc);
  std::vector<MPI_Request> requests;
  int numReqs = 0;
  int numRecv = 0;

  std::unordered_set<graphNode> fnbors;
  std::vector<MPI_Request> recv_reqs;
  recv_reqs.resize(totalNodes);

  std::unordered_map<graphNode, std::vector<graphNode>> uncolored;
  std::unordered_map<graphNode, std::vector<graphNode>> smaller;

  for (int i = startNode; i < endNode; i++) {
    uncolored[i] = {};
    smaller[i] = {};
    for (const auto &nbor : graph[i]) {
      if (nbor >= endNode) {
        fnbors.insert(nbor);
        uncolored[i].emplace_back(nbor);
      } else if (nbor < startNode) {
        smaller[i].emplace_back(nbor);
      }
    }
    if (uncolored[i].empty()) {
      uncolored.erase(i);
    }
  }

  for (const auto &fnbor : fnbors) {
    MPI_Irecv(&colors[fnbor], 1, MPI_INT, nodeToProc(fnbor, totalNodes, nproc), 
                    fnbor, MPI_COMM_WORLD, &recv_reqs[numRecv++]);
  }

  int needColor = endNode - startNode;

  for (int node = startNode; node < endNode; node++) {
    if (uncolored.find(node) != uncolored.end()) continue; 
    // wait for colors to come in, then assign colors 
    int color = firstAvailableColor(node, graph, colors);
    colors[node] = color;
    
    std::vector<bool> sent_pids;
    sent_pids.resize(nproc, false);

    // send color to nbor that needs it and is smaller
    for (const auto &nbor : smaller[node]) {
      int send_pid = nodeToProc(nbor, totalNodes, nproc);
      if (!sent_pids[send_pid]) {
	requests.emplace_back();
        MPI_Isend(&colors[node], 1, MPI_INT, send_pid, node, MPI_COMM_WORLD, &requests[numReqs++]);
        sent_pids[send_pid] = true;
      }
    }
    needColor--;
  }

  while (needColor > 0) {
    std::vector<graphNode> toRemove;
    for (auto &node : uncolored) {
      auto fnbor = node.second.begin();
      while (fnbor != node.second.end()) {
        if (colors[*fnbor] != -1) {
          node.second.erase(fnbor);
        } else {
          fnbor++;
        }
      }
      if (node.second.empty()) {
        int color = firstAvailableColor(node.first, graph, colors);
        colors[node.first] = color;
    
        std::vector<bool> sent_pids;
        sent_pids.resize(nproc, false);

        for (const auto &nbor : smaller[node.first]) {
          int send_pid = nodeToProc(nbor, totalNodes, nproc);
          if (!sent_pids[send_pid]) {
	    requests.emplace_back();
            MPI_Isend(&colors[node.first], 1, MPI_INT, send_pid, node.first, MPI_COMM_WORLD, &requests[numReqs++]);
            sent_pids[send_pid] = true;
          }
        }
        needColor--;
        toRemove.emplace_back(node.first);
      }
    }
    for (auto &node : toRemove) {
      uncolored.erase(node);
    }
    int flag;
    MPI_Testall(numRecv, &recv_reqs[0], &flag, MPI_STATUS_IGNORE);
  }
  MPI_Waitall(numRecv, &recv_reqs[0], MPI_STATUS_IGNORE);
  MPI_Waitall(numReqs, &requests[0], MPI_STATUS_IGNORE);
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
  
  MPI_Datatype MPI_PAIR;
  MPI_Type_contiguous(2, MPI_INT, &MPI_PAIR);
  MPI_Type_commit(&MPI_PAIR);

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

  // Broadcast nodes and pairs from head to task processes
  MPI_Bcast(&nodes[0], numNodes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pairs[0], numPairs, MPI_PAIR, 0, MPI_COMM_WORLD);
  // Now all processors know the nodes and the pairs/edges

  // Each processor builds the full graph
  std::unordered_map<graphNode, std::vector<graphNode>> graph;
  buildGraph(nodes, pairs, graph);

  
  std::unordered_map<graphNode, color> colors;
  for (const auto &node : nodes) {
    colors[node] = - 1;
  }
  
  // Start timer
  MPI_Barrier(MPI_COMM_WORLD);
  Timer t;
  t.reset();


  // Start head and task processes
  // Head allocates a chunk of nodes that each task process colors
  colorNodes(graph, colors, pid, nproc, numNodes);
  
  int startNode = pid * (((float) numNodes) / nproc);
  int endNode = (pid + 1) * (((float) numNodes) / nproc);
  int colorArray[(endNode - startNode) * 2];

  for (int i = startNode; i < endNode; i++) {
    int index = 2 * (i - startNode);
    colorArray[index] = i;
    colorArray[index + 1] = colors[i];
  }

  int chunkSizes[nproc];
  int startPos[nproc];
  startPos[0] = 0;
  chunkSizes[0] = ((float) numNodes) / nproc;
  chunkSizes[0] *= 2;

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
