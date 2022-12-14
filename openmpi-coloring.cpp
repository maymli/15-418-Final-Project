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

void colorGraph(int start, int end, std::vector<graphNode> &nodes, std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                     std::unordered_map<graphNode, color> &colors) {

    // std::cerr << "coloring graph!" << "\n";
    for (int i = start; i < end; i ++) {
      auto node = nodes[i];
      int color = firstAvailableColor(node, graph, colors);
      colors[node] = color;
      // std::cerr << "node " << node << ", color " << colors[node] << "\n";
    }
}

//
void receiveAndCopy(std::vector<graphNode> &nodes, MPI_Status &status, const int &chunkSize, std::unordered_map<graphNode, color> &colors) {
  std::vector<int> nodeColors;
  nodeColors.resize(chunkSize);

  MPI_Recv(&nodeColors[0], chunkSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // any tag is fine
  
  int recvChunk_i = status.MPI_TAG; // store chunk_i in sender's tag
  int start = chunkSize * recvChunk_i;
  int end = chunkSize * (recvChunk_i + 1);

  for (int i = start; i < end; i ++) {
    int node = nodes[i];
    int color = nodeColors[i - start];
   
    // update to be max
    int prev_color = colors[node];
    colors[node] = std::max(prev_color, color);
  }
}

void headProc(std::vector<graphNode> &nodes, std::unordered_map<graphNode, color> &colors, const int &nproc,
                  const int chunkSize) {
  int numChunks = (float) nodes.size() / chunkSize;
  int chunk_i = 0; // index of chunk of work to be done
  int maxTaskPid = 1;

  // initial distribution of chunks to task processors
  for (int pid = 1; pid < nproc; pid ++) {
      if (chunk_i < numChunks) {
          // send the chunk_i
          MPI_Request chunkReq;
          MPI_Isend(&chunk_i, 1, MPI_INT, pid, 0, MPI_COMM_WORLD, &chunkReq); // tag is 0, doesn't matter
          chunk_i += 1;
          maxTaskPid += 1;
      }
      else
          break;
  }

  // assign remainder of chunks when task processors free up
  while (chunk_i < numChunks) {
      MPI_Status status;
      receiveAndCopy(nodes, status, chunkSize, colors); // status.MPI_SOURCE should be a newly freed processor

      MPI_Request chunkReq;
      MPI_Isend(&chunk_i, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &chunkReq); // tag is 0, doesn't matter
      chunk_i += 1;
  }
  
  // finish processing task responses
  for (int pid = 1; pid < maxTaskPid; pid++) {
      MPI_Status status;
      receiveAndCopy(nodes, status, chunkSize, colors);

      MPI_Request doneReq;
      MPI_Isend(0, 0, MPI_INT, pid, chunkSize, MPI_COMM_WORLD, &doneReq); // tag is chunkSize, indicates completed
  }
}

void headCorrection(std::unordered_map<graphNode, std::vector<graphNode>> graph, const int numPairs, std::vector<graphNode> nodes, std::unordered_map<graphNode, color> &colors, const int &nproc,
                  const int chunkSize) {

  // can maybe partition nodes (ids) and parallelize checking and correcting? for now have P0 do all though
  int numColors = 0;
  int maxColor = 0;
  for (size_t i = 0; i < nodes.size(); i++) {
      int node = nodes[i];
      maxColor = std::max(maxColor, colors[node]);
      // std::cerr << "node " << node << " color " << colors[node] << "\n";
  }
  numColors = maxColor + 1;
  std::cerr << "num colors " << numColors << "\n";

  for (size_t i = 0 ; i < nodes.size(); i++) {
    int node = nodes[i];
    int color = colors[node];

    for (auto &nbor : graph[node]) {
        if (color == colors[nbor]) {
          std::cerr << "(node " << node << ", " << nbor << "): color " << color << " and new " << numColors << "\n";
          colors[node] = numColors++;
          //std::cerr << "num colors " << numColors << "\n";
          break;
        }
    }
  }
}
  
void taskProc(std::unordered_map<graphNode, std::vector<graphNode>> &graph, std::vector<graphNode> &nodes, const int &chunkSize) {
    int chunk_i;
    MPI_Status status;

    while (true) {
      MPI_Recv(&chunk_i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (status.MPI_TAG == chunkSize) {
          return;
      }

      // find a coloring for these specific nodes
      std::unordered_map<graphNode, color> colors;
      // if (chunk_i == 0)
      int start = chunk_i * chunkSize;
      int end = (chunk_i + 1) * chunkSize;
      colorGraph(start, end, nodes, graph, colors);

      std::vector<int> colorsToSend; // even position is node, next odd number is color of prev even node
      for (int i = start; i < end; i++) {
        auto node = nodes[i];
        colorsToSend.push_back(colors[node]); // color of node
        
        /*
        if (chunk_i == 1)
          std::cerr << "(node " << node << ", color " << colorsToSend[i] << "), ";
        */
      }

      // send the nodes and colors back to Head 
      MPI_Request chunkReq;
      MPI_Isend(&colorsToSend[0], chunkSize, MPI_INT, 0, chunk_i, MPI_COMM_WORLD, &chunkReq); // tag with chunk_i
    }
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
  
  std::unordered_map<graphNode, color> colors;

  // Start timer
  MPI_Barrier(MPI_COMM_WORLD);
  Timer t;
  t.reset();

  MPI_Datatype MPI_PAIR;
  MPI_Type_contiguous(2, MPI_INT, &MPI_PAIR);
  MPI_Type_commit(&MPI_PAIR);

  // Broadcast nodes and pairs from head to task processes
  MPI_Bcast(&nodes[0], numNodes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pairs[0], numPairs, MPI_PAIR, 0, MPI_COMM_WORLD);

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

  // Start head and task processes
  // Head allocates a chunk of nodes that each task process colors
  int chunkSize = 2500;
  if (pid == 0) {
      headProc(nodes, colors, nproc, chunkSize);
  }
  else {
      taskProc(graph, nodes, chunkSize);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  double totalSimulationTime = t.elapsed();

  // Return simulation time + check correctness
  if (pid == 0) {

    // Correct any mis-colored nodes
    headCorrection(graph, numPairs, nodes, colors, nproc, chunkSize);

    // Finish timing and print information
    printf("total simulation time: %.6fs\n", totalSimulationTime);
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
}
