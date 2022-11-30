#include "graph.h"
#include "mpi.h"
#include "timing.h"

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

StartupOptions parseOptions(int argc, const char **argv) {
  StartupOptions so;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-f") == 0) {
      so.inputFile = argv[i+1];
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

void receiveAndCopy(MPI_Status &status, const int &chunkSize, std::unordered_map<graphNode, color> &colors) {
    std::vector<int> nodesAndColors;
    nodesAndColors.resize(2 * chunkSize);

    MPI_Recv(&nodesAndColors[0], 2 * chunkSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // any tag is fine
    int recvChunk_i = status.MPI_TAG; // store chunk_i in sender's tag

    for (int i = 0; i < chunkSize; i ++) {
      int node = nodesAndColors[2 * i];
      int color = nodesAndColors[2 * i + 1];
      int prev_color = colors[node];
      // update to be max
      colors[node] = std::max(prev_color, color);
    }
}

void headProc(std::unordered_map<graphNode, std::vector<graphNode>> &graph, std::vector<graphNode> nodes,
              std::unordered_map<graphNode, color> &colors, const int &nproc,
                  const int chunkSize) {
    int numChunks = (float) graph.size() / chunkSize;
    int chunk_i = 0; // index of chunk of work to be done
    int maxTaskPid = 1;

    // initial distribution of chunks to task processors
    for (int pid = 1; pid < nproc; pid ++) {
        if (chunk_i < numChunks) {
            // send the chunk_i and the chunk_size of the particles
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
        receiveAndCopy(status, chunkSize, colors);
        // status.MPI_SOURCE should be a newly freed processor

        MPI_Request chunkReq;
        MPI_Isend(&chunk_i, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &chunkReq); // tag is 0, doesn't matter
        chunk_i += 1;
    }

    // finish processing task responses
    for (int pid = 1; pid < maxTaskPid; pid++) {
        MPI_Status status;
        receiveAndCopy(status, chunkSize, colors);
        MPI_Request doneReq;
        MPI_Isend(0, 0, MPI_INT, pid, chunkSize, MPI_COMM_WORLD, &doneReq); // tag is chunkSize, indicates completed
    }

    // ASSUMING chunk size perfectly divides number of edges!! / pairs.size(), otherwise, probably will have issue
    // now, each node should be colored, and we do a check (and correction as appropriate)

    // can partition nodes (ids) and parallelize checking and correcting, for now have P0 do all though
    int numColors = 0;
    for (size_t i = 0; i < nodes.size(); i++) {
        int node = nodes[i];
        numColors = std::max(numColors, colors[node] + 1);
    }

    for (size_t i = 0 ; i < nodes.size(); i++) {
        int node = nodes[i];
        int color = colors[node];
        for (auto &nbor : graph[node]) {
            if (color == colors[nbor]) {
            colors[node] = numColors++;
            break;
            }
        }
    }
}

void buildGraph(const std::vector<graphNode> &nodes, const std::vector<graphNode> &first, const std::vector<graphNode> &second,
                std::unordered_map<graphNode, std::vector<graphNode>> &graph) {
    for (auto &node : nodes) {
        graph[node] = {};
    }

    for (size_t i = 0; i < first.size(); i++) {
        graph[first[i]].push_back(second[i]);
        graph[second[i]].push_back(first[i]);
    }
}

// builds a partial graph from subset of edges
void buildPartialGraph(const int &chunk_i, const int &chunkSize, 
                const std::vector<graphNode> &first, const std::vector<graphNode> &second,
                std::unordered_map<graphNode, std::vector<graphNode>> &graph) {
    
    /*
    for (int i = chunk_i * chunkSize; i < (chunk_i + 1) * chunkSize; i++) {
    for (auto &node : nodes) {
      graph[node] = {};
    }
    is the initialization necessary?
    */
  
    // only add chunk of edges to graph
    for (int i = chunk_i * chunkSize; i < (chunk_i + 1) * chunkSize; i++) {
      graph[first[i]].push_back(second[i]);
      graph[second[i]].push_back(first[i]);
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
    for (const auto &node : graph) {
        int color = firstAvailableColor(node.first, graph, colors);
        colors[node.first] = color;
    }
}
  
void taskProc(const std::vector<graphNode> &first, const std::vector<graphNode> &second, const int &chunkSize) {
    int chunk_i;
    MPI_Status status;

    while (true) {
        MPI_Recv(&chunk_i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == chunkSize) {
            return;
        }

        // create a new graph from the chunk of edges
        std::unordered_map<graphNode, std::vector<graphNode>> graph;
        buildPartialGraph(chunk_i, chunkSize, first, second, graph);

        // find a coloring for that graph
        std::unordered_map<graphNode, color> colors;
        colorGraph(graph, colors);

        std::vector<int> colorsToSend; // even position is node, next odd number is color of prev even node
        colorsToSend.resize(2 * chunkSize);
        for (const auto &c : colors) {
            colorsToSend.push_back(c.first); // node
            colorsToSend.push_back(c.second); // color of node
        }

        // send the nodes and colors back to Head 
        MPI_Request chunkReq;
        MPI_Isend(&colorsToSend[0], 2 * chunkSize, MPI_INT, 0, chunk_i, MPI_COMM_WORLD, &chunkReq); // tag with chunk_i
    }
}   


int main(int argc, char *argv[]) {
    int pid;
    int nproc;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    StartupOptions options = parseOptions(argc, argv);

    std::vector<graphNode> nodes;
    std::vector<std::pair<graphNode, graphNode>> pairs;
    // first and second elements of pairs
    std::vector<graphNode> first;
    std::vector<graphNode> second;

    if (pid == 0) {
        readGraphFromFile(options.inputFile, nodes, pairs);

        // using first and second as new vectors to not have to create new MPI type,
        // and pass it around
        // but it might be cleaner to just create new type?
        for (size_t i = 0; i < pairs.size(); i ++) {
            first[i] = pairs[i].first;
            second[i] = pairs[i].second;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    Timer t;
    t.reset();

    // broadcast nodes and num of nodes so everyone is aware
    // and pairs and num of pairs
    int numNodes = (int) nodes.size();
    MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    int numPairs = (int) pairs.size();
    MPI_Bcast(&numPairs, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    if (pid != 0) {
        nodes.resize(numNodes);
        first.resize(numPairs);
        second.resize(numPairs);
    }
    MPI_Bcast(&nodes[0], numNodes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&first[0], numPairs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&second[0], numPairs, MPI_INT, 0, MPI_COMM_WORLD);
    // now all processors know the nodes and the edges: ith edge =  (first[i], second[i])!

    if (pid == 0) {
        headProc(graph, nodes, colors, nproc, chunkSize);
    }
    else {
        taskProc(first, second, chunk_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double totalSimulationTime = t.elapsed();

    // return simulation time + check correctness
    if (pid == 0) {
        std::unordered_map<graphNode, std::vector<graphNode>> graph;
        buildGraph(nodes, first, second, graph);

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
        for (auto &color : colors) {
            max = std::max(max, color.second);
        }
        std::cout << max + 1 << " colors\n"; 
        }
    }

    MPI_Finalize();
    return 0;
}

/*
idea:
assume each particle has implicit id
given P = number of processors 
given C = chunk size (or granularity)


for each of N/C chunks do about C^2 work => O(N/C * C^2) = O(N*C) [vs. O(N^2)]
lower C better for amount of work here, but if C = 1 for instance, then will get very unoptimal coloring in correction step

head:
- responsible for giving other particles a chunk to work with
task:
- responsible for finding a greedy coloring for its given chunk

say there are N number of particles
correction phase:
head:
- receive all colorings from each task processor
- check all neighbors BFS style
- if there's something incorrect, use a new color

- return head's final coloring
*/