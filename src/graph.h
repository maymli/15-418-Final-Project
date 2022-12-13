#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>

typedef int color;
typedef int graphNode;

class ColorGraph {
  public:
    
    
    virtual void buildGraph(std::vector<graphNode> &nodes,
                            std::vector<std::pair<graphNode, graphNode>> &pairs,
                            std::unordered_map<graphNode, std::vector<graphNode>> &graph) = 0;

    virtual void colorGraph(std::unordered_map<graphNode, std::vector<graphNode>> &graph,
                            std::unordered_map<graphNode, color> &colors) = 0;
};

std::unique_ptr<ColorGraph> createSeqColorGraph();
std::unique_ptr<ColorGraph> createOpenMPColorGraph();
std::unique_ptr<ColorGraph> createJPOpenMPColorGraph();
