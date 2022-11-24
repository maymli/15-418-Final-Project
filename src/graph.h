#include <vector>
#include <unordered_map>
#include <utility>

typedef int color;
typedef int graphNode;

class ColorGraph {
  public:
    std::unordered_map<graphNode, std::vector<graphNode>> graph;
    
    virtual void buildGraph(std::vector<graphNode> &nodes,
                            std::vector<std::pair<graphNode, graphNode>> &pairs);

    virtual void colorGraph(std::unordered_map<graphNode, color> &colors);
};


