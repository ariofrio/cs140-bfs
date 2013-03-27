#include <cilk.h>
#include <cilkview.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <deque>
#include <list>
#include <set>
using namespace std;

typedef list< pair<int, int> > list_pair;

class Graph {
public:
  Graph(string filename);
  void print();
  int degree(int node) {
    if(node == indexes.size() - 1)
      return neighbors.size() - indexes[node];
    else
      return indexes[node+1] - indexes[node];
  }
  int neighbor(int node, int index) {
    return neighbors[indexes[node] + index];
  }
  int size() {
    return indexes.size();
  }


private:
  vector<int> neighbors;
  vector<int> indexes;
  vector< pair<int, int> > edges; // might be empty if cache was used

  void read_cache(istream &file);
  void write_cache(ofstream &file);
  void read_edges(istream &file);
  void build_graph();
};

Graph::Graph(string filename) {
  if(filename == "-") {
    read_edges(cin);
    build_graph();
  } else {
    ifstream cache((filename + ".cache").c_str());
    if(cache) {
      read_cache(cache);
    } else {
      cache.close();

      ifstream file(filename.c_str());
      read_edges(file);
      build_graph();

      ofstream cache_out((filename + ".cache").c_str());
      write_cache(cache_out);
    }
  }
}

void Graph::read_cache(istream &file) {
  int neighbor_count;
  file >> neighbor_count;
  neighbors.resize(neighbor_count);
  for(int i=0; i<neighbor_count; i++)
    file >> neighbors[i];

  int index_count;
  file >> index_count;
  indexes.resize(index_count);
  for(int i=0; i<index_count; i++)
    file >> indexes[i];
}

void Graph::write_cache(ofstream &file) {
  file << neighbors.size() << endl;
  for(int i=0; i<neighbors.size(); i++)
    file << neighbors[i] << endl;

  file << indexes.size() << endl;
  for(int i=0; i<indexes.size(); i++)
    file << indexes[i] << endl;
}

void Graph::read_edges(istream &file) {
  int _, node_count, edge_count;
  file >> _ >> node_count >> edge_count;
  if(_ != 0) {
    cerr << "Invalid input format: first line should be '0 NODES EDGES'" << endl;
    exit(1);
  }
  indexes.resize(node_count);
  edges.reserve(edge_count);

  while(!file.eof()) {
    pair<int, int> edge;
    file >> edge.first;
    file >> edge.second;
    if(file.fail()) break;
    edges.push_back(edge);
  }
}

void Graph::build_graph() {
  for(int i=0; i<edges.size(); i++)
    indexes.at(edges[i].first)++;
  for(int i=1; i<indexes.size(); i++)
    indexes[i] += indexes[i-1];

  neighbors.resize(indexes.back());
  
  for(int i=0; i<edges.size(); i++) {
    int index = --indexes.at(edges[i].first);
    neighbors.at(index) = edges[i].second;
  }
}

void Graph::print() {
  printf("\nGraph has %d vertices and %d edges\n",
      indexes.size(), neighbors.size());
}

class BFS {
public:
  Graph& graph;
  int root;
  BFS(Graph& graph, int root) : graph(graph), root(root) {
    node_level.resize(graph.size(), -1);
    node_parent.resize(graph.size(), -1);
    node_queued.resize(graph.size(), false);
  };

  vector<int> node_level;
  vector<int> node_parent;
  deque<bool> node_queued; // vector<bool> is specialized and causes data races
  void run();
  void print_statistics();
  void print_results();
  bool validate();

private:
  bool validation_failed;
  void validate_bfs_node_points_root(int node);
  void validate_bfs_node_points_root(int node, list<int> &offspring);
  void validate_bfs_edge_level(int child, int parent);
  void validate_graph_edge_level(int child, int parent);
  void validate_graph_edges_in_bfs(int parent);
  void validate_bfs_edge_in_graph(int node, int parent);
#if method == 0
  list_pair queue;
  list_pair next;
#elif method == 1
  void process_queue(list_pair &queue, list_pair &next, 
      int grainsize, int level);
#elif method == 2
  void process_queue(list_pair &queue, list_pair &next, int grainsize,
      int level, int strands, int first_strand, int last_strand);
  int count_strands(int grainsize, int size);
#endif
};

#if method == 0
void BFS::run() {
  queue.push_back(make_pair(-1, root));
  node_queued[root] = true;
  for(int level=0; !queue.empty(); level++) {
    for(list_pair::iterator edge = queue.begin();
        edge != queue.end(); edge++) {
      int node = edge->second;
      node_level[node] = level;
      node_parent[node] = edge->first;
      for(int i=0; i<graph.degree(node); i++) {
        int neighbor = graph.neighbor(node, i);
        if(!node_queued[neighbor]) {
          node_queued[neighbor] = true;
          next.push_back(make_pair(node, neighbor));
        }
      }
    }
    queue.swap(next);
    next.clear();
  }
}
#elif method == 1
void BFS::run() {
  list_pair queue, next;
  queue.push_back(make_pair(-1, root));
  node_queued[root] = true;
  for(int level=0; !queue.empty(); level++) {
    int grainsize = min((long unsigned) 2048,
        queue.size() / (8*cilk::current_worker_count()));
    process_queue(queue, next, grainsize, level);
    queue.swap(next);
    next.clear();
  }
}

void BFS::process_queue(list_pair &queue, list_pair &next, int grainsize, int level) {
  if(queue.size() <= grainsize || queue.size() <= 1) {
    for(list_pair::iterator edge = queue.begin();
        edge != queue.end(); edge++) {
      int node = edge->second;
      node_level[node] = level;
      node_parent[node] = edge->first;
      for(int i=0; i<graph.degree(node); i++) {
        int neighbor = graph.neighbor(node, i);
        if(!node_queued[neighbor]) {
          node_queued[neighbor] = true;
          next.push_back(make_pair(node, neighbor));
        }
      }
    }
  } else {
    // Find the middle.
    list_pair::iterator middle = queue.begin();
    advance(middle, queue.size() / 2);

    // Split the queue.
    list_pair left, right;
    left.splice(left.end(), queue, queue.begin(), middle);
    right.splice(right.end(), queue);

    // Run the job!
    list_pair left_next, right_next;
    cilk_spawn process_queue(left, left_next, grainsize, level);
    /*      */ process_queue(right, right_next, grainsize, level);
    cilk_sync;

    // Join the queues.
    next.splice(next.end(), left_next);
    next.splice(next.end(), right_next);
  }
}

#elif method == 2
void BFS::run() {
  list_pair queue, next;
  queue.push_back(make_pair(-1, root));
  node_queued[root] = true;
  for(int level=0; !queue.empty(); level++) {
    int grainsize = min((long unsigned) 2048,
        queue.size() / (8*cilk::current_worker_count()));
    int strands = count_strands(grainsize, queue.size());
    process_queue(queue, next, grainsize, level, strands, 0, strands-1);
    queue.swap(next);
    next.clear();
  }
}

int BFS::count_strands(int grainsize, int size) {
  if(size <= grainsize || size <= 1) {
    return 1;
  } else {
    return count_strands(grainsize, size/2) +
      count_strands(grainsize, size - size/2);
  }
}

void BFS::process_queue(list_pair &queue, list_pair &next, int grainsize,
    int level, int strands, int first_strand, int last_strand) {
  if(first_strand == last_strand) {
    list_pair::iterator edge = queue.begin();
    for(int i=0; i<first_strand; i++) edge++;
    while(edge != queue.end()) {
      int node = edge->second;
      node_level[node] = level;
      node_parent[node] = edge->first;
      for(int i=0; i<graph.degree(node); i++) {
        int neighbor = graph.neighbor(node, i);
        if(!node_queued[neighbor]) {
          node_queued[neighbor] = true;
          next.push_back(make_pair(node, neighbor));
        }
      }
      for(int i=0; i<strands && edge != queue.end(); i++) edge++;
    }
  } else {
    // Find the middle: NOOP.
    // Split the queue: NOOP.

    // Run the job!
    int middle_strand = (first_strand + last_strand) / 2;
    list_pair left_next, right_next;
    cilk_spawn process_queue(queue, left_next, grainsize, level, 
        strands, first_strand, middle_strand);
    /*      */ process_queue(queue, right_next, grainsize, level,
        strands, middle_strand + 1, last_strand);
    cilk_sync;

    // Join the queues.
    next.splice(next.end(), left_next);
    next.splice(next.end(), right_next);
  }
}
#endif

bool BFS::validate() {
  set<int> needs_parent;
  validation_failed = false;

  for(int node=0; node<graph.size(); node++) {
    int parent = node_parent[node];
    if(parent == -1) continue;

    validate_bfs_node_points_root(node); // (1)
    validate_bfs_edge_level(node, parent); // (2)
    validate_bfs_edge_in_graph(node, parent); // (5)
  }

  for(int node=0; node<graph.size(); node++) {
    for(int i=0; i<graph.degree(node); i++) {
      int child = graph.neighbor(node, i);
      validate_graph_edge_level(child, node); // (3)
    }
  }

  validate_graph_edges_in_bfs(root); // (4)
  return !validation_failed;
}

// (1) the BFS tree is a tree and does not contain cycles
void BFS::validate_bfs_node_points_root(int node) {
  list<int> offspring;
  validate_bfs_node_points_root(node, offspring);
}

// Offspring contains one of the current node's children as the last element,
// one of the current node's grandchildren as the second to last element, etc.
void BFS::validate_bfs_node_points_root(int node, list<int> &offspring) {
  static vector<bool> checked;
  checked.resize(graph.size(), false);
  
  if(node == -1 && offspring.back() != root) {
      cerr << "FAILED: node " << offspring.back() << 
        " has no parent but is not root (" << root << ")" << endl;
      validation_failed = true;
  } else if(node == -1 || checked[node]) {
    for(list<int>::iterator it = offspring.begin(); it != offspring.end(); it++)
      checked[*it] = true;
  } else {
    list<int>::iterator it = find(offspring.begin(), offspring.end(), node);
    if(it != offspring.end()) {
      cerr << "FAILED: found a cycle ";
      for(; it != offspring.end(); it++) cerr << (*it) << " -> ";
      cerr << node << endl;
      validation_failed = true;
    } else {
      offspring.push_back(node);
      validate_bfs_node_points_root(node_parent[node], offspring);
    }
  }
}

// (2) each tree edge connects vertices whose BFS levels differ by exactly one
void BFS::validate_bfs_edge_level(int child, int parent) {
  int child_level = node_level[child];
  int parent_level = node_level[parent];
  if(child_level - parent_level != 1) {
    cerr << "FAILED: parent of node " << child << " (level " <<
      child_level << ") is " << parent << " and has level " << 
      parent_level << ", expected level " << child_level - 1 << endl;
    validation_failed = true;
  }
}

// (3) every edge in the input list has vertices with levels that differ by at
// most one or that both are not in the BFS tree
void BFS::validate_graph_edge_level(int child, int parent) {
  int child_level = node_level[child];
  int parent_level = node_level[parent];
  if(parent_level == -1) return;
  if(child_level - parent_level > 1) {
    cerr << "FAILED: node " << child << " has parent " << parent <<
      " with level " << parent_level << " but has level " << child_level <<
      ", expected level " << parent_level + 1 << " or lower" << endl;
    validation_failed = true;
  }
}

// (4) the BFS tree spans an entire connected component's vertices
void BFS::validate_graph_edges_in_bfs(int parent) {
  static vector<bool> checked;
  checked.resize(graph.size(), false);

  if(checked[parent]) return;
  for(int i=0; i<graph.degree(parent); i++) {
    int node = graph.neighbor(parent, i);
    if((node != root && node_parent[node] == -1) && node_level[node] == -1) {
      cerr << "FAILED: node " << node << " is reachable in the original " <<
        "graph through parent node " << parent << " but no parent or level information " <<
        "was found in the BFS tree" << endl;
      validation_failed = true;
    } else if(node != root && node_parent[node] == -1) {
      cerr << "FAILED: node " << node << " is reachable in the original " <<
        "graph through parent node " << parent << " but no parent information " <<
        "was found in the BFS tree" << endl;
      validation_failed = true;
    } else if(node_level[node] == -1) {
      cerr << "FAILED: node " << node << " is reachable in the original " <<
        "graph through parent node " << parent << " but no level information " <<
        "was found in the BFS tree" << endl;
      validation_failed = true;
    }
    checked[parent] = true;

    validate_graph_edges_in_bfs(node);
  }
}

// (5) a node and its parent are joined by an edge of the original graph
void BFS::validate_bfs_edge_in_graph(int node, int parent) {
  bool found = false;
  for(int i=0; i<graph.degree(parent); i++) {
    if(node == graph.neighbor(parent, i)) {
      found = true;
      break;
    }
  }
  if(!found) {
    cerr << "FAILED: parent of node " << node << " is supposedly " <<
      parent << " but there is no such edge in the original graph" << endl;
    validation_failed = true;
  }
}

void BFS::print_statistics() {
  cout << "Search reached ";
  int levels = *max_element(node_level.begin(), node_level.end()) + 1;
  cout << levels << " levels and ";
  int nodes = node_level.size() -
    count(node_level.begin(), node_level.end(), -1);
  cout << nodes << " vertices" << endl;

  double average_degree = 0;
  for(int i=0; i<node_level.size(); i++)
    if(node_level[i] != -1) average_degree += graph.degree(i);
  average_degree /= nodes;
  cout << "Average degree of reached nodes " <<
    "(including duplicate edges) is " << average_degree << endl;

  for(int i=0; i<levels; i++) {
    cout << "level " << i << " vertices: " <<
      count(node_level.begin(), node_level.end(), i) << endl;
  }
}

void BFS::print_results() {
  printf("\n  vertex parent  level\n");
  for(int i=0; i<node_level.size(); i++) {
    printf("%6d%7d%7d\n", i, node_parent[i], node_level[i]);
  }
  cout << endl;
}


int cilk_main(int argc, char** argv) {
  if(argc < 2 || argc > 4) {
    cerr << "Usage: " << argv[0] << "[FILE] NODE [-v]" << endl;
    return 1;
  }

  string filename; int node; int verbose = 0;
  if(argv[1][0] >= '0' && argv[1][0] <= '9') {
    filename = "-";
    node = atoi(argv[1]);
    if(argc == 3 && argv[2][0] == '-' && argv[2][1] == 'v')
      verbose = strlen(argv[2]) - 1;
  } else {
    filename = argv[1];
    node = atoi(argv[2]);
    if(argc == 4 && argv[3][0] == '-' && argv[3][1] == 'v')
      verbose = strlen(argv[3]) - 1;
  }

  Graph graph(filename);
  graph.print();
  if(graph.size() < node+1) {
    cerr << "The graph does not contain starting vertex " << node << endl;
    exit(1);
  }
  cout << "Starting vertex: " << node << endl << endl;
  BFS bfs(graph, node);

  cilk::cilkview cv;
  cv.start();
  bfs.run();
  cv.stop();
  cout << "Time: " << cv.accumulated_milliseconds() << endl;
  cv.dump("parallel.profile");

  if(verbose >= 1) {
    cout << endl;
    bfs.print_statistics();
    if(verbose >= 2) bfs.print_results();
    cout << endl;
  }

  bool success = bfs.validate();
  if(success) {
    cout << "Results validated successfully" << endl;
    return 0;
  } else {
    cout << "Results were invalid" << endl;
    return 1;
  }
}

