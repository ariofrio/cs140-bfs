#include <cilk.h>
#include <cilkview.h>
#include <algorithm>
#include <iostream>
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
  Graph();
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

  vector< pair<int, int> > edges;

private:
  vector<int> neighbors;
  vector<int> indexes;

  void read_edges();
  void build_graph();
};

Graph::Graph() {
  read_edges();
  build_graph();
}

void Graph::read_edges() {
  while(!cin.eof()) {
    pair<int, int> edge;
    cin >> edge.first;
    cin >> edge.second;
    if(cin.fail()) break;
    edges.push_back(edge);
  }
}

void Graph::build_graph() {
  for(int i=0; i<edges.size(); i++) {
    if(indexes.size() < edges[i].first+1)
      indexes.resize(edges[i].first+1);
    if(indexes.size() < edges[i].second+1)
      indexes.resize(edges[i].second+1);
  }

  for(int i=0; i<edges.size(); i++)
    indexes.at(edges[i].first)++;
  for(int i=1; i<indexes.size(); i++)
    indexes[i] += indexes[i-1];

  neighbors.resize(indexes.back()+1);
  
  for(int i=0; i<edges.size(); i++) {
    int index = --indexes.at(edges[i].first);
    neighbors.at(index) = edges[i].second;
  }
}

void Graph::print() {
  printf("\nGraph has %d vertices and %d edges\n",
      indexes.size(), edges.size());
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
  void print();
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

  for(vector< pair<int, int> >::iterator it = graph.edges.begin();
      it != graph.edges.end(); it++) {
    int parent = it->first;
    int child = it->second;
    validate_graph_edge_level(child, parent); // (3)
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

void BFS::print() {
  int levels = *max_element(node_level.begin(), node_level.end()) + 1;
  int nodes = node_level.size() -
    count(node_level.begin(), node_level.end(), -1);
  cout << "Search reached "
    << levels << " levels and "
    << nodes << " vertices" << endl;
  for(int i=0; i<levels; i++) {
    cout << "level " << i << " vertices: " <<
      count(node_level.begin(), node_level.end(), i) << endl;
  }
  printf("\n  vertex parent  level\n");
  for(int i=0; i<node_level.size(); i++) {
    printf("%6d%7d%7d\n", i, node_parent[i], node_level[i]);
  }
  cout << endl;
}

int cilk_main(int argc, char** argv) {
  if(argc != 2 && argc != 3) {
    cerr << "Usage: ./parallel NODE [-v]" << endl;
    return 1;
  }
  int node = atoi(argv[1]);

  Graph graph;
  BFS bfs(graph, node);
  graph.print();
  cout << "Starting vertex: " << node << endl << endl;

  cilk::cilkview cv;
  cv.start();
  bfs.run();
  cv.stop();
  cv.dump("parallel.profile");

  if(argc == 3 && strcmp(argv[2], "-v") == 0) {
    bfs.print();
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

