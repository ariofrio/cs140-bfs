#ifdef __cilkplusplus
#include <cilk.h>
#include <cilkview.h>
#define current_worker_count() cilk::current_worker_count()
#else
#define CILK_STUB 1

/* Compile away the Cilk++ keywords. */
#define __cilk
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#define cilk_run

/* cilk_main and cilk_wmain are not keywords, but are more-or-less treated as them. */
#define cilk_main main
#define cilk_wmain wmain
#define current_worker_count() 4
#endif

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstring>
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
      return size_edges() - indexes[node];
    else
      return indexes[node+1] - indexes[node];
  }
  int neighbor(int node, int index) {
    return neighbors[indexes[node] + index];
  }
  int order() {
    return indexes.size();
  }
  int size_edges() {
    return neighbors.size();
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
      order(), size_edges());
}

#if method == 3
template<class T>
class bag {
 
  // A node inside a pennant. If this node doesn't have a right child and its
  // left subtree is a perfect complete binary tree (or is nonexistent), then
  // this node is also a pennant. This means that a pennant contains exactly
  // 2^k nodes (including the root).
  class node {
  public:

    class iterator : public std::iterator<input_iterator_tag, T> {
    public:
      iterator(node* x) : x(x) { if(x != NULL) x->parent = NULL; }
      iterator(const iterator& other) : x(other.x) {}
      T& operator* () { return **x; }
      T* operator-> () { return &(**x); }
      bool operator==(const iterator& other) { return x == other.x; }
      bool operator!=(const iterator& other) { return x != other.x; }
      iterator operator++(int) { iterator temp(*this); ++(*this); return temp; }
      iterator& operator++() {
        if(x->left != NULL) {
          x->left->parent = x;
          x = x->left;
        } else if(x->right != NULL) {
          x->right->parent = x;
          x = x->right;
        } else {
          node* previous; do {
            previous = x;
            x = x->parent;
            if(x == NULL) return *this;
          } while(previous == x->right);
          if(x->right != NULL)
            x->right->parent = x;
          x = x->right;
        }
        return *this;
      }
    private:
      node* x;
    };

    node(T value) : value(value), left(NULL), right(NULL) {}
    T& operator*() { return value; }
    
    // Merges a pennant into this pennant. After this call, the other node
    // is no longer a pennant and should not be used as one.
    //
    // Iterators: all become invalid.
    //
    // For convenience, returns this pennant.
    node* merge(node* other) {
      other->right = left;
      left = other;
      return this;
    }

    // Splits off half the elements of this pennant, assuming this pennant
    // contains at least 2 elements.
    //
    // Iterators: all become invalid.
    node* split() {
      node* other = left;
      left = other->right;
      other->right = NULL;
      return other;
    }

    iterator begin() { return iterator(this); }
    iterator end() { return iterator(NULL); }

    friend void bag<T>::clear_pennant_bare(node* pennant);

  private:
    T value;
    node* left;
    node* right;
    node* parent; // not always accurate, used by iterator

  };

public:

  class iterator : public std::iterator<input_iterator_tag, T> {
  public:
    iterator(typename vector<node*>::iterator b, typename vector<node*>::iterator end)
      : b(b), end(end), p(NULL) {
      while(b != end && *b == NULL) b++;
      if(b != end) p = (*b)->begin();
    }
    iterator(const iterator& other) : b(other.b), end(other.end), p(other.p) {}
    T& operator* () { return *p; }
    T* operator-> () { return &(*p); }
    bool operator==(const iterator& other) { return b != end ? p == other.p : other.b == end; }
    bool operator!=(const iterator& other) { return b != end ? p != other.p : other.b != end; }
    iterator operator++(int) { iterator temp(*this); ++(*this); return temp; }
    iterator& operator++() {
      if(++p == (*b)->end()) {
        do { b++; } while(b != end && *b == NULL);
        if(b != end) p = (*b)->begin();
      }
      return *this;
    }
  private:
    typename vector<node*>::iterator b;
    typename vector<node*>::iterator end;
    typename node::iterator p;
  };

  bag(size_t capacity) : _capacity(capacity), _size(0) {
    int r = 0; while(pow(double(2), r+1) <= capacity) r++;
    backbone.resize(r+1, NULL);
  }
  ~bag() { clear(); }

  size_t capacity() { return _capacity; }
  size_t size() { 
    //return _size;
    size_t ret = 0;
    for(int k=0; k<backbone.size(); k++) ret |= (backbone[k] != NULL) << k;
    return ret;
  }
  bool empty() { return _size == 0; }

  void insert(T value) {
    insert(new node(value));
  }

  // Merges a bag into this bag. After this call, the other bag is left empty.
  void merge(bag<T>& other) {
    node* y = NULL; // the "carry" bit
    for(int k=0; k<backbone.size(); k++) {
      full_adder(backbone[k], y, backbone[k], other.backbone[k], y);
      other.backbone[k] = NULL;
    }
    _size += other._size;
    other._size = 0;
  }

  // Splits off half of the elements of this bag into another empty bag.
  void split(bag<T>& other) {
    node* y = backbone[0];
    backbone[0] = NULL;
    for(int k=1; k<backbone.size(); k++) {
      if(backbone[k] != NULL) {
        other.backbone[k-1] = backbone[k]->split(); // safe since k >= 1
        backbone[k-1] = backbone[k];
        backbone[k] = NULL;
      }
    }
    _size >>= 1;
    other._size = _size;
    if(y != NULL) insert(y);
  }

  void swap(bag<T>& x) {
    std::swap(backbone, x.backbone); // constant time
    std::swap(_capacity, x._capacity);
    std::swap(_size, x._size);
  }

  void clear() {
    for(int k=0; k<backbone.size(); k++)
      clear_pennant(backbone[k]);
    _size = 0;
  }

  iterator begin() { return iterator(backbone.begin(), backbone.end()); }
  iterator end() { return iterator(backbone.end(), backbone.end()); }

private:
  // Each entry backbone[k] contains either a NULL pointer or a pointer to a
  // pennant of size 2^k.
  vector<node*> backbone;
  size_t _capacity;
  size_t _size;

  void insert(node* x) {
    int k = 0;
    while(backbone[k] != NULL) {
      x = backbone[k]->merge(x);
      backbone[k++] = NULL;
    }
    backbone[k] = x;
    _size++;
  }

  void full_adder(node* &s, node* &c, node* x, node* y, node* z) {
    char bits = (x != NULL) << 6 | (y != NULL) << 3 | (z != NULL);
    switch(bits) {
      case 0000: s = NULL; c = NULL; break;
      case 0100: s = x; c = NULL; break;
      case 0010: s = y; c = NULL; break;
      case 0001: s = z; c = NULL; break;
      case 0110: s = NULL; c = x->merge(y); break;
      case 0101: s = NULL; c = x->merge(z); break;
      case 0011: s = NULL; c = y->merge(z); break;
      case 0111: s = x; c = y->merge(z); break;
    }
  }

  void clear_pennant(node* &pennant) {
    clear_pennant_bare(pennant);
    pennant = NULL;
  }

  void clear_pennant_bare(node* pennant) {
    if(pennant == NULL) return;
    clear_pennant_bare(pennant->left);
    clear_pennant_bare(pennant->right);
    delete pennant;
  }

};

template<class T> void std::swap(bag<T> &a, bag<T> &b) {
  a.swap(b);
}

typedef bag< pair<int, int> > bag_pair;
#endif

class BFS {
public:
  Graph& graph;
  int root;
  BFS(Graph& graph, int root) : graph(graph), root(root) {
    node_level.resize(graph.order(), -1);
    node_parent.resize(graph.order(), -1);
    node_queued.resize(graph.order(), false);
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
#elif method == 3
  void process_queue(bag_pair &queue, bag_pair &next, 
      int grainsize, int level);
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
        queue.size() / (8*current_worker_count()));
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
        queue.size() / (8*current_worker_count()));
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
#elif method == 3
void BFS::run() {
  bag_pair queue(graph.size_edges()), next(graph.size_edges());
  queue.insert(make_pair(-1, root));
  node_queued[root] = true;
  for(int level=0; !queue.empty(); level++) {
    int grainsize = min((long unsigned) 2048,
        queue.size() / (8*current_worker_count()));
    cout << "level " << level << " has queue.size()=" << queue.size();
    cout << " and grain size " << grainsize << endl;
    process_queue(queue, next, grainsize, level);
    queue.swap(next);
    next.clear();
  }
}

void BFS::process_queue(bag_pair &queue, bag_pair &next, int grainsize, int level) {
  if(queue.size() <= grainsize || queue.size() <= 1) {
    cout << "processing...";
    for(bag_pair::iterator edge = queue.begin(); edge != queue.end(); edge++) {
      int node = edge->second;
      node_level[node] = level;
      node_parent[node] = edge->first;
      for(int i=0; i<graph.degree(node); i++) {
        int neighbor = graph.neighbor(node, i);
        if(!node_queued[neighbor]) {
          node_queued[neighbor] = true;
          next.insert(make_pair(node, neighbor));
        }
      }
    }
    cout << "DONE" << endl;
  } else {
    cout << "splitting...";
    // Split the queue.
    bag_pair right_queue(graph.size_edges());
    queue.split(right_queue);
    cout << "DONE" << endl;

    // Run the job!
    bag_pair right_next(graph.size_edges());
    /*cilk_spawn*/ process_queue(queue, next, grainsize, level);
    /*      */ process_queue(right_queue, right_next, grainsize, level);
    //cilk_sync;

    // Join the queues.
    next.merge(right_next);
  }
}

#endif

bool BFS::validate() {
  set<int> needs_parent;
  validation_failed = false;

  for(int node=0; node<graph.order(); node++) {
    int parent = node_parent[node];
    if(parent == -1) continue;

    validate_bfs_node_points_root(node); // (1)
    validate_bfs_edge_level(node, parent); // (2)
    validate_bfs_edge_in_graph(node, parent); // (5)
  }

  for(int node=0; node<graph.order(); node++) {
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
  checked.resize(graph.order(), false);
  
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
  checked.resize(graph.order(), false);

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
  if(graph.order() < node+1) {
    cerr << "The graph does not contain starting vertex " << node << endl;
    exit(1);
  }
  cout << "Starting vertex: " << node << endl << endl;
  BFS bfs(graph, node);

#ifdef __cilkplusplus
  cilk::cilkview cv;
  cv.start();
#endif
  bfs.run();
#ifdef __cilkplusplus
  cv.stop();
  cout << "Time: " << cv.accumulated_milliseconds() << endl;
  cv.dump("parallel.profile");
#endif

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

