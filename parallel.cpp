#include <cilk.h>
#include <cilkview.h>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <vector>
#include <list>
using namespace std;

typedef list< pair<int, int> > list_pair;

class Graph {
public:
  Graph();
  void print();
  int degree(int node) {
    if(node == indexes.size() - 1)
      return indexes.size() - indexes[node];
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
  vector< pair<int, int> > edges;
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
  BFS(Graph& graph) : graph(graph) {
    node_level.resize(graph.size(), -1);
    node_parent.resize(graph.size(), -1);
    node_todo.resize(graph.size(), false);
  };

  vector<int> node_level;
  vector<int> node_parent;
  vector<bool> node_todo;
  void run(int root);
  void print();

private:
#if method == 0
  list_pair queue;
  list_pair next;
#elif method == 1
  void process_queue(list_pair &queue, list_pair &next, 
      int grainsize, int level);
#endif
};

#if method == 0
void BFS::run(int root) {
  queue.push_back(make_pair(-1, root));
  for(int level=0; !queue.empty(); level++) {
    for(list_pair::iterator edge = queue.begin();
        edge != queue.end(); edge++) {
      int node = edge->second;
      node_level[node] = level;
      node_parent[node] = edge->first;
      for(int i=0; i<graph.degree(node); i++) {
        int neighbor = graph.neighbor(node, i);
        if(!node_todo[neighbor]) {
          node_todo[neighbor] = true;
          next.push_back(make_pair(node, neighbor));
        }
      }
    }
    queue.swap(next);
    next.clear();
  }
}
#elif method == 1
void BFS::run(int root) {
  list_pair queue, next;
  queue.push_back(make_pair(-1, root));
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
        if(!node_todo[neighbor]) {
          next.push_back(make_pair(node, neighbor));
          node_todo[neighbor] = true;
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
    process_queue(grainsize);

void BFS::process_queue(int step,
    list_pair::iterator begin,
    list_pair::iterator end) {
  if(begin == end) {
    queue.
  }
}
#endif

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
  if(argc != 2) {
    cerr << "Usage: ./parallel NODE" << endl;
    return 1;
  }
  int node = atoi(argv[1]);

  Graph graph;
  BFS bfs(graph);
  graph.print();
  cout << "Starting vertex: " << node << endl << endl;

  cilk::cilkview cv;
  cv.start();
  bfs.run(node);
  cv.stop();

  bfs.print();
  cv.dump("parallel.profile");

  return 0;
}

