#include "GraphColor.h"
#include <algorithm>
#include <cstdlib>

typedef struct {
  int index; 
  int deg; 
  int color;
} Node;

bool compare(Node n0, Node n1) {
  return n0.deg > n1.deg;
}

bool compare1(Node n0, Node n1) {
  return n0.index < n1.index;
}

bool check_ok(int n, bool **M, Node *nodes, int *cid, int i, int j)
{
  if (M[nodes[i].index][nodes[j].index] != 0 || cid[j] != 0) return false;
  for (int k=0; k<n; k++) {
    if (M[nodes[i].index][nodes[k].index] == 0 && M[nodes[k].index][nodes[j].index] != 0)
      return false;
  }
  return true;
}

int welsh_powell(int n, bool **M, int *cid)
{
  Node *nodes = (Node*)malloc(sizeof(Node)*n);
  for (int i=0; i<n; i++) {
    int deg = 0;
    for (int j=0; j<n; j++) 
      if (M[i][j]) deg ++;
    nodes[i].index = i;
    nodes[i].deg = deg;
    nodes[i].color = 0;
  }

  std::stable_sort(nodes, nodes+n, compare);
  int k=0;
  
  while (1) {
    k++;
    int i;
    for (i=0; i<n; i++) {
      if (nodes[i].color == 0) {
        nodes[i].color = k;
        break;
      }
    }
    if (i==n) break;
    for (int j=0; j<n; j++)
      if (i!=j && check_ok(n, M, nodes, cid, i, j))
        nodes[j].color = k;
  }

  std::stable_sort(nodes, nodes+n, compare1);
  for (int i=0; i<n; i++) 
    cid[i] = nodes[i].color - 1;

  return k-1;
}
