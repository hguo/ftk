#include "FeatureTransitionMatrix.h"
#include <sstream>
#include <cstdio>
#include <climits>
#include <cmath>
#include <cassert>

FeatureTransitionMatrix::FeatureTransitionMatrix() :
  _n0(INT_MAX), _n1(INT_MAX)
{
}

FeatureTransitionMatrix::FeatureTransitionMatrix(int n0, int n1) :
  _n0(n0), _n1(n1)
{
  _match.resize(n0*n1);
}

FeatureTransitionMatrix::FeatureTransitionMatrix(int t0, int t1, int n0, int n1) :
  _interval(std::make_pair(t0, t1)), 
  _n0(n0), _n1(n1)
{
  _match.resize(_n0*_n1);
}

FeatureTransitionMatrix::~FeatureTransitionMatrix()
{
}

#if 0
std::string FeatureTransitionMatrix::MatrixFileName(const std::string& dataname, int t0, int t1) const
{
  std::stringstream ss;
  ss << dataname << ".match." << t0 << "." << t1;
  return ss.str();
}
#endif

int FeatureTransitionMatrix::operator()(int i, int j) const
{
  return _match[i*n1() + j];
}

int& FeatureTransitionMatrix::operator()(int i, int j)
{
  return _match[i*n1() + j];
}

int FeatureTransitionMatrix::at(int i, int j) const
{
  return _match[i*n1() + j];
}

int& FeatureTransitionMatrix::at(int i, int j)
{
  return _match[i*n1() + j];
}

int FeatureTransitionMatrix::colsum(int j) const
{
  int sum = 0;
  for (int i=0; i<n0(); i++)
    sum += _match[i*n1() + j];
  return sum;
}

int FeatureTransitionMatrix::rowsum(int i) const 
{
  int sum = 0;
  for (int j=0; j<n1(); j++) 
    sum += _match[i*n1() + j];
  return sum;
}

void FeatureTransitionMatrix::GetModule(int i, std::set<int> &lhs, std::set<int> &rhs, int &event) const
{
  if (i>=_lhss.size()) return;

  lhs = _lhss[i];
  rhs = _rhss[i];
  event = _events[i];
}

void FeatureTransitionMatrix::Normalize()
{
  if (NModules() == 0) Modularize();
  for (int i=0; i<NModules(); i++) {
    const std::set<int> &lhs = _lhss[i];
    const std::set<int> &rhs = _rhss[i];
    for (std::set<int>::const_iterator it0=lhs.begin(); it0!=lhs.end(); it0++) 
      for (std::set<int>::const_iterator it1=rhs.begin(); it1!=rhs.end(); it1++) {
        int l = *it0, r = *it1;
        at(l, r) = 1;
      }
  }
}

void FeatureTransitionMatrix::Modularize()
{
  const int n = n0() + n1();

  _lhss.clear();
  _rhss.clear();
  _events.clear();

  std::set<int> unvisited; 
  for (int i=0; i<n; i++) 
    unvisited.insert(i);

  while (!unvisited.empty()) {
    std::set<int> lhs, rhs; 
    std::vector<int> Q;
    Q.push_back(*unvisited.begin());

    while (!Q.empty()) {
      int v = Q.back();
      Q.pop_back();
      unvisited.erase(v);
      if (v<n0()) lhs.insert(v);
      else rhs.insert(v-n0());

      if (v<n0()) { // left hand side
        for (int j=0; j<n1(); j++) 
          if (at(v, j)>0 && unvisited.find(j+n0()) != unvisited.end())
            Q.push_back(j+n0());
      } else { // right hand side
        for (int i=0; i<n0(); i++) 
          if (at(i, v-n0())>0 && unvisited.find(i) != unvisited.end())
            Q.push_back(i);
      }
    }

    int event; 
    if (lhs.size() == 1 && rhs.size() == 1) {
      event = FEATURE_EVENT_DUMMY;
    } else if (lhs.size() == 0 && rhs.size() == 1) {
      event = FEATURE_EVENT_BIRTH;
    } else if (lhs.size() == 1 && rhs.size() == 0) {
      event = FEATURE_EVENT_DEATH;
    } else if (lhs.size() == 1 && rhs.size() == 2) {
      event = FEATURE_EVENT_SPLIT;
    } else if (lhs.size() == 2 && rhs.size() == 1) { 
      event = FEATURE_EVENT_MERGE;
    } else if (lhs.size() == 2 && rhs.size() == 2) { 
      event = FEATURE_EVENT_RECOMBINATION;
    } else {
      event = FEATURE_EVENT_COMPOUND;
    }

    _lhss.push_back(lhs);
    _rhss.push_back(rhs);
    _events.push_back(event);
  }

  moving_speeds.resize(_n0, NAN);
}

void FeatureTransitionMatrix::Print() const
{
  fprintf(stderr, "Interval={%d, %d}, n0=%d, n1=%d\n", 
      _interval.first, _interval.second, _n0, _n1);

  for (int i=0; i<_n0; i++) {
    for (int j=0; j<_n1; j++) {
      fprintf(stderr, "%d\t", at(i, j));
    }
    fprintf(stderr, "\n");
  }
}

void FeatureTransitionMatrix::SaveAscii(const std::string& filename) const 
{
  FILE *fp = fopen(filename.c_str(), "w");
  
  fprintf(fp, "Interval={%d, %d}, n0=%d, n1=%d\n", 
      _interval.first, _interval.second, _n0, _n1);

  for (int i=0; i<_n0; i++) {
    for (int j=0; j<_n1; j++) {
      fprintf(fp, "%d\t", at(i, j));
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}
