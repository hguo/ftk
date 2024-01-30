#ifndef _FTK_KD_LITE_HH
#define _FTK_KD_LITE_HH

// kd impl that can be used in a cuda kernel

#include <ftk/config.hh>
#include <ftk/numeric/vector_norm.hh>

namespace ftk {

// implementation of https://www.nvidia.com/content/gtc-2010/pdfs/2140_gtc2010.pdf
template <int nd, typename I=int, typename F=double>
__host__
inline void kdlite_build_recursive(
    const I n,
    const I current,
    const F *X, // coordinates
    const I level, // the current level
    const I offset, // the current offset, 0 for the root pass
    const I length, // the current length
    I *heap, // out: heap
    I *ids) // out: pre-allocated array for ids
{
  const I axis = level % nd;
  const I h = std::ceil(std::log2(length+1));
  const I half = std::pow(2, h-2);
  const I lastRow = length - 2 * half + 1;
  const I lbm = half + std::min(half, lastRow);
  // const I half = length / 2;
  
  // fprintf(stderr, "current=%d, offset=%d, length=%d\n", 
  //     current, offset, length);

  if (length == 1) {
    heap[current] = ids[offset];
    // fprintf(stderr, "current=%d, offset=%d, length=%d, lbm=%d, median=%d\n", current, offset, length, lbm, heap[current]);
  } else {
    std::nth_element(
        ids + offset, 
        ids + offset + lbm-1, // median
        ids + offset + length, 
        [X, axis](I i, I j) {
          return X[i*nd+axis] < X[j*nd+axis];
        });

    heap[current] = ids[offset + lbm-1]; // the median
    // fprintf(stderr, "current=%d, offset=%d, length=%d, lbm=%d, median=%d\n", current, offset, length, lbm, heap[current]);

    if (lbm - 1 >= 1) 
      kdlite_build_recursive<nd, I, F>(n, current*2+1, X, level+1, offset,      lbm-1,     heap, ids); // left
    if (length - lbm >= 1)  
      kdlite_build_recursive<nd, I, F>(n, current*2+2, X, level+1, offset+lbm, length-lbm, heap, ids); // right
  }
}

template <int nd, typename I=int, typename F=double>
__host__
inline void kdlite_build(
    const I n, // number of points
    const F *X, // coordinates
    I *heap) // out: pre-allocated heap
{
  std::vector<I> ids;
  ids.resize(n);
  for (int i = 0; i < n; i ++)
    ids[i] = i;

  kdlite_build_recursive<nd, I, F>(n, 0, X, 0, 0, n, heap, ids.data());

  // for (int i = 0; i < n; i ++) 
  //   fprintf(stderr, "i=%d, heap=%d\n", i, heap[i]);
}

template <int nd, typename I=int, typename F=double>
__device__ __host__
inline I kdlite_nearest(I n, const F *X, const I *heap, const F *x)
{
  static constexpr size_t max_stack_size = 32; // TODO

  I S[max_stack_size];
  I top = 0;

  S[top++] = 0; // push root //  S[top].depth = 0; // root // depth = log2(i+1);
  
  I best = -1; // no best yet
  F best_d2 = 1e32; // no best distance yet

  while (top != 0) { // stack is not empty
    const I i = S[--top]; // pop stack

    const I xid = heap[i];
#ifdef __CUDACC__
    const I depth = 31 - __clz(i+1); // see https://forums.developer.nvidia.com/t/integer-logarithm-to-the-base-2-how/3722/6
#else
    const I depth = std::log2(i+1);
#endif
    const I axis = depth % nd;
    I next, other;

    if (x[axis] < X[nd*xid+axis]) {
      next = i * 2 + 1; // left child
      other = i * 2 + 2; // right child
    } else {
      next = i * 2 + 2; // right child
      other = i * 2 + 1; // left child
    }

    const F d2 = vector_dist_2norm2<F>(nd, x, X + nd*xid); // distance to the current node
    if (d2 < best_d2) {
      best = xid;
      best_d2 = d2;

      // fprintf(stderr, "current_best=%d, d2=%f, X=%f, %f, %f\n", 
      //     best, best_d2, 
      //     X[nd*xid], X[nd*xid+1], X[nd*xid+2]);
    }
      
    // const F dp = x[axis] - X[nd*xid+axis]; // distance to the median
    // const F dp2 = dp * dp;

    if (next < n) { // the next node exists
      assert(top < max_stack_size);
      S[top++] = next; // push stack
    }

    if (other < n) {
      const F dp = x[axis] - X[nd*xid+axis];
      const F dp2 = dp * dp;

      if (dp2 <= best_d2) {
        assert(top < max_stack_size);
        S[top++] = other; // push stack
      }
    }
  }

  return best;
}

template <int nd, typename I=int, typename F=double>
__device__ __host__
I kdlite_nearest_bfs(I n, const F *X, const I *heap, const F *x)
{
  static size_t max_queue_size = 32768; // TODO
  typedef struct {
    I i, depth;
  } qnode_t;

  qnode_t Q[max_queue_size];
  I iq = 0, jq = 0; // start and end of the queue

  Q[jq].i = 0; Q[jq].depth = 0;  // push the root node
  jq = (jq+1) % max_queue_size;
  
  I best = -1; // no best yet
  F best_d2 = 1e32; // no best distance yet

  while (iq != jq) { // Q is not empty
    // fprintf(stderr, "iq=%d, jq=%d\n", iq, jq);

    qnode_t current = Q[iq]; // Q.pop
    iq = (iq+1) % max_queue_size;

    const I i = current.i;
    const I xid = heap[i];
    // fprintf(stderr, "current i %d, xid=%d\n", i, xid);
    const I axis = current.depth % nd;
    I next, other;

    if (x[axis] < X[nd*xid+axis]) {
      next = i * 2 + 1; // left child
      other = i * 2 + 2; // right child
    } else {
      next = i * 2 + 2; // right child
      other = i * 2 + 1; // left child
    }

    const F d2 = vector_dist_2norm2<F>(nd, x, X + nd*xid); // distance to the current node
    if (d2 < best_d2) {
      best = xid;
      best_d2 = d2;

      // fprintf(stderr, "current_best=%d, d2=%f, X=%f, %f, %f\n", 
      //     best, best_d2, 
      //     X[nd*xid], X[nd*xid+1], X[nd*xid+2]);
    }
      
    // const F dp = x[axis] - X[nd*xid+axis]; // distance to the median
    // const F dp2 = dp * dp;

    if (heap[next] >= 0 && next < 2*n) { // the next node exists
      Q[jq].i = next;
      // fprintf(stderr, "adding next %d\n", next);
      Q[jq].depth = current.depth + 1;
      jq = (jq+1) % max_queue_size;
    }

    if (heap[other] >= 0 && other < 2*n) {
      const F dp = x[axis] - X[nd*xid+axis];
      const F dp2 = dp * dp;

      if (dp2 <= best_d2) {
        Q[jq].i = other;
        // fprintf(stderr, "adding other %d\n", other);
        Q[jq].depth = current.depth + 1;
        jq = (jq+1) % max_queue_size;
      }
    }
  }

  return best;
}

} // namespace

#endif
