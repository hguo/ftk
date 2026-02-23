# FTK Thread Safety Documentation

This document describes the thread-safety guarantees and considerations for the FTK (Feature Tracking Kit) library.

## Overview

FTK is designed to support parallel and distributed computation through multiple threading backends and MPI. However, not all components are thread-safe, and users must understand the safety guarantees to use the library correctly in multi-threaded environments.

## Thread Safety Levels

### Level 1: Thread-Safe (Multiple readers AND writers)
These components can be safely accessed from multiple threads simultaneously:

- **`tracking_graph`** (`include/ftk/tracking_graph/tracking_graph.hh`)
  - Uses `std::mutex` for internal synchronization
  - All public methods are thread-safe
  - Safe to add nodes/edges from multiple threads

- **TBB-based trackers** (when built with `FTK_HAVE_TBB`)
  - `xgc_blob_filament_tracker` - Uses `tbb::concurrent_hash_map` and `tbb::concurrent_unordered_set`
  - Concurrent access to feature detection is safe

- **Disjoint Union-Find** (`include/ftk/basic/duf.hh`)
  - Uses `std::mutex` for synchronization
  - Thread-safe union-find operations

### Level 2: Read-Only Thread-Safe (Multiple readers OK, single writer only)
These components are safe for concurrent reading but require external synchronization for writes:

- **Mesh classes** (when not being modified)
  - `simplicial_unstructured_2d_mesh`
  - `simplicial_unstructured_3d_mesh`
  - `simplicial_regular_mesh`
  - `mpas_mesh`
  - `simplicial_xgc_2d_mesh` / `simplicial_xgc_3d_mesh`
  - Safe to query from multiple threads after construction
  - **NOT safe** to modify from multiple threads

- **Feature structures** (when not being modified)
  - `feature_point`
  - `feature_curve`
  - `feature_surface`
  - `feature_curve_set`

### Level 3: Not Thread-Safe (Requires external synchronization)
These components require external synchronization for ANY concurrent access:

- **Tracker field data snapshots**
  - `push_field_data_snapshot()` methods are **NOT thread-safe**
  - Only one thread should call push methods at a time
  - Use external mutex if needed

- **Filter configuration**
  - Setting tracker parameters (thresholds, options, etc.) is **NOT thread-safe**
  - Configure before parallel execution begins

- **I/O operations**
  - File reading/writing methods are **NOT thread-safe**
  - Serialize I/O operations or use thread-local files

## Threading Backends

FTK supports multiple threading backends configured via `set_thread_backend()`:

### FTK_THREAD_PTHREAD (Default)
- Uses standard C++ threads (`std::thread`)
- Portable across all platforms
- Used in `object::parallel_for()`

### FTK_THREAD_OPENMP
- Requires `FTK_HAVE_OPENMP`
- Uses OpenMP parallel for loops
- Good for CPU-bound workloads
- Control thread count with `OMP_NUM_THREADS`

### FTK_THREAD_TBB
- Requires `FTK_HAVE_TBB`
- Uses Intel TBB for work stealing
- Best for irregular/dynamic workloads
- Automatic load balancing

## Parallel Execution Patterns

### Pattern 1: Parallel Feature Detection (Thread-Safe)

```cpp
// Configure tracker BEFORE parallel execution
tracker.set_domain(domain);
tracker.set_threshold(threshold);
tracker.initialize();

// Parallel detection (if tracker supports it internally)
// Example: critical_point_tracker with TBB backend
tracker.set_thread_backend(FTK_THREAD_TBB);
tracker.detect_features(); // Internally parallel and thread-safe
```

### Pattern 2: Parallel Timestep Processing (Requires Synchronization)

```cpp
std::mutex push_mutex;

// Process timesteps in parallel
#pragma omp parallel for
for (int t = 0; t < n_timesteps; t++) {
    ndarray<double> field = load_field(t);

    // Synchronize push operations
    {
        std::lock_guard<std::mutex> lock(push_mutex);
        tracker.push_field_data_snapshot(field);
    }
}

// Finalize after all pushes complete
tracker.finalize();
```

### Pattern 3: Distributed Processing with MPI

```cpp
diy::mpi::communicator comm;
tracker.set_communicator(comm);

// Each MPI rank processes its subdomain
tracker.set_local_domain(get_rank_subdomain(comm.rank()));
tracker.initialize();

// Independent processing per rank
for (int t = 0; t < n_timesteps; t++) {
    tracker.push_field_data_snapshot(field[t]);
}

// Collective finalization
tracker.finalize(); // Includes MPI collectives
```

## MPI and Thread Interaction

When using both MPI and threading:

1. **MPI Thread Safety Level**
   - FTK assumes `MPI_THREAD_SERIALIZED` or higher
   - MPI calls are made from a single thread per rank
   - If using `MPI_THREAD_MULTIPLE`, ensure your MPI implementation supports it

2. **Collective Operations**
   - `finalize()` typically includes MPI collectives
   - All ranks must call collectively
   - Do NOT call from within a parallel region

3. **DIY Framework**
   - FTK uses DIY for distributed parallelism
   - DIY operations are thread-safe within a single rank
   - Cross-rank communication requires synchronization

## Hardware Accelerators

### CUDA Support (FTK_HAVE_CUDA)
- CUDA kernels are launched asynchronously
- Synchronization happens at tracker boundaries
- Device memory management is **NOT thread-safe**
- Only one CPU thread should manage GPU operations

### SYCL/hipSYCL Support
- Similar to CUDA thread-safety model
- Queue submissions should be serialized
- Device-side operations are independent

## Best Practices

### Do:
- ✅ Configure trackers before parallel execution
- ✅ Use tracker's built-in parallelism when available
- ✅ Read mesh data concurrently after construction
- ✅ Use TBB-based trackers for automatic parallelism
- ✅ Protect push operations with mutexes if needed
- ✅ Initialize MPI with `MPI_THREAD_SERIALIZED` or higher

### Don't:
- ❌ Modify tracker configuration during parallel execution
- ❌ Call push methods from multiple threads without synchronization
- ❌ Mix threading backends within the same tracker
- ❌ Modify meshes from multiple threads
- ❌ Call MPI collectives from parallel regions
- ❌ Share CUDA/SYCL contexts across threads

## Thread-Safety by Component

| Component | Thread-Safe? | Notes |
|-----------|--------------|-------|
| `tracking_graph` | ✅ Yes | Uses internal mutex |
| `duf` (union-find) | ✅ Yes | Uses internal mutex |
| Mesh classes | ⚠️ Read-only | Safe for queries, not modifications |
| Feature structures | ⚠️ Read-only | Safe for queries, not modifications |
| Tracker configuration | ❌ No | Configure before parallel execution |
| `push_field_data_snapshot()` | ❌ No | Requires external synchronization |
| File I/O | ❌ No | Serialize I/O operations |
| TBB containers | ✅ Yes | When built with TBB support |

## Files Using Threading Primitives

The following files contain explicit threading synchronization:

### Mutexes (std::mutex):
- `include/ftk/basic/duf.hh`
- `include/ftk/tracking_graph/tracking_graph.hh`
- `include/ftk/filters/filter.hh`
- `include/ftk/filters/critical_point_tracker.hh`
- `include/ftk/filters/critical_point_tracker_2d_regular.hh`
- `include/ftk/filters/critical_point_tracker_3d_regular.hh`
- `include/ftk/filters/critical_point_tracker_2d_unstructured.hh`
- `include/ftk/filters/critical_point_tracker_3d_unstructured.hh`
- `include/ftk/filters/contour_tracker_2d_regular.hh`
- `include/ftk/filters/contour_tracker_3d_regular.hh`
- `include/ftk/filters/critical_line_tracker_3d_regular.hh`
- `include/ftk/filters/critical_line_tracker_3d_unstructured.hh`
- `include/ftk/filters/tdgl_vortex_tracker_3d_regular.hh`
- `include/ftk/filters/xgc_blob_filament_tracker.hh`
- `include/ftk/filters/xgc_blob_threshold_tracker.hh`

### TBB Concurrent Containers:
- `include/ftk/filters/xgc_blob_filament_tracker.hh` - `tbb::concurrent_hash_map`, `tbb::concurrent_unordered_set`
- `include/ftk/mesh/simplicial_unstructured_2d_mesh.hh` - `tbb::concurrent_vector`
- `include/ftk/utils/gather.hh` - TBB parallel algorithms
- `include/ftk/utils/serialization.hh` - TBB support

## Debugging Thread Issues

### Tools:
- **ThreadSanitizer (TSan)**: Compile with `-fsanitize=thread` to detect data races
- **Helgrind**: Valgrind tool for detecting synchronization errors
- **Intel Inspector**: For TBB-based code analysis

### Common Issues:
1. **Race conditions in push operations**: Always protect with mutex
2. **Deadlocks**: Ensure consistent lock ordering across threads
3. **False sharing**: Performance degradation from cache line contention
4. **MPI thread level**: Verify with `MPI_Query_thread()`

## Future Improvements

Areas for potential thread-safety enhancements:
- Fine-grained locking in trackers (reduce contention)
- Lock-free data structures for feature collections
- Better thread-safe I/O buffering
- Thread-local storage for temporary computations
- Concurrent push operations with internal queuing

## Questions?

If you encounter thread-safety issues or have questions:
1. Check this document first
2. Review the source code for mutex usage
3. Enable ThreadSanitizer for testing
4. Report issues at: https://github.com/hguo/ftk/issues

---

**Last Updated**: February 2026
**FTK Version**: Post ndarray migration (no_ndarray branch)
