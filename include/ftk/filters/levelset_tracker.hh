#ifndef _FTK_LEVELSET_TRACKER
#define _FTK_LEVELSET_TRACKER

#include <ftk/filters/connected_component_tracker.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>

namespace ftk {

enum {
  FTK_COMPARE_GE, // greater or equal to
  FTK_COMPARE_GT, // greater than
  FTK_COMPARE_LE, // less or equal to
  FTK_COMPARE_LT // less than
};

template <typename TimeIndexType=size_t, typename LabelIdType=size_t>
struct levelset_tracker : public connected_component_tracker<TimeIndexType, LabelIdType>
{
  levelset_tracker() {}
  virtual ~levelset_tracker() {};

  void set_threshold(double threshold, int mode=FTK_COMPARE_GE);
  void set_input_shape(const lattice& shape);

  template <typename FloatType>
  void push_scalar_field_data_snapshot(const ndarray<FloatType>&);

  ndarray<LabelIdType> get_last_labeled_array_snapshot() const;

protected:
  double threshold = 0.0;
  int mode = FTK_COMPARE_GE;

  std::vector<size_t> input_shape;
};

///////////////
template <typename TimeIndexType, typename LabelIdType>
void levelset_tracker<TimeIndexType, LabelIdType>::set_threshold(double t, int m)
{
  threshold = t;
  mode = m;
}

template <typename TimeIndexType, typename LabelIdType>
template <typename FloatType>
void levelset_tracker<TimeIndexType, LabelIdType>::push_scalar_field_data_snapshot(const ndarray<FloatType>& array)
{
  input_shape = array.shape();

  ndarray<LabelIdType> labels; 
  labels.reshape(array);
  for (auto i = 0; i < array.nelem(); i ++) {
    switch (mode) {
    case FTK_COMPARE_GE: labels[i] = array[i] >= threshold; break;
    case FTK_COMPARE_GT: labels[i] = array[i] > threshold; break;
    case FTK_COMPARE_LE: labels[i] = array[i] <= threshold; break;
    case FTK_COMPARE_LT: labels[i] = array[i] < threshold; break;
    default: break;
    }
  }

  // relabel w/ ccl
  if (array.nd() == 2) hoshen_kopelman_2d(labels);
  else assert(false); // not yet implemented

  this->push_labeled_data_snapshot(labels.std_vector());
}

template <typename TimeIndexType, typename LabelIdType>
ndarray<LabelIdType> levelset_tracker<TimeIndexType, LabelIdType>::get_last_labeled_array_snapshot() const
{
  ndarray<LabelIdType> array(input_shape);
  array.from_vector(this->get_last_labeled_data_snapshot());
  return array;
}

}

#endif
