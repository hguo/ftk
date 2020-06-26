#ifndef _FTK_LEVELSET_TRACKER
#define _FTK_LEVELSET_TRACKER

#include <ftk/filters/connected_component_tracker.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>

namespace ftk {

template <typename TimeIndexType=size_t, typename LabelIdType=size_t>
struct levelset_tracker : public connected_component_tracker<TimeIndexType, LabelIdType>
{
  levelset_tracker() {}
  virtual ~levelset_tracker() {};

  void set_threshold(double threshold, bool above=true);
  void set_input_shape(const lattice& shape);

  template <typename FloatType>
  void push_scalar_field_data_snapshot(const ndarray<FloatType>&);

protected:
  double threshold = 0.0;
  bool above = true;
};

///////////////
template <typename TimeIndexType, typename LabelIdType>
void levelset_tracker<TimeIndexType, LabelIdType>::set_threshold(double t, bool a)
{
  threshold = t;
  above = a;
}

template <typename TimeIndexType, typename LabelIdType>
template <typename FloatType>
void levelset_tracker<TimeIndexType, LabelIdType>::push_scalar_field_data_snapshot(const ndarray<FloatType>& array)
{
  ndarray<LabelIdType> labels; 
  labels.reshape(array);
  for (auto i = 0; i < array.nelem(); i ++) {
    if ((above && array[i] >= threshold) || (!above && array[i] <= threshold))
      labels[i] = 1;
  }

  // relabel
  if (array.nd() == 2) hoshen_kopelman_2d(labels);
  else assert(false); // not yet implemented

  this->push_labeled_data_snapshot(labels.std_vector());
}

}

#endif
