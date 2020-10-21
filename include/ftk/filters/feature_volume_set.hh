#ifndef _FEATURE_VOLUME_SET_H
#define _FEATURE_VOLUME_SET_H

#include <ftk/filters/feature_volume.hh>

#if FTK_HAVE_VTK
#include <vtkAppendFilter.h>
#endif

namespace ftk {

struct feature_volume_set_t : public std::map<int, feature_volume_t>
{
  int add(const feature_volume_t& f);
  std::vector<int> add(const std::vector<feature_volume_t>&);

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
#endif
};

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkUnstructuredGrid> feature_volume_set_t::to_vtu() const
{
  vtkSmartPointer<vtkAppendFilter> append_filter = vtkAppendFilter::New();
  for (const auto &kv : *this)
    append_filter->AddInputData(kv.second.to_vtu());
  append_filter->Update();

  vtkSmartPointer<vtkUnstructuredGrid> combined = append_filter->GetOutput();
  return combined;
}
#endif

}

#endif
