#ifndef _FTK_FEATURE_POINT_SET_HH
#define _FTK_FEATURE_POINT_SET_HH

#include <ftk/features/feature_point.hh>
#include <list>

#if FTK_HAVE_VTK
#include <vtkPolyData.h>
#endif

namespace ftk {

struct feature_point_set_t : public std::list<feature_point_t>
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> to_vtp(
      const std::vector<std::string>& scalar_components) const;
#endif
};

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> feature_point_set_t::to_vtp(
    const std::vector<std::string>& scalar_components) const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &cp : *this) {
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  poly->SetPoints(points);
  poly->SetVerts(vertices);
  
  if (1) {
    vtkSmartPointer<vtkUnsignedIntArray> id_array = vtkSmartPointer<vtkUnsignedIntArray>::New();
    id_array->SetNumberOfValues(size());
    int i = 0;
    for (const auto &cp : *this)
      id_array->SetValue(i++, cp.id);
    id_array->SetName("id");
    poly->GetPointData()->AddArray(id_array);
  }
 
  if (1) {
    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetNumberOfValues(size());
    int i = 0;
    for (const auto &cp : *this)
      time_array->SetValue(i++, cp.t);
    time_array->SetName("time");
    poly->GetPointData()->AddArray(time_array);
  }
  
  if (1) { 
    vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
    types->SetNumberOfValues(size());
    int i = 0;
    for (const auto &cp : *this)
      types->SetValue(i++, cp.type);
    types->SetName("type");
    poly->GetPointData()->AddArray(types);
  }
  
  for (auto k = 0; k < scalar_components.size(); k ++) {
    vtkSmartPointer<vtkFloatArray> scalar = vtkSmartPointer<vtkFloatArray>::New();
    scalar->SetNumberOfValues(size());
    size_t i = 0;
    for (const auto &cp : *this)
      scalar->SetValue(i++, cp.scalar[k]);
    scalar->SetName(scalar_components[k].c_str());
    poly->GetPointData()->AddArray(scalar);
  }

  return poly;
}
#endif

}

#endif
