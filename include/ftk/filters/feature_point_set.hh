#ifndef _FTK_FEATURE_POINT_SET_HH
#define _FTK_FEATURE_POINT_SET_HH

#include <ftk/filters/feature_point_set.hh>

namespace ftk {

struct feature_point_set_t : public std::set<feature_point_t>
{

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> to_vtp(int fpdims) const;
#endif
};

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> feature_point_set::to_vtp(int fpdims) const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &cp : get_intersections()) {
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    if (fpdims == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  
  vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
  time_array->SetNumberOfValues(get_intersections().size());
  int i = 0;
  for (const auto &cp : get_intersections())
    time_array->SetValue(i++, cp.t);
  time_array->SetName("time");
  polyData->GetPointData()->AddArray(time_array);

  return polyData;

}
#endif

}

#endif
