#ifndef _FTK_FEATURE_SURFACE_HH
#define _FTK_FEATURE_SURFACE_HH

#include <ftk/filters/feature_curve.hh>

#if FTK_HAVE_VTK
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkXMLPolyDataWriter.h>
#endif

namespace ftk {

struct feature_surface_t {
  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 3>> conn;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> to_vtp() const;
#endif
};

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> feature_surface_t::to_vtp() const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  vtkSmartPointer<vtkDataArray> array_id = vtkUnsignedIntArray::New();
  array_id->SetName("id");
  array_id->SetNumberOfComponents(1);
  array_id->SetNumberOfTuples(pts.size());

  vtkSmartPointer<vtkDataArray> array_time = vtkDoubleArray::New();
  array_time->SetName("time");
  array_time->SetNumberOfComponents(1);
  array_time->SetNumberOfTuples(pts.size());

  for (int i = 0; i < pts.size(); i ++) {
    const auto &p = pts[i];
    points->InsertNextPoint(p.x[0], p.x[1], p.x[2]);
    array_id->SetTuple1(i, p.id);
    array_time->SetTuple1(i, p.t);
  }
  poly->SetPoints(points);
  poly->GetPointData()->AddArray(array_id);
  poly->GetPointData()->AddArray(array_time);

  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    vtkSmartPointer<vtkTriangle> tri = vtkTriangle::New();
    tri->GetPointIds()->SetId(0, std::get<0>(c));
    tri->GetPointIds()->SetId(1, std::get<1>(c));
    tri->GetPointIds()->SetId(2, std::get<2>(c));
    cells->InsertNextCell(tri);
  }
  poly->SetPolys(cells);

  return poly;
}
#endif

}

#endif
