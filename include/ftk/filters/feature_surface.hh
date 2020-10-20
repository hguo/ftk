#ifndef _FTK_FEATURE_SURFACE_HH
#define _FTK_FEATURE_SURFACE_HH

#include <ftk/filters/feature_curve.hh>

#if FTK_HAVE_VTK
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#endif

namespace ftk {

struct feature_surface_t {
  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 3>> conn; // triangles

  void reorient();

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const; // preferred
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

  vtkSmartPointer<vtkDataArray> array_grad = vtkDoubleArray::New();
  array_grad->SetNumberOfComponents(3);
  array_grad->SetNumberOfTuples(pts.size());

  for (int i = 0; i < pts.size(); i ++) {
    const auto &p = pts[i];
    points->InsertNextPoint(p.x[0], p.x[1], p.x[2]);
    array_id->SetTuple1(i, p.id);
    array_time->SetTuple1(i, p.t);
    array_grad->SetTuple3(i, p.v[0], p.v[1], p.v[2]);
  }
  poly->SetPoints(points);
  poly->GetPointData()->AddArray(array_id);
  poly->GetPointData()->AddArray(array_time);
  poly->GetPointData()->SetNormals(array_grad);

  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    vtkSmartPointer<vtkTriangle> tri = vtkTriangle::New();
    tri->GetPointIds()->SetId(0, c[0]);
    tri->GetPointIds()->SetId(1, c[1]);
    tri->GetPointIds()->SetId(2, c[2]);
    cells->InsertNextCell(tri);
  }
  poly->SetPolys(cells);

  return poly;
}

inline vtkSmartPointer<vtkUnstructuredGrid> feature_surface_t::to_vtu() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();

  vtkSmartPointer<vtkDataArray> array_id = vtkUnsignedIntArray::New();
  array_id->SetName("id");
  array_id->SetNumberOfComponents(1);
  array_id->SetNumberOfTuples(pts.size());

  vtkSmartPointer<vtkDataArray> array_time = vtkDoubleArray::New();
  array_time->SetName("time");
  array_time->SetNumberOfComponents(1);
  array_time->SetNumberOfTuples(pts.size());
  
  vtkSmartPointer<vtkDataArray> array_scalar = vtkDoubleArray::New();
  array_scalar->SetName("scalar");
  array_scalar->SetNumberOfComponents(1);
  array_scalar->SetNumberOfTuples(pts.size());
  
  vtkSmartPointer<vtkDataArray> array_grad = vtkDoubleArray::New();
  array_grad->SetNumberOfComponents(3);
  array_grad->SetNumberOfTuples(pts.size());

  for (int i = 0; i < pts.size(); i ++) {
    const auto &p = pts[i];
    points->InsertNextPoint(p.x[0], p.x[1], p.x[2]);
    array_id->SetTuple1(i, p.id);
    array_time->SetTuple1(i, p.t);
    array_scalar->SetTuple1(i, p.scalar[0]);
    array_grad->SetTuple3(i, p.v[0], p.v[1], p.v[2]);
  }
  grid->SetPoints(points);
  grid->GetPointData()->AddArray(array_id);
  grid->GetPointData()->AddArray(array_time);
  grid->GetPointData()->AddArray(array_scalar);
  grid->GetPointData()->SetNormals(array_grad);

  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    vtkIdType ids[3] = {c[0], c[1], c[2]};
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  }

  return grid;
}
#endif

inline void feature_surface_t::reorient()
{
  fprintf(stderr, "reorienting, #pts=%zu, #tris=%zu\n", pts.size(), conn.size());
  auto edge = [](int i, int j) {
    if (i > j) std::swap(i, j);
    return std::make_tuple(i, j);
  };

  // 1. build triangle-triangle graph
  std::map<std::tuple<int, int>, std::set<int>> edge_triangle;
  for (int i = 0; i < conn.size(); i ++) {
    auto tri = conn[i];
    edge_triangle[edge(tri[0], tri[1])].insert(i);
    edge_triangle[edge(tri[1], tri[2])].insert(i);
    edge_triangle[edge(tri[0], tri[2])].insert(i);
  }

  std::vector<std::set<int>> neighbors(conn.size());
  auto add_neighbors = [&neighbors](int i, const std::set<int>& set) {
    for (const auto j : set)
      if (i != j) 
        neighbors[i].insert(j);
  };
 
  for (int i = 0; i < conn.size(); i ++) {
    auto tri = conn[i];
    add_neighbors(i, edge_triangle[edge(tri[0], tri[1])]);
    add_neighbors(i, edge_triangle[edge(tri[1], tri[2])]);
    add_neighbors(i, edge_triangle[edge(tri[0], tri[2])]);
  }

  // for (int i = 0; i < neighbors.size(); i ++)
    // fprintf(stderr, "%d, %zu\n", i, neighbors[i].size());

  // 2. reorientate triangles via bfs
}

}

#endif
