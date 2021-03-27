#ifndef _FTK_FEATURE_LINE_HH
#define _FTK_FEATURE_LINE_HH

#include <ftk/config.hh>
#include <ftk/features/feature_curve.hh>
#include <ftk/features/feature_curve_set.hh>

#if FTK_HAVE_VTK
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#endif

namespace ftk {

struct feature_line_t {
  // feature_line_t is different from feature_curve_t and feature_curve_set_t
  //  - feature_line_t consists of vertices and their connectivities
  //  - feature_curve_t represents a single curve; vertices are linearly
  //  - feature_curve_set_t is a set of feature curves
  // a feature_line_t can be converted to feature_curve_set_t; in principle, vice versa

  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 2>> edges;

  void clear();
  void relabel();

  feature_curve_set_t to_curve_set() const;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const; // to be removed
#endif
};

/////
void feature_line_t::clear()
{
  pts.clear();
  edges.clear();
}

void feature_line_t::relabel()
{
  simple_union_find<int> uf(pts.size());
  for (auto edge : edges)
    uf.unite(edge[0], edge[1]);

  for (int i = 0; i < pts.size(); i ++)
    pts[i].id = uf.find(i);
}

feature_curve_set_t feature_line_t::to_curve_set() const
{
  feature_curve_set_t curve_set;
  std::set<int> nodes;
  std::map<int, std::set<int>> links;

  for (int i = 0; i < pts.size(); i ++)
    nodes.insert(i);

  for (int i = 0; i < edges.size(); i ++) {
    const int v0 = edges[i][0], v1 = edges[i][1];
    links[v0].insert(v1);
    links[v1].insert(v0);
  }
  
  auto cc = extract_connected_components<int, std::set<int>>(
      [&](int i) { return links[i]; }, 
      nodes);
  
  for (const auto &component : cc) {
    auto linear_graphs = connected_component_to_linear_components<int>(component, [&](int i) {return links[i];});
    for (int j = 0; j < linear_graphs.size(); j ++) {
      feature_curve_t traj;
      traj.loop = is_loop<int>(linear_graphs[j], [&](int i) {return links[i];});
      
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(pts[linear_graphs[j][k]]);
      // curve_set.add(traj); // , traj[0].id);
      curve_set.add(traj, traj[0].id);
      // fprintf(stderr, "curve_set.size=%zu\n", curve_set.size());
    }
  }

  return curve_set;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkUnstructuredGrid> feature_line_t::to_vtu() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();

  vtkSmartPointer<vtkDataArray> array_id = vtkUnsignedIntArray::New();
  array_id->SetName("id");
  array_id->SetNumberOfComponents(1);
  array_id->SetNumberOfTuples(pts.size());

  vtkSmartPointer<vtkDataArray> array_time = vtkFloatArray::New();
  array_time->SetName("time");
  array_time->SetNumberOfComponents(1);
  array_time->SetNumberOfTuples(pts.size());
  
  vtkSmartPointer<vtkDataArray> array_scalar = vtkFloatArray::New();
  array_scalar->SetName("scalar");
  array_scalar->SetNumberOfComponents(1);
  array_scalar->SetNumberOfTuples(pts.size());
  
  vtkSmartPointer<vtkDataArray> array_grad = vtkFloatArray::New();
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

  for (int i = 0; i < edges.size(); i ++) {
    vtkIdType ids[2] = {edges[i][0], edges[i][1]};
    grid->InsertNextCell(VTK_LINE, 2, ids);
  }

  return grid;
}
#endif

}

#endif
