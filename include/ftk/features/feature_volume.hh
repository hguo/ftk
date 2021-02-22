#ifndef _FTK_FEATURE_VOLUME_HH
#define _FTK_FEATURE_VOLUME_HH

#include <ftk/features/feature_surface.hh>
#include <ftk/basic/simple_union_find.hh>

#if FTK_HAVE_VTK
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

namespace ftk {

struct feature_volume_t {
  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 4>> conn;

  void clear() { pts.clear(); conn.clear(); }
  void relabel();

  feature_surface_t slice(std::function<bool(const feature_point_t&)>) const;
  feature_surface_t slice_time(int t) const;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
#endif
};

inline void feature_volume_t::relabel()
{
  simple_union_find<int> uf(pts.size());
  for (auto tet : conn) {
    uf.unite(tet[0], tet[1]);
    uf.unite(tet[0], tet[2]);
    uf.unite(tet[0], tet[3]);
  }

  for (int i = 0; i < pts.size(); i ++)
    pts[i].id = uf.find(i);
}

inline feature_surface_t feature_volume_t::slice(std::function<bool(const feature_point_t&)> f) const 
{
  feature_surface_t surf;

  std::map<int, int> map; // id, new_id
  int j = 0;
  for (int i = 0; i < pts.size(); i ++)
    if (f(pts[i])) {
      map[i] = j ++;
      surf.pts.push_back(pts[i]);
    }

  std::set<std::array<int, 3>> tris;
  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    int count = 0;
    std::array<int, 4> q;
    
    for (int j = 0; j < 4; j ++)
      if (map.find(c[j]) != map.end())
        q[count ++] = map[c[j]];

    if (count == 3) {
      std::array<int, 3> q3 = {q[0], q[1], q[2]};
      std::sort(q3.begin(), q3.end());
      tris.insert(q3);
      // surf.conn.push_back({q[0], q[1], q[2]});
    }
  }

  for (auto tri : tris)
    surf.tris.push_back(tri);

  return surf;
}

inline feature_surface_t feature_volume_t::slice_time(int t) const
{
  feature_surface_t surf;

  std::map<int, int> map; // id, new_id
  int j = 0;
  for (int i = 0; i < pts.size(); i ++)
    if (pts[i].timestep == t && pts[i].ordinal) {
      map[i] = j ++;
      surf.pts.push_back(pts[i]);
    }

  std::set<std::array<int, 3>> tris;
  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    int count = 0;
    std::array<int, 4> q;
 
    bool valid = true;
    for (int j = 0; j < 4; j ++)
      if (pts[c[j]].timestep != t)
        valid = false;

    if (valid) {
      for (int j = 0; j < 4; j ++)
        if (map.find(c[j]) != map.end()) {
          q[count ++] = map[c[j]];
        }

      if (count == 3) {
        std::array<int, 3> q3 = {q[0], q[1], q[2]};
        std::sort(q3.begin(), q3.end());
        tris.insert(q3);
        // surf.conn.push_back({q[0], q[1], q[2]});
      }
    }
  }

  for (auto tri : tris)
    surf.tris.push_back(tri);

  fprintf(stderr, "sliced time=%d, #pts=%zu, #tri=%zu\n", t, surf.pts.size(), surf.tris.size());
  return surf;
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkUnstructuredGrid> feature_volume_t::to_vtu() const
{
  // fprintf(stderr, "%zu, %zu\n", pts.size(), conn.size());
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

  for (int i = 0; i < pts.size(); i ++) {
    const auto &p = pts[i];
    points->InsertNextPoint(p.x[0], p.x[1], p.x[2]);
    array_id->SetTuple1(i, p.id);
    array_time->SetTuple1(i, p.t);
  }
  grid->SetPoints(points);
  grid->GetPointData()->AddArray(array_id);
  grid->GetPointData()->AddArray(array_time);

  for (int i = 0; i < conn.size(); i ++) {
    const auto &c = conn[i];
    vtkIdType ids[4] = {c[0], c[1], c[2], c[3]};
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  }

  // grid->PrintSelf(std::cerr, vtkIndent(2));
  return grid;
}
#endif

}

#endif
