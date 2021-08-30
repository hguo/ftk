#ifndef _FTK_FEATURE_SURFACE_HH
#define _FTK_FEATURE_SURFACE_HH

#include <ftk/config.hh>
#include <ftk/features/feature_curve.hh>
#include <ftk/features/feature_curve_set.hh>
#include <ftk/basic/simple_union_find.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/utils/serialization.hh>
#include <ftk/io/util.hh>

#if FTK_HAVE_VTK
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#endif

#if FTK_HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#endif

namespace ftk {

struct feature_surface_t {
  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 3>> tris; 
  std::vector<std::array<int, 4>> quads;
  std::vector<std::array<int, 5>> pentagons;

  void clear();
  void relabel();

  feature_curve_set_t slice_time(int t) const;

  void triangulate(); // WIP
  void reorient(); // WIP

  void discard(std::function<bool(const feature_point_t&)> f);

public: // IO
  void save(const std::string& filename, std::string format="") const;
  void load(const std::string& filename, std::string format="");

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
  vtkSmartPointer<vtkPolyData> to_vtp(bool generate_normal=false, bool spacetime=false, double z_time_scale=1.0) const;
#endif
};

} // namespace ftk

//////////////////////// serialization
namespace diy {
  template <> struct Serialization<ftk::feature_surface_t> {
    static void save(diy::BinaryBuffer& bb, const ftk::feature_surface_t& s) {
      diy::save(bb, s.pts); 
      diy::save(bb, s.tris); 
      diy::save(bb, s.quads); 
      diy::save(bb, s.pentagons); 
    }

    static void load(diy::BinaryBuffer& bb, ftk::feature_surface_t &s) {
      diy::load(bb, s.pts); 
      diy::load(bb, s.tris); 
      diy::load(bb, s.quads); 
      diy::load(bb, s.pentagons); 
    }
  };
}
////////////////////////

namespace ftk {

inline void feature_surface_t::clear()
{
  pts.clear();
  tris.clear();
  quads.clear();
  pentagons.clear();
}

inline void feature_surface_t::relabel()
{
  simple_union_find<int> uf(pts.size());
  for (auto tet : tris)
    for (int j = 1; j < 3; j ++)
      uf.unite(tet[0], tet[j]);

  for (auto quad : quads) 
    for (int j = 1; j < 4; j ++)
      uf.unite(quad[0], quad[j]);

  for (auto pentagon : pentagons)
    for (int j = 1; j < 4; j ++)
      uf.unite(pentagon[0], pentagon[j]);

  for (int i = 0; i < pts.size(); i ++)
    pts[i].id = uf.find(i);
}

inline void feature_surface_t::triangulate()
{
#if FTK_HAVE_CGAL
  typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
  typedef K::FT                                                   FT;
  typedef K::Point_3                                              Point_3;
  typedef CGAL::Surface_mesh<Point_3>                             Mesh;
  typedef std::array<FT, 3>                                       Custom_point;
  typedef std::vector<std::size_t>                                CGAL_Polygon;

  namespace PMP = CGAL::Polygon_mesh_processing;
#endif
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> feature_surface_t::to_vtp(bool generate_normal, bool spacetime, double z_time_scale) const
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

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
  
  vtkSmartPointer<vtkDataArray> array_type = vtkUnsignedIntArray::New();
  array_type->SetName("type");
  array_type->SetNumberOfComponents(1);
  array_type->SetNumberOfTuples(pts.size());

  vtkSmartPointer<vtkDataArray> array_grad = vtkFloatArray::New();
  array_grad->SetName("grad");
  array_grad->SetNumberOfComponents(3);
  array_grad->SetNumberOfTuples(pts.size());

  for (int i = 0; i < pts.size(); i ++) {
    const auto &p = pts[i];
    if (spacetime) 
      points->InsertNextPoint(p.x[0], p.x[1], p.t * z_time_scale);
    else 
      points->InsertNextPoint(p.x[0], p.x[1], p.x[2]);

    array_id->SetTuple1(i, p.id);
    array_time->SetTuple1(i, p.t);
    array_scalar->SetTuple1(i, p.scalar[0]);
    array_type->SetTuple1(i, p.type);
    array_grad->SetTuple3(i, p.v[0], p.v[1], p.v[2]);
  }
  poly->SetPoints(points);
  poly->GetPointData()->AddArray(array_id);
  poly->GetPointData()->AddArray(array_time);
  poly->GetPointData()->AddArray(array_scalar);
  poly->GetPointData()->AddArray(array_type);
  poly->GetPointData()->SetNormals(array_grad);

  for (int i = 0; i < tris.size(); i ++) {
    const auto &c = tris[i];
    vtkSmartPointer<vtkTriangle> tri = vtkTriangle::New();
    for (int j = 0; j < 3; j ++)
      tri->GetPointIds()->SetId(j, c[j]);
    cells->InsertNextCell(tri);
  }
  for (int i = 0; i < quads.size(); i ++) {
    const auto &q = quads[i];
    vtkSmartPointer<vtkQuad> quad = vtkQuad::New();
    for (int j = 0; j < 4; j ++)
      quad->GetPointIds()->SetId(j, q[j]);
    cells->InsertNextCell(quad);
  }
  poly->SetPolys(cells);

  if (generate_normal) {
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(poly);
    normalGenerator->ConsistencyOn();
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    // normalGenerator->SetFlipNormals(true);
    normalGenerator->AutoOrientNormalsOn();
    normalGenerator->Update();
    return normalGenerator->GetOutput();
  } else 
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

  for (int i = 0; i < tris.size(); i ++) {
    const auto &c = tris[i];
    vtkIdType ids[3] = {c[0], c[1], c[2]};
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  }

  return grid;
}
#endif

inline void feature_surface_t::reorient()
{
  fprintf(stderr, "reorienting, #pts=%zu, #tris=%zu\n", pts.size(), tris.size());
  auto edge = [](int i, int j) {
    if (i > j) std::swap(i, j);
    return std::make_tuple(i, j);
  };

  // 1. build triangle-triangle graph
  std::map<std::tuple<int, int>, std::set<int>> edge_triangle;
  for (int i = 0; i < tris.size(); i ++) {
    auto tri = tris[i];
    for (int j = 0; j < 3; j ++)
      edge_triangle[edge(tri[j], tri[(j+1)%3])].insert(i);
  }

  auto chirality = [](std::array<int, 3> a, std::array<int, 3> b) {
    std::vector<std::tuple<int, int>> ea, eb;
    for (int i = 0; i < 3; i ++) {
      ea.push_back(std::make_tuple(a[i], a[(i+1)%3]));
      eb.push_back(std::make_tuple(b[i], b[(i+1)%3]));
    }

    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++) {
        if (ea[i] == eb[i]) return 1;
        else if (ea[i] == std::make_tuple(std::get<1>(eb[j]), std::get<0>(eb[j]))) 
          return -1;
      }

    assert(false);
    return 0;
  };

  // 2. reorientate triangles with bfs
  std::set<int> visited;
  std::queue<int> Q;
  Q.push(0);
  visited.insert(0);

  while (!Q.empty()) {
    auto current = Q.front();
    Q.pop();

    const int i = current;

    for (int j = 0; j < 3; j ++) {
      auto e = edge(tris[i][j], tris[i][(j+1)%3]);
      auto neighbors = edge_triangle[e];
      neighbors.erase(i);

      // fprintf(stderr, "#neighbors=%zu\n", neighbors.size());
      for (auto k : neighbors)
        if (visited.find(k) == visited.end()) {
          // fprintf(stderr, "pushing %d, chi=%d\n", k, chirality(tris[i], tris[k]));
          if (chirality(tris[i], tris[k]) < 0) {
            // fprintf(stderr, "flipping %d\n", k);
            std::swap(tris[k][0], tris[k][1]);
          }
          Q.push(k); // std::make_tuple(k, chirality(tris[i], tris[k])));
          visited.insert(k);
        }
    }
  }
}

inline feature_curve_set_t feature_surface_t::slice_time(int t) const
{
  feature_curve_set_t curve_set;
  
  std::vector<feature_point_t> mypts;
  std::set<int> nodes;
  std::map<int, std::set<int>> links;

  std::map<int, int> map; // id, new_id
  int j = 0;
  for (int i = 0; i < pts.size(); i ++)
    if (pts[i].timestep == t && pts[i].ordinal) {
      nodes.insert(j);
      map[i] = j;
      mypts.push_back(pts[i]);
      j ++;
    }

  fprintf(stderr, "#pts=%zu\n", map.size());

  for (int i = 0; i < tris.size(); i ++) {
    const auto &c = tris[i];

    std::set<int> cc;
    for (int j = 0; j < 3; j ++)
      if (map.find(c[j]) != map.end())
        cc.insert(map[c[j]]);

    if (!cc.empty()) {
      for (const auto k : cc)
        links[k].insert(cc.begin(), cc.end());
    }
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
        traj.push_back(mypts[linear_graphs[j][k]]);
      // curve_set.add(traj); // , traj[0].id);
      curve_set.add(traj, traj[0].id);
      // fprintf(stderr, "curve_set.size=%zu\n", curve_set.size());
    }
  }

  // fprintf(stderr, "curve_set.size=%zu\n", curve_set.size());
  size_t count = 0;
  for (const auto &kv : curve_set)
    count += kv.second.size();
  fprintf(stderr, "pts_count=%zu, curve_count=%zu\n", count, curve_set.size());

  return curve_set;
}

inline void feature_surface_t::load(const std::string& filename, std::string format)
{
  const int fmt = file_extension(filename, format);
  if (fmt == FILE_EXT_BIN) {
    diy::unserializeFromFile(filename, *this);
  } else {
    fprintf(stderr, "unsupported file format\n");
    assert(false);
  }
}

inline void feature_surface_t::save(const std::string& filename, std::string format) const
{
  const int fmt = file_extension(filename, format);
  if (fmt == FILE_EXT_BIN || fmt == FILE_EXT_NULL) {
    diy::serializeToFile(*this, filename);
  } else {
#if FTK_HAVE_VTK
    fprintf(stderr, "writing polydata, #pts=%zu, #tris=%zu\n", pts.size(), tris.size());
    write_polydata(filename, to_vtp(), format);
#else
    fprintf(stderr, "unsupported file format\n");
    assert(false);
    // fatal("unsupported file format");
#endif
  }
}

inline void feature_surface_t::discard(std::function<bool(const feature_point_t&)> f)
{
  std::vector<int> parent(pts.size());
  for (int i = 0; i < pts.size(); i ++) 
    parent[i] = i;
  std::vector<std::set<int>> children(pts.size());

  // build node-edge-node links
  std::vector<std::set<int>> neighbors(pts.size());
  for (int i = 0; i < tris.size(); i ++) {
    const auto tri = tris[i];

    for (int j = 0; j < 3; j ++)
      for (int k = 0; k < 3; k ++)
        if (j != k)
          neighbors[tri[j]].insert(tri[k]);
  }

  // def find_nearest(i):
  //   find nearest j among neighbors of i
  //   while parent[j] != j:   // j is deleted
  //     j <-- parent[j]
  auto nearest = [&](int i) {
    double mindist = std::numeric_limits<double>::max();
    int mindistj = -1;
    for (auto j : neighbors[i]) {
      double d = dist(pts[i], pts[j]);
      if (d < mindist) {
        mindist = d;
        mindistj = j;
      }
    }

    int j = mindistj;
    while (j != -1 && parent[j] != j) 
      j = parent[j];

    return j;
  };

  // foreach i : verts
  //   if f(i):
  //     if j <-- find_nearest(i):
  //       parent[i] <-- j
  //       children[j].add(i)
  //       children[j].unite(children[i])
  //       foreach k : children[i]: 
  //         parent[k] <-- j
  //       children[i] <-- empty
  //       mark i as deleted, i.e. parent[i] != i
  //     else: 
  //       mark i as permenantly deleted, or parent[i] = -1;
  for (int i = 0; i < pts.size(); i ++) {
    if (f(pts[i])) {
      const int j = nearest(i);
      if (j != -1) {
        children[j].insert(i);
        children[j].insert( children[i].begin(), children[i].end() );
        for (auto k : children[i])
          parent[k] = j;
        children[i].clear();
        parent[i] = j;
      } else 
        parent[i] = -1;
    }
  }

  // build new pts and tris
  std::map<int, int> ptsmap; // old id to new id
  std::vector<feature_point_t> newpts;
  std::vector<std::array<int, 3>> newtris;
  for (int i = 0; i < pts.size(); i ++) {
    if (parent[i] == i) {
      ptsmap[i] = newpts.size();
      newpts.push_back(pts[i]);
    }
  }

  for (int i = 0; i < tris.size(); i ++) {
    const auto tri = tris[i];
    bool b[3];
    for (int k = 0; k < 3; k ++)
      b[k] = ptsmap.find(tri[k]) != ptsmap.end();

    if (b[0] && b[1] && b[2]) {
      newtris.push_back( {ptsmap[tri[0]], ptsmap[tri[1]], ptsmap[tri[2]]} );
    }
  }

  pts = newpts;
  tris = newtris;
}

}

#endif
