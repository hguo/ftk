#ifndef _FTK_CRITICAL_POINT_TRAJ_SET_HH
#define _FTK_CRITICAL_POINT_TRAJ_SET_HH

#include <ftk/features/feature_curve.hh>
#include <map>

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkFloatArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#endif

namespace ftk {

struct feature_curve_set_t : public std::multimap<int, feature_curve_t>
{
  int add(const feature_curve_t&);
  void add(const feature_curve_t&, int label);
  std::vector<int> add(const std::vector<feature_curve_t>&);

  std::list<feature_curve_t> to_list() const;
  void from_list(const std::list<feature_curve_t>&);

  // std::vector<int> split(int);
  void split_all();
  
  void foreach(std::function<void(const feature_curve_t&)> f) const {for (const auto& kv : *this) f(kv.second);}
  void foreach(std::function<void(feature_curve_t&)> f) {for (auto& kv : *this) f(kv.second);}
  void foreach(std::function<void(int, const feature_curve_t&)> f) const {for (const auto& kv : *this) f(kv.first, kv.second);}
  void foreach(std::function<void(int, feature_curve_t&)> f) {for (auto& kv : *this) f(kv.first, kv.second);}

  void filter(std::function<bool(const feature_curve_t&)> f);

  std::vector<feature_point_t> slice(int t) const;
  feature_curve_set_t intercept(int t0, int t1) const;

public: // IO
  void write(const std::string& format, const std::string& filename) const;
  void read(const std::string& format, const std::string& filename);

  void write_json(const std::string& filename, int indent=0) const;
  void read_json(const std::string& filename);

  void write_binary(const std::string& filename) const;
  void read_binary(const std::string& filename);

  void write_text(std::ostream& os, const int cpdims, const std::vector<std::string>& scalar_components) const;
  void write_text(const std::string& filename, const int cpdims, const std::vector<std::string>& scalar_components) const;

  void write_vtk(const std::string& filename) const;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> to_vtp(const int ncdims, 
      const std::vector<std::string> &scalar_components, 
      double tfactor=1.0) const;
#endif

protected:
  int get_new_id() const {
    if (empty()) return 0; 
    else return rbegin()->first + 1;
  }
};

}

namespace nlohmann
{
  using namespace ftk;
  template <>
  struct adl_serializer<feature_curve_set_t> {
    static void to_json(json& j, const feature_curve_set_t& s) {
      j = {{"trajs", s.to_list()}};
    }
   
    // TODO FIXME: json i/o w/ multimap has problem...
    static void from_json(const json&j, feature_curve_set_t& s) {
      std::list<feature_curve_t> list = j["trajs"];
      s.from_list(list);
    }
  };
}


// serialization
namespace diy {
  static void save(diy::BinaryBuffer& bb, const ftk::feature_curve_set_t &s) {
    diy::save(bb, s.size());
    for (const auto &kv : s) {
      diy::save(bb, kv.first);  // identifer
      diy::save(bb, kv.second); // traj
    }
  }

  static void load(diy::BinaryBuffer& bb, ftk::feature_curve_set_t &s) {
    size_t size;
    diy::load(bb, size);
    for (auto i = 0; i < size; i ++) {
      int id;
      diy::load(bb, id);
      
      ftk::feature_curve_t traj;
      diy::load(bb, traj);
      traj.relabel(id);

      s.insert({id, traj});
      // s[id] = traj;
    }
  }
} // namespace diy

//////
namespace ftk {

inline std::list<feature_curve_t> feature_curve_set_t::to_list() const
{
  std::list<feature_curve_t> list;
  for (const auto &kv : *this)
    list.push_back(kv.second);
  return list;
}

inline void feature_curve_set_t::from_list(const std::list<feature_curve_t>& list)
{
  for (const auto &c : list)
    this->insert({get_new_id(), c});
}

inline void feature_curve_set_t::write(const std::string& format, const std::string& filename) const
{
  // TODO
}

inline void feature_curve_set_t::read(const std::string& format, const std::string& filename)
{
  // TODO
}

inline void feature_curve_set_t::write_json(const std::string& filename, int indent) const
{
  // TODO
}

inline void feature_curve_set_t::read_json(const std::string& filename)
{
  // TODO
}

inline void feature_curve_set_t::write_binary(const std::string& filename) const
{
  // TODO
}

inline void feature_curve_set_t::read_binary(const std::string& filename)
{
  // TODO
}


inline void feature_curve_set_t::write_text(std::ostream& os, const int cpdims, const std::vector<std::string>& scalar_components) const
{
  os << "#trajectories=" << size() << std::endl;
  for (const auto &kv : *this) {
    const auto i = kv.first;
    const auto &curve = kv.second;
    
    os << "--trajectory " << i << ", ";
   
    if (scalar_components.size() > 0) {
      os << "min=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.min[k] << ", ";
        else os << curve.min[k] << "), ";
      
      os << "max=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.max[k] << ", ";
        else os << curve.max[k] << "), ";
      
      os << "persistence=(";
      for (int k = 0; k < scalar_components.size(); k ++)
        if (k < scalar_components.size()-1) os << curve.persistence[k] << ", ";
        else os << curve.persistence[k] << "), ";
    }

    os << "bbmin=(";
    for (int k = 0; k < cpdims; k ++)
      if (k < cpdims) os << curve.bbmin[k] << ", ";
      else os << curve.bbmin[k] << "), ";
    
    os << "bbmax=(";
    for (int k = 0; k < cpdims; k ++)
      if (k < cpdims) os << curve.bbmax[k] << ", ";
      else os << curve.bbmax[k] << "), ";
    
    os << "tmin=" << curve.tmin << ", tmax=" << curve.tmax << ", ";

    // os << "consistent_type=" << critical_point_type_to_string(cpdims, curve.consistent_type, scalar_components.size()) << ", ";
    os << "consistent_type=" << curve.consistent_type << ", ";
    os << "loop=" << curve.loop;
    os << std::endl;

    for (int k = 0; k < curve.size(); k ++) {
      os << "---";
      curve[k].print(os, cpdims, scalar_components) << std::endl;
    }
  }
}

inline void feature_curve_set_t::write_text(const std::string& filename, const int cpdims, const std::vector<std::string>& scalar_components) const
{
  // TODO
}


inline void feature_curve_set_t::write_vtk(const std::string& filename) const
{
  // TODO
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> feature_curve_set_t::to_vtp(const int cpdims, const std::vector<std::string> &scalar_components, double tfactor) const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> lines = vtkCellArray::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();

  auto pt2coords = [tfactor, cpdims](const feature_point_t& p) {
    return std::array<double, 3>{p.x[0], p.x[1], cpdims == 2 ? p.t * tfactor : p.x[2]};
  };

  foreach([&](const feature_curve_t& curve) {
    if (curve.empty()) return;
    for (auto i = 0; i < curve.size(); i ++) {
      const auto p = pt2coords(curve[i]);
      points->InsertNextPoint(p[0], p[1], p[2]);
    }
    if (curve.loop) {
      const auto p = pt2coords(curve[0]);
      points->InsertNextPoint(p[0], p[1], p[2]);
    }
  });

  size_t nv = 0;
  foreach([&](const feature_curve_t& curve) {
    if (curve.empty()) return;
    const auto npts = curve.loop ? curve.size()+1 : curve.size();
    
    vtkSmartPointer<vtkPolyLine> obj = vtkPolyLine::New();
    obj->GetPointIds()->SetNumberOfIds(npts);
    for (int i = 0; i < npts; i ++)
      obj->GetPointIds()->SetId(i, i+nv);
      
    if (curve.size() < 2) // isolated vertex
      verts->InsertNextCell(obj);
    else 
      lines->InsertNextCell(obj);
    
    nv += npts;
  });
 
  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->SetVerts(verts);

  // point data for types
  if (1) { // if (type_filter) {
    vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
    types->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        types->SetValue(i ++, curve[j].type);
      if (curve.loop)
        types->SetValue(i ++, curve[0].type);
    });
    types->SetName("type");
    polyData->GetPointData()->AddArray(types);
  }

  if (1) { // ids
    vtkSmartPointer<vtkUnsignedIntArray> ids = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ids->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](int id, const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        ids->SetValue(i ++, curve[j].id);
      if (curve.loop)
        ids->SetValue(i ++, curve[0].id);
    });
    ids->SetName("id");
    polyData->GetPointData()->AddArray(ids);
  }
  
  if (1) { // velocities
    vtkSmartPointer<vtkFloatArray> vels = vtkSmartPointer<vtkFloatArray>::New();
    vels->SetNumberOfComponents(3);
    vels->SetNumberOfTuples(nv);
    size_t i = 0;
    foreach([&](const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        vels->SetTuple3(i ++, curve[j].v[0], curve[j].v[1], curve[j].v[2]);
      if (curve.loop)
        vels->SetTuple3(i ++, curve[0].v[0], curve[0].v[1], curve[0].v[2]);
    });
    vels->SetName("velocity");
    polyData->GetPointData()->AddArray(vels);
  }
  
  if (1) { // time
    vtkSmartPointer<vtkFloatArray> time = vtkSmartPointer<vtkFloatArray>::New();
    time->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        time->SetValue(i ++, curve[j].t);
      if (curve.loop)
        time->SetValue(i ++, curve[0].t);
    });
    time->SetName("time");
    polyData->GetPointData()->AddArray(time);
  }
  
  if (1) { // condition numbers
    vtkSmartPointer<vtkFloatArray> conds = vtkSmartPointer<vtkFloatArray>::New();
    conds->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        conds->SetValue(i ++, curve[j].cond);
      if (curve.loop)
        conds->SetValue(i ++, curve[0].cond);
    });
    conds->SetName("cond");
    polyData->GetPointData()->AddArray(conds);
  }
  
  // point data for scalars
  // if (has_scalar_field) {
  for (auto k = 0; k < scalar_components.size(); k ++) {
    vtkSmartPointer<vtkFloatArray> scalar = vtkSmartPointer<vtkFloatArray>::New();
    scalar->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](const feature_curve_t& curve) {
      if (curve.empty()) return;
      for (auto j = 0; j < curve.size(); j ++)
        scalar->SetValue(i ++, curve[j].scalar[k]);
      if (curve.loop)
        scalar->SetValue(i ++, curve[0].scalar[k]);
    });
    scalar->SetName(scalar_components[k].c_str());
    polyData->GetPointData()->AddArray(scalar);
  }

  return polyData;
}
#endif

inline void feature_curve_set_t::filter(std::function<bool(const feature_curve_t&)> f)
{
  std::vector<std::multimap<int, feature_curve_t>::iterator> to_erase;

  for (auto it = this->begin(); it != this->end(); it ++)
    if (!f(it->second)) to_erase.push_back(it);
  
  for (const auto k : to_erase)
    erase(k);
}

inline int feature_curve_set_t::add(const feature_curve_t& t)
{
  const int id = get_new_id();
  auto it = insert(std::pair<int, feature_curve_t>(id, t));
  it->second.relabel(id);
  // at(id).relabel(id);
  return id;
}

inline void feature_curve_set_t::add(const feature_curve_t& t, int label)
{
  // fprintf(stderr, "inserting new curve, id=%d\n", label);
  auto it = insert(std::pair<int, feature_curve_t>(label, t));
  it->second.relabel(label);
}

inline std::vector<int> feature_curve_set_t::add(const std::vector<feature_curve_t>& trajs)
{
  std::vector<int> ids;
  for (const auto &traj : trajs)
    ids.push_back(add(traj));
  return ids;
}

#if 0
inline std::vector<int> feature_curve_set_t::split(int i)
{
  std::vector<int> result;
  const auto subtrajs = at(i).split();
  for (const auto &traj : subtrajs)
    result.push_back( add(traj) );
  erase(i);
  return result;
}
#endif

inline void feature_curve_set_t::split_all()
{
  std::vector<feature_curve_t> result;
  std::vector<std::multimap<int, feature_curve_t>::iterator> to_erase;

  for (auto it = this->begin(); it != this->end(); it ++) {
    if (!it->second.consistent_type) {
      const auto subtrajs = it->second.split();
      result.insert(result.end(), subtrajs.begin(), subtrajs.end());
      to_erase.push_back(it);
    }
  }

  for (const auto k : to_erase)
    erase(k);

  for (const auto &traj : result)
    add(traj, traj[0].id);
}

inline feature_curve_set_t feature_curve_set_t::intercept(int t0, int t1) const
{
  feature_curve_set_t result;
  for (const auto &kv : * this) {
    auto traj = kv.second.intercept(t0, t1);
    if (!traj.empty())
      // result[kv.first] = traj;
      result.insert({kv.first, traj});
  }
  return result;
}

} // namespace ftk

#endif
