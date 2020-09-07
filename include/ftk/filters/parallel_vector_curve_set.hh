#ifndef _FTK_PARALLEL_VECTOR_CURVE_SET_T_HH
#define _FTK_PARALLEL_VECTOR_CURVE_SET_T_HH

#include <ftk/filters/parallel_vector_curve.hh>
#include <map>

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#endif

namespace ftk {

struct parallel_vector_curve_set_t : public std::map<int, parallel_vector_curve_t>
{
  void foreach(std::function<void(int, const parallel_vector_curve_t&)> f) const {for (const auto& kv : *this) f(kv.first, kv.second); }
  void foreach(std::function<void(int, parallel_vector_curve_t&)> f) {for (auto& kv : *this) f(kv.first, kv.second); }
  
  int add(const parallel_vector_curve_t&);

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> to_vtp() const;
#endif

protected:
  int get_new_id() const {
    if (empty()) return 0; 
    else return rbegin()->first + 1;
  }
};

////
int parallel_vector_curve_set_t::add(const parallel_vector_curve_t& curve)
{
  const int id = get_new_id();
  insert(std::pair<int, parallel_vector_curve_t>(id, curve));
  // at(id).relabel(id);
  return id;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> parallel_vector_curve_set_t::to_vtp() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> lines = vtkCellArray::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();
  
  foreach([&](int, const parallel_vector_curve_t& curve) {
    for (auto i = 0; i < curve.size(); i ++) {
      double p[3] = {curve[i].x[0], curve[i].x[1], curve[i].x[2]};
      points->InsertNextPoint(p);
    }
  });

  size_t nv = 0;
  foreach([&](int, const parallel_vector_curve_t& curve) {
    if (curve.empty()) { // skip
    } else if (curve.size() < 2) { // isolated vertex
      vtkSmartPointer<vtkVertex> obj = vtkVertex::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      verts->InsertNextCell(obj);
    } else { // lines
      vtkSmartPointer<vtkPolyLine> obj = vtkPolyLine::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      lines->InsertNextCell(obj);
    }
    nv += curve.size();
  });

  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->SetVerts(verts);

  if (1) { // ids
    vtkSmartPointer<vtkUnsignedIntArray> ids = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ids->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](int id, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        ids->SetValue(i ++, curve[j].id);
    });
    ids->SetName("id");
    polyData->GetPointData()->AddArray(ids);
  }
  
  if (1) { // lambda
    vtkSmartPointer<vtkDoubleArray> lambdas = vtkSmartPointer<vtkDoubleArray>::New();
    lambdas->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](int, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        lambdas->SetValue(i ++, curve[j].lambda);
    });
    lambdas->SetName("lambda");
    polyData->GetPointData()->AddArray(lambdas);
  }
  
  if (1) { // v
    vtkSmartPointer<vtkDoubleArray> v = vtkSmartPointer<vtkDoubleArray>::New();
    v->SetNumberOfComponents(3);
    v->SetNumberOfTuples(nv);
    size_t i = 0;
    foreach([&](int, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        v->SetTuple3(i ++, curve[j].v[0], curve[j].v[1], curve[j].v[2]);
    });
    v->SetName("v");
    polyData->GetPointData()->AddArray(v);
  }
  
  if (1) { // w
    vtkSmartPointer<vtkDoubleArray> w = vtkSmartPointer<vtkDoubleArray>::New();
    w->SetNumberOfComponents(3);
    w->SetNumberOfTuples(nv);
    size_t i = 0;
    foreach([&](int, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        w->SetTuple3(i ++, curve[j].w[0], curve[j].w[1], curve[j].w[2]);
    });
    w->SetName("w");
    polyData->GetPointData()->AddArray(w);
  }

  if (1) { // condition numbers
    vtkSmartPointer<vtkDoubleArray> conds = vtkSmartPointer<vtkDoubleArray>::New();
    conds->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](int, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        conds->SetValue(i ++, curve[j].cond);
    });
    conds->SetName("cond");
    polyData->GetPointData()->AddArray(conds);
  }
  
  if (1) { // time
    vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
    time->SetNumberOfValues(nv);
    size_t i = 0;
    foreach([&](int, const parallel_vector_curve_t& curve) {
      for (auto j = 0; j < curve.size(); j ++)
        time->SetValue(i ++, curve[j].t);
    });
    time->SetName("time");
    polyData->GetPointData()->AddArray(time);
  }

  return polyData;
}
#endif

}

#endif
