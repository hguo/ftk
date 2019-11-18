#ifndef _FTK_CURVE2VTK_HH
#define _FTK_CURVE2VTK_HH

#include <ftk/ftk_config.hh>

#if FTK_HAVE_VTK

#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>

namespace ftk {

template <typename T>
vtkSmartPointer<vtkPolyData> 
curves2vtk(const std::vector<std::vector<T>>& curves, int skip = 3)
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  for (const auto &curve : curves) {
    for (int i = 0; i < curve.size() / skip; i ++) {
      T p[3] = {curve[i*skip], curve[i*skip+1], curve[i*skip+2]};
      points->InsertNextPoint(p);
    }
  }

  size_t nv = 0;
  for (const auto &curve :curves) {
    vtkSmartPointer<vtkPolyLine> polyLine = vtkPolyLine::New();
    polyLine->GetPointIds()->SetNumberOfIds(curve.size() / skip);
    for (int i = 0; i < curve.size() / skip; i ++)
      polyLine->GetPointIds()->SetId(i, i+nv);

    cells->InsertNextCell(polyLine);
    nv += curve.size() / skip;
  }

  polyData->SetPoints(points);
  polyData->SetLines(cells);

  return polyData;
}

template <typename T>
vtkSmartPointer<vtkPolyData>
curves2vtk(const std::vector<std::vector<std::vector<T>>>& curve_groups, int skip = 3)
{
  std::vector<std::vector<T>> curves;
  for (const auto &g : curve_groups)
    for (const auto &c : g)
      curves.push_back(c);

  return curves2vtk(curves, skip);
}

template <typename T>
void write_curves_vtk(const std::vector<std::vector<std::vector<T>>>& curve_groups, const std::string& filename, int skip = 3)
{
  auto poly = curves2vtk(curve_groups, skip);
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(poly);
  writer->Write();
}

template <typename T>
void write_curves_vtk(const std::vector<std::vector<T>>& curves, const std::string& filename, int skip = 3)
{
  auto poly = curves2vtk(curves, skip);
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(poly);
  writer->Write();
}

}

#endif

#endif
