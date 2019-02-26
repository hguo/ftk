#ifndef _FTK_CURVE2VTK_HH
#define _FTK_CURVE2VTK_HH

#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>

namespace ftk {

template <typename T>
vtkSmartPointer<vtkPolyData> 
curves2vtk(const std::vector<std::vector<T>>& curves)
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  for (const auto &curve : curves) {
    for (int i = 0; i < curve.size() / 3; i ++) {
      T p[3] = {curve[i*3], curve[i*3+1], curve[i*3+2]};
      points->InsertNextPoint(p);
    }
  }

  size_t nv = 0;
  for (const auto &curve :curves) {
    vtkSmartPointer<vtkPolyLine> polyLine = vtkPolyLine::New();
    polyLine->GetPointIds()->SetNumberOfIds(curve.size() / 3);
    for (int i = 0; i < curve.size() / 3; i ++)
      polyLine->GetPointIds()->SetId(i, i+nv);

    cells->InsertNextCell(polyLine);
    nv += curve.size() / 3;
  }

  polyData->SetPoints(points);
  polyData->SetLines(cells);

  return polyData;
}

template <typename T>
vtkSmartPointer<vtkPolyData>
curves2vtk(const std::vector<std::vector<std::vector<float>>>& curve_groups)
{
  std::vector<std::vector<T>> curves;
  for (const auto &g : curve_groups)
    curves.push_back(g);
}

template <typename T>
void write_curves_vtk(const std::vector<std::vector<std::vector<float>>>& curve_groups, const std::string& filename)
{
  auto poly = curves2vtk(curve_groups);
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(poly);
  writer->Write();
}

template <typename T>
void write_curves_vtk(const std::vector<std::vector<float>>& curves, const std::string& filename)
{
  auto poly = curves2vtk(curves);
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(poly);
  writer->Write();
}

}

#endif
