#ifndef _FTK_CURVE2VTK_HH
#define _FTK_CURVE2VTK_HH

#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>

namespace ftk {

vtkSmartPointer<vtkPolyData> 
curves2vtk(const std::vector<std::vector<float>>& curves)
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPoints::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  for (const auto &curve : curves) {
    for (int i = 0; i < curves.size() / 3; i ++) {
      float p[3] = {curve[i*3], curve[i*3+1], curve[i*3+2]};
      points->InsertNextPoint(p);
    }
  }

  size_t nv = 0;
  for (const auto &curve :curves) {
    vtkSmartPointer<vtkPolyLine> = vtkPolyLine::New();
    polyLine->GetPointIds()->SetNumberOfIds(curve.size() / 3);
    for (int i = 0; i < curve.size() / 3; i ++)
      polyLine->GetPointIds()->SetId(j, j+nv);

    cells->InsertNextCell(polyLine);
    nv += curve.size() / 3;
  }

  polyData->SetPoints(points);
  polyData->SetLines(cells);

  return polyData;
}

}

#endif
