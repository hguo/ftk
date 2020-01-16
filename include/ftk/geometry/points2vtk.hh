#ifndef _FTK_POINTS2VTK_HH
#define _FTK_POINTS2VTK_HH

#include <ftk/ftk_config.hh>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>

namespace ftk {

template <typename Array>
vtkSmartPointer<vtkPolyData>
points2vtk(const Array &array)
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();

  vtkIdType pid[1];

  for (auto i = 0; i < array.size(); i ++) {
    double p[3] = {array[i][0], array[i][1], array[i][2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  return polyData;
}

template <typename T>
vtkSmartPointer<vtkPolyData>
points2vtk(const std::vector<T>& array, int skip = 3)
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();

  vtkIdType pid[1];

  for (auto i = 0; i < array.size() / skip; i ++) {
    T p[3] = {array[i*skip], array[i*skip+1],
      skip == 2 ? T(0) : array[i*skip+2]};

    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  return polyData;
}

template <typename T>
vtkSmartPointer<vtkPolyData>
points2vtk(const std::vector<T>& array, const std::vector<T>& scalar, int skip = 3)
{
  auto polyData = points2vtk(array, skip);
  vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
  
  weights->SetNumberOfValues(scalar.size());
  for (auto i = 0; i < scalar.size(); i ++)
    weights->SetValue(i, scalar[i]);

  polyData->GetPointData()->SetScalars(weights);
  return polyData;
}

inline void write_vtp(const std::string& filename, vtkSmartPointer<vtkPolyData> polydata)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(polydata);
  writer->Write();
}

}
#endif
#endif
