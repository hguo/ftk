#ifndef __ftkCriticalPointTracker2DUnstructured_h
#define __ftkCriticalPointTracker2DUnstructured_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>

class vtkDataSet;

class ftkCriticalPointTracker2DUnstructured : public vtkUnstructuredGridAlgorithm
{
public:
  static ftkCriticalPointTracker2DUnstructured *New();
  vtkTypeMacro(ftkCriticalPointTracker2DUnstructured, vtkUnstructuredGridAlgorithm);

  vtkSetMacro(ZTimeScale, double);
  vtkGetMacro(ZTimeScale, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  ftkCriticalPointTracker2DUnstructured();
  ~ftkCriticalPointTracker2DUnstructured();

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkCriticalPointTracker2DUnstructured(const ftkCriticalPointTracker2DUnstructured&);
  void operator=(const ftkCriticalPointTracker2DUnstructured&);

private:
  double ZTimeScale;
  std::string InputVariable;

  int currentTimestep;
  
  std::shared_ptr<ftk::critical_point_tracker_2d_unstructured> tracker; 
  ftk::simplicial_unstructured_2d_mesh<> m2;
};

#endif
