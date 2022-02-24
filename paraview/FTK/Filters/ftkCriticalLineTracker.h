#ifndef __ftkCriticalLineTracker_h
#define __ftkCriticalLineTracker_h

#include "vtkAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkRectilinearGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include <ftk/filters/critical_line_tracker_3d_regular.hh>
// #include <ftk/filters/critical_line_tracker_3d_unstructured.hh>

class vtkDataSet;

class ftkCriticalLineTracker : public vtkAlgorithm
{
public:
  static ftkCriticalLineTracker *New();
  vtkTypeMacro(ftkCriticalLineTracker, vtkAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);
  
  vtkSetMacro(ComputeDegrees, bool);
  vtkGetMacro(ComputeDegrees, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);
  
  vtkSetMacro(ZTimeScale, double);
  vtkGetMacro(ZTimeScale, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);
  
protected:
  ftkCriticalLineTracker();
  ~ftkCriticalLineTracker();

  // int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  // int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;
  int FillInputPortInformation(int, vtkInformation*) override;

private:
  int RequestData_vti(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int RequestData_vtu(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
  ftkCriticalLineTracker(const ftkCriticalLineTracker&);
  void operator=(const ftkCriticalLineTracker&);

private:
  bool UseGPU;
  bool ComputeDegrees;
  double GaussianKernelSize;
  double ZTimeScale;
  std::string InputVariable, InputVariable2;

  int currentTimestep;
  int inputDataComponents;
  
  std::shared_ptr<ftk::critical_line_tracker_3d_regular> tlpr;
  std::shared_ptr<ftk::critical_line_tracker> tlp;
  
  // ftk::simplicial_unstructured_2d_mesh<> m2u;
};

#endif
