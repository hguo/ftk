#ifndef __ftkCriticalPointTracker_h
#define __ftkCriticalPointTracker_h

#include "vtkAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkRectilinearGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>

class vtkDataSet;

class ftkCriticalPointTracker : public vtkAlgorithm
{
public:
  static ftkCriticalPointTracker *New();
  vtkTypeMacro(ftkCriticalPointTracker, vtkAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);
  
  vtkSetMacro(ZTimeScale, double);
  vtkGetMacro(ZTimeScale, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  ftkCriticalPointTracker();
  ~ftkCriticalPointTracker();

  // int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  // int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;
  int FillInputPortInformation(int, vtkInformation*) override;

private:
  int RequestData_vti(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int RequestData_vtr(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int RequestData_vts(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int RequestData_vtu(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
  ftkCriticalPointTracker(const ftkCriticalPointTracker&);
  void operator=(const ftkCriticalPointTracker&);

private:
  bool UseGPU;
  double GaussianKernelSize;
  double ZTimeScale;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  std::shared_ptr<ftk::critical_point_tracker_regular> tcpr;
  std::shared_ptr<ftk::critical_point_tracker_2d_unstructured> tcp2du;
  
  ftk::simplicial_unstructured_2d_mesh<> m2u;
};

#endif
