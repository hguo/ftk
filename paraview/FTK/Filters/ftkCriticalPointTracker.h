#ifndef __ftkCriticalPointTracker_h
#define __ftkCriticalPointTracker_h

#include "vtkAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_2d_regular.hh>

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

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*); // override;
  int FillOutputPortInformation(int, vtkInformation*) override;

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
  
  ftk::critical_point_tracker_2d_regular tracker; 
};

#endif
