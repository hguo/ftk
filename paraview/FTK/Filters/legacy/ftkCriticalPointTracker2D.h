#ifndef __ftkCriticalPointTracker2D_h
#define __ftkCriticalPointTracker2D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_2d_regular.hh>

class vtkDataSet;

class ftkCriticalPointTracker2D : public vtkImageAlgorithm
{
public:
  static ftkCriticalPointTracker2D *New();
  vtkTypeMacro(ftkCriticalPointTracker2D, vtkImageAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);
  
  vtkSetMacro(ZTimeScale, double);
  vtkGetMacro(ZTimeScale, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  ftkCriticalPointTracker2D();
  ~ftkCriticalPointTracker2D();

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkCriticalPointTracker2D(const ftkCriticalPointTracker2D&);
  void operator=(const ftkCriticalPointTracker2D&);

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
