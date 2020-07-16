#ifndef __ftkCriticalPointTracker3D_h
#define __ftkCriticalPointTracker3D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_3d_regular.hh>

class vtkDataSet;

class ftkCriticalPointTracker3D : public vtkImageAlgorithm
{
public:
  static ftkCriticalPointTracker3D *New();
  vtkTypeMacro(ftkCriticalPointTracker3D, vtkImageAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  ftkCriticalPointTracker3D();
  ~ftkCriticalPointTracker3D();

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkCriticalPointTracker3D(const ftkCriticalPointTracker3D&);
  void operator=(const ftkCriticalPointTracker3D&);

private:
  bool UseGPU;
  double GaussianKernelSize;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  ftk::critical_point_tracker_3d_regular tracker; 
};

#endif
