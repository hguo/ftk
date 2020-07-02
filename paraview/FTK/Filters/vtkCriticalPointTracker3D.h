#ifndef __vtkCriticalPointTracker3D_h
#define __vtkCriticalPointTracker3D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_3d_regular.hh>

class vtkDataSet;

class vtkCriticalPointTracker3D : public vtkImageAlgorithm
{
public:
  static vtkCriticalPointTracker3D *New();
  vtkTypeMacro(vtkCriticalPointTracker3D, vtkImageAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  vtkCriticalPointTracker3D();
  ~vtkCriticalPointTracker3D();

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkCriticalPointTracker3D(const vtkCriticalPointTracker3D&);
  void operator=(const vtkCriticalPointTracker3D&);

private:
  bool UseGPU;
  double GaussianKernelSize;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  ftk::critical_point_tracker_3d_regular tracker; 
};

#endif
