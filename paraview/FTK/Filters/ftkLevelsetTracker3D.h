#ifndef __ftkLevelsetTracker3D_h
#define __ftkLevelsetTracker3D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/critical_point_tracker_3d_regular.hh>

class vtkDataSet;

class ftkLevelsetTracker3D : public vtkImageAlgorithm
{
public:
  static ftkLevelsetTracker3D *New();
  vtkTypeMacro(ftkLevelsetTracker3D, vtkImageAlgorithm);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

  vtkSetMacro(GaussianKernelSize, double);
  vtkGetMacro(GaussianKernelSize, double);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);

protected:
  ftkLevelsetTracker3D();
  ~ftkLevelsetTracker3D();

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkLevelsetTracker3D(const ftkLevelsetTracker3D&);
  void operator=(const ftkLevelsetTracker3D&);

private:
  bool UseGPU;
  double GaussianKernelSize;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  ftk::critical_point_tracker_3d_regular tracker; 
};

#endif
