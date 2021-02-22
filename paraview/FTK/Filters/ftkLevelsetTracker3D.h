#ifndef __ftkLevelsetTracker3D_h
#define __ftkLevelsetTracker3D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/contour_tracker_3d_regular.hh>

class vtkDataSet;

class ftkLevelsetTracker3D : public vtkImageAlgorithm
{
public:
  static ftkLevelsetTracker3D *New();
  vtkTypeMacro(ftkLevelsetTracker3D, vtkImageAlgorithm);
  
  vtkSetMacro(Threshold, double);
  vtkGetMacro(Threshold, double);

  vtkSetMacro(UseGPU, bool);
  vtkGetMacro(UseGPU, bool);

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
  double Threshold;
  bool UseGPU;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
 
private:
  ftk::contour_tracker_3d_regular tracker; 
  bool tracked = false;
};

#endif
