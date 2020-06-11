#ifndef __vtkCriticalPointTracker2D_h
#define __vtkCriticalPointTracker2D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkCriticalPointTracker2D : public vtkImageAlgorithm
{
public:
  static vtkCriticalPointTracker2D *New();
  vtkTypeMacro(vtkCriticalPointTracker2D, vtkImageAlgorithm);

  void SetUseGPU(bool);
  void SetGaussianKernelSize(double);

protected:
  vtkCriticalPointTracker2D();
  ~vtkCriticalPointTracker2D();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  int TrackCriticalPoints2DSpaceTime(vtkImageData*, vtkPolyData*); // the input is a stack of 2D slices

private:
  vtkCriticalPointTracker2D(const vtkCriticalPointTracker2D&);
  void operator=(const vtkCriticalPointTracker2D&);

private:
  bool bUseGPU;
  double dGaussianKernelSize;


  int step;
};

#endif
