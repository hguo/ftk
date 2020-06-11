#ifndef __vtkCriticalPointTracker2DSpaceTime_h
#define __vtkCriticalPointTracker2DSpaceTime_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkCriticalPointTracker2DSpaceTime : public vtkImageAlgorithm
{
public:
  static vtkCriticalPointTracker2DSpaceTime *New();
  vtkTypeMacro(vtkCriticalPointTracker2DSpaceTime, vtkImageAlgorithm);

  void SetUseGPU(bool);
  void SetGaussianKernelSize(double);

protected:
  vtkCriticalPointTracker2DSpaceTime();
  ~vtkCriticalPointTracker2DSpaceTime();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  int TrackCriticalPoints2DSpaceTime(vtkImageData*, vtkPolyData*); // the input is a stack of 2D slices

private:
  vtkCriticalPointTracker2DSpaceTime(const vtkCriticalPointTracker2DSpaceTime&);
  void operator=(const vtkCriticalPointTracker2DSpaceTime&);

private:
  bool bUseGPU;
  double dGaussianKernelSize;
};

#endif
