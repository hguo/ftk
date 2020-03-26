#ifndef __vtkCriticalPoint2DTracker_h
#define __vtkCriticalPoint2DTracker_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkCriticalPoint2DTracker : public vtkImageAlgorithm
{
public:
  static vtkCriticalPoint2DTracker *New();
  vtkTypeMacro(vtkCriticalPoint2DTracker, vtkImageAlgorithm);

  void SetUseGPU(bool);
  void SetGaussianKernelSize(double);

protected:
  vtkCriticalPoint2DTracker();
  ~vtkCriticalPoint2DTracker();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  int TrackCriticalPoints2D(vtkImageData*, vtkPolyData*);

private:
  vtkCriticalPoint2DTracker(const vtkCriticalPoint2DTracker&);
  void operator=(const vtkCriticalPoint2DTracker&);

private:
  bool bUseGPU;
  double dGaussianKernelSize;
};

#endif
