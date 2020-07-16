#ifndef __ftkCriticalPointTracker2DSpacetime_h
#define __ftkCriticalPointTracker2DSpacetime_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkCriticalPointTracker2DSpacetime : public vtkImageAlgorithm
{
public:
  static ftkCriticalPointTracker2DSpacetime *New();
  vtkTypeMacro(ftkCriticalPointTracker2DSpacetime, vtkImageAlgorithm);
  
  vtkSetMacro(UseGPU, bool);
  vtkSetMacro(GaussianKernelSize, double);

protected:
  ftkCriticalPointTracker2DSpacetime();
  ~ftkCriticalPointTracker2DSpacetime();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  int TrackCriticalPoints2DSpacetime(vtkImageData*, vtkPolyData*); // the input is a stack of 2D slices

private:
  ftkCriticalPointTracker2DSpacetime(const ftkCriticalPointTracker2DSpacetime&);
  void operator=(const ftkCriticalPointTracker2DSpacetime&);

private:
  bool UseGPU;
  double GaussianKernelSize;
};

#endif
