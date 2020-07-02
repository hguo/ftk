#ifndef __vtkCriticalPointTracker2DSpacetime_h
#define __vtkCriticalPointTracker2DSpacetime_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkCriticalPointTracker2DSpacetime : public vtkImageAlgorithm
{
public:
  static vtkCriticalPointTracker2DSpacetime *New();
  vtkTypeMacro(vtkCriticalPointTracker2DSpacetime, vtkImageAlgorithm);
  
  vtkSetMacro(UseGPU, bool);
  vtkSetMacro(GaussianKernelSize, double);

protected:
  vtkCriticalPointTracker2DSpacetime();
  ~vtkCriticalPointTracker2DSpacetime();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  int TrackCriticalPoints2DSpacetime(vtkImageData*, vtkPolyData*); // the input is a stack of 2D slices

private:
  vtkCriticalPointTracker2DSpacetime(const vtkCriticalPointTracker2DSpacetime&);
  void operator=(const vtkCriticalPointTracker2DSpacetime&);

private:
  bool UseGPU;
  double GaussianKernelSize;
};

#endif
