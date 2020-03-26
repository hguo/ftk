#ifndef __vtkSpiral2DSource_h
#define __vtkSpiral2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkSpiral2DSource : public vtkImageAlgorithm
{
public:
  static vtkSpiral2DSource *New();
  vtkTypeMacro(vtkSpiral2DSource, vtkImageAlgorithm);

protected:
  vtkSpiral2DSource();
  ~vtkSpiral2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkSpiral2DSource(const vtkSpiral2DSource&);
  void operator=(const vtkSpiral2DSource&);

private:
  bool bUseGPU;
  double dGaussianKernelSize;
};

#endif
