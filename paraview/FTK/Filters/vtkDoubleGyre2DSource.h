#ifndef __vtkDoubleGyre2DSource_h
#define __vtkDoubleGyre2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkDoubleGyre2DSource : public vtkImageAlgorithm
{
public:
  static vtkDoubleGyre2DSource *New();
  vtkTypeMacro(vtkDoubleGyre2DSource, vtkImageAlgorithm);

  vtkSetMacro(A, double);
  vtkSetMacro(Omega, double);
  vtkSetMacro(Epsilon, double);
  vtkSetMacro(StartTime, double);
  vtkSetMacro(TimeScale, double);
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);

protected:
  vtkDoubleGyre2DSource();
  ~vtkDoubleGyre2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkDoubleGyre2DSource(const vtkDoubleGyre2DSource&);
  void operator=(const vtkDoubleGyre2DSource&);

private:
  double A, Omega, Epsilon;
  double StartTime, TimeScale;
  int DW, DH, DT;
};

#endif
