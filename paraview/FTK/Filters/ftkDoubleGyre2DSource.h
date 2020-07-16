#ifndef __ftkDoubleGyre2DSource_h
#define __ftkDoubleGyre2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkDoubleGyre2DSource : public vtkImageAlgorithm
{
public:
  static ftkDoubleGyre2DSource *New();
  vtkTypeMacro(ftkDoubleGyre2DSource, vtkImageAlgorithm);

  vtkSetMacro(A, double);
  vtkSetMacro(Omega, double);
  vtkSetMacro(Epsilon, double);
  vtkSetMacro(StartTime, double);
  vtkSetMacro(TimeScale, double);
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);

protected:
  ftkDoubleGyre2DSource();
  ~ftkDoubleGyre2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkDoubleGyre2DSource(const ftkDoubleGyre2DSource&);
  void operator=(const ftkDoubleGyre2DSource&);

private:
  double A, Omega, Epsilon;
  double StartTime, TimeScale;
  int DW, DH, DT;
};

#endif
