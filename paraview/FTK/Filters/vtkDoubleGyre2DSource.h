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

  void setA(double a) {A = a;}
  void setOmega(double omega) {Omega = omega;}
  void setEpsilon(double epsilon) {Epsilon = epsilon;}
  void setStartTime(double t) {StartTime = t;}
  void setTimeScale(double ts) {TimeScale = ts;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}

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
