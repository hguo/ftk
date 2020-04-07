#ifndef __vtkDoubleGyreVectorField2DSource_h
#define __vtkDoubleGyreVectorField2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkDoubleGyreVectorField2DSource : public vtkImageAlgorithm
{
public:
  static vtkDoubleGyreVectorField2DSource *New();
  vtkTypeMacro(vtkDoubleGyreVectorField2DSource, vtkImageAlgorithm);

  void setA(double a) {A = a;}
  void setOmega(double omega) {Omega = omega;}
  void setEpsilon(double epsilon) {Epsilon = epsilon;}
  void setStartTime(double t) {StartTime = t;}
  void setTimeScale(double ts) {TimeScale = ts;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}

protected:
  vtkDoubleGyreVectorField2DSource();
  ~vtkDoubleGyreVectorField2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkDoubleGyreVectorField2DSource(const vtkDoubleGyreVectorField2DSource&);
  void operator=(const vtkDoubleGyreVectorField2DSource&);

private:
  double A, Omega, Epsilon;
  double StartTime, TimeScale;
  int DW, DH, DT;
};

#endif
