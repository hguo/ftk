#ifndef __vtkABCFlow3DSource_h
#define __vtkABCFlow3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkABCFlow3DSource : public vtkImageAlgorithm
{
public:
  static vtkABCFlow3DSource *New();
  vtkTypeMacro(vtkABCFlow3DSource, vtkImageAlgorithm);

  void setA(double a) {A = a;}
  void setB(double b) {B = b;}
  void setC(double c) {C = c;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setD(int d) {DD = d;}
  // void setScalingFactor(double s) {scale = s;}

protected:
  vtkABCFlow3DSource();
  ~vtkABCFlow3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkABCFlow3DSource(const vtkABCFlow3DSource&);
  void operator=(const vtkABCFlow3DSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
