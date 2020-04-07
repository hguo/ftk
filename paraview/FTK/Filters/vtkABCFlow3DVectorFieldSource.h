#ifndef __vtkABCFlow3DVectorFieldSource_h
#define __vtkABCFlow3DVectorFieldSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkABCFlow3DVectorFieldSource : public vtkImageAlgorithm
{
public:
  static vtkABCFlow3DVectorFieldSource *New();
  vtkTypeMacro(vtkABCFlow3DVectorFieldSource, vtkImageAlgorithm);

  void setA(double a) {A = a;}
  void setB(double b) {B = b;}
  void setC(double c) {C = c;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setD(int d) {DD = d;}
  // void setScalingFactor(double s) {scale = s;}

protected:
  vtkABCFlow3DVectorFieldSource();
  ~vtkABCFlow3DVectorFieldSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkABCFlow3DVectorFieldSource(const vtkABCFlow3DVectorFieldSource&);
  void operator=(const vtkABCFlow3DVectorFieldSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
