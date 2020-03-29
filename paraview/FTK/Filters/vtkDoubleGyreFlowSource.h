#ifndef __vtkDoubleGyreFlowSource_h
#define __vtkDoubleGyreFlowSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkDoubleGyreFlowSource : public vtkImageAlgorithm
{
public:
  static vtkDoubleGyreFlowSource *New();
  vtkTypeMacro(vtkDoubleGyreFlowSource, vtkImageAlgorithm);

  void setA(double a) {A = a;}
  void setB(double b) {B = b;}
  void setC(double c) {C = c;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setD(int d) {DD = d;}
  // void setScalingFactor(double s) {scale = s;}

protected:
  vtkDoubleGyreFlowSource();
  ~vtkDoubleGyreFlowSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkDoubleGyreFlowSource(const vtkDoubleGyreFlowSource&);
  void operator=(const vtkDoubleGyreFlowSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
