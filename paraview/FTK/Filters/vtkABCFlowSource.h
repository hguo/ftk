#ifndef __vtkABCFlowSource_h
#define __vtkABCFlowSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkABCFlowSource : public vtkImageAlgorithm
{
public:
  static vtkABCFlowSource *New();
  vtkTypeMacro(vtkABCFlowSource, vtkImageAlgorithm);

  void setA(double a) {A = a;}
  void setB(double b) {B = b;}
  void setC(double c) {C = c;}
  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setD(int d) {DD = d;}
  // void setScalingFactor(double s) {scale = s;}

protected:
  vtkABCFlowSource();
  ~vtkABCFlowSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkABCFlowSource(const vtkABCFlowSource&);
  void operator=(const vtkABCFlowSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
