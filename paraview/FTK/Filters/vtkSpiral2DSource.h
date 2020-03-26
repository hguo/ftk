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

  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}
  void setScalingFactor(double s) {scale = s;}

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
  int DW, DH, DT;
  double scale;
};

#endif
