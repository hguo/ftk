#ifndef __vtkSpiral2DStackSource_h
#define __vtkSpiral2DStackSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkSpiral2DStackSource : public vtkImageAlgorithm
{
public:
  static vtkSpiral2DStackSource *New();
  vtkTypeMacro(vtkSpiral2DStackSource, vtkImageAlgorithm);

  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}
  void setScalingFactor(double s) {scale = s;}

protected:
  vtkSpiral2DStackSource();
  ~vtkSpiral2DStackSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkSpiral2DStackSource(const vtkSpiral2DStackSource&);
  void operator=(const vtkSpiral2DStackSource&);

private:
  int DW, DH, DT;
  double scale;
};

#endif
