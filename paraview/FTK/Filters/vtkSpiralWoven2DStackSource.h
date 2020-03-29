#ifndef __vtkSpiralWoven2DStackSource_h
#define __vtkSpiralWoven2DStackSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkSpiralWoven2DStackSource : public vtkImageAlgorithm
{
public:
  static vtkSpiralWoven2DStackSource *New();
  vtkTypeMacro(vtkSpiralWoven2DStackSource, vtkImageAlgorithm);

  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}
  void setScalingFactor(double s) {scale = s;}

protected:
  vtkSpiralWoven2DStackSource();
  ~vtkSpiralWoven2DStackSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkSpiralWoven2DStackSource(const vtkSpiralWoven2DStackSource&);
  void operator=(const vtkSpiralWoven2DStackSource&);

private:
  int DW, DH, DT;
  double scale;
};

#endif
