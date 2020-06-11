#ifndef __vtkSpiralWoven2DSource_h
#define __vtkSpiralWoven2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkSpiralWoven2DSource : public vtkImageAlgorithm
{
public:
  static vtkSpiralWoven2DSource *New();
  vtkTypeMacro(vtkSpiralWoven2DSource, vtkImageAlgorithm);

  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}
  void setScalingFactor(double s) {scale = s;}
  void setStartTime(double t) {StartTime = t;}
  void setTimeScale(double ts) {TimeScale = ts;}

protected:
  vtkSpiralWoven2DSource();
  ~vtkSpiralWoven2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkSpiralWoven2DSource(const vtkSpiralWoven2DSource&);
  void operator=(const vtkSpiralWoven2DSource&);

private:
  int DW, DH, DT;
  double scale;
  double StartTime, TimeScale;
};

#endif
