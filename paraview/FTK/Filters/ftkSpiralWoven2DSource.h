#ifndef __ftkSpiralWoven2DSource_h
#define __ftkSpiralWoven2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkSpiralWoven2DSource : public vtkImageAlgorithm
{
public:
  static ftkSpiralWoven2DSource *New();
  vtkTypeMacro(ftkSpiralWoven2DSource, vtkImageAlgorithm);

  vtkSetMacro(DW, int);
  vtkGetMacro(DW, int);

  vtkSetMacro(DH, int);
  vtkGetMacro(DH, int);

  vtkSetMacro(DT, int);
  vtkGetMacro(DT, int);

  vtkSetMacro(ScalingFactor, double);
  vtkGetMacro(ScalingFactor, double);

  vtkSetMacro(StartTime, double);
  vtkGetMacro(StartTime, double);

  vtkSetMacro(TimeScale, double);
  vtkGetMacro(TimeScale, double);

  vtkSetMacro(NoiseInjection, double);
  vtkGetMacro(NoiseInjection, double);

protected:
  ftkSpiralWoven2DSource();
  ~ftkSpiralWoven2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkSpiralWoven2DSource(const ftkSpiralWoven2DSource&);
  void operator=(const ftkSpiralWoven2DSource&);

private:
  int DW, DH, DT;
  double ScalingFactor;
  double StartTime, TimeScale;
  double NoiseInjection;
};

#endif
