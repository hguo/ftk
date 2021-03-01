#ifndef __ftkCappedWovenGradient2DSource_h
#define __ftkCappedWovenGradient2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkCappedWovenGradient2DSource : public vtkImageAlgorithm
{
public:
  static ftkCappedWovenGradient2DSource *New();
  vtkTypeMacro(ftkCappedWovenGradient2DSource, vtkImageAlgorithm);

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
  
  vtkSetMacro(Capping, double);
  vtkGetMacro(Capping, double);

protected:
  ftkCappedWovenGradient2DSource();
  ~ftkCappedWovenGradient2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkCappedWovenGradient2DSource(const ftkCappedWovenGradient2DSource&);
  void operator=(const ftkCappedWovenGradient2DSource&);

private:
  int DW, DH, DT;
  double ScalingFactor;
  double StartTime, TimeScale;
  double NoiseInjection;
  double Capping;
};

#endif
