#ifndef __vtkSpiralWoven2DSpacetimeSource_h
#define __vtkSpiralWoven2DSpacetimeSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkSpiralWoven2DSpacetimeSource : public vtkImageAlgorithm
{
public:
  static vtkSpiralWoven2DSpacetimeSource *New();
  vtkTypeMacro(vtkSpiralWoven2DSpacetimeSource, vtkImageAlgorithm);

  void setW(int w) {DW = w;}
  void setH(int h) {DH = h;}
  void setT(int t) {DT = t;}
  void setScalingFactor(double s) {scale = s;}

protected:
  vtkSpiralWoven2DSpacetimeSource();
  ~vtkSpiralWoven2DSpacetimeSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkSpiralWoven2DSpacetimeSource(const vtkSpiralWoven2DSpacetimeSource&);
  void operator=(const vtkSpiralWoven2DSpacetimeSource&);

private:
  int DW, DH, DT;
  double scale;
};

#endif
