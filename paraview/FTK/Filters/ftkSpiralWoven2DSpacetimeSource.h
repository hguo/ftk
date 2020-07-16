#ifndef __ftkSpiralWoven2DSpacetimeSource_h
#define __ftkSpiralWoven2DSpacetimeSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkSpiralWoven2DSpacetimeSource : public vtkImageAlgorithm
{
public:
  static ftkSpiralWoven2DSpacetimeSource *New();
  vtkTypeMacro(ftkSpiralWoven2DSpacetimeSource, vtkImageAlgorithm);
  
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);
  vtkSetMacro(ScalingFactor, double);

protected:
  ftkSpiralWoven2DSpacetimeSource();
  ~ftkSpiralWoven2DSpacetimeSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkSpiralWoven2DSpacetimeSource(const ftkSpiralWoven2DSpacetimeSource&);
  void operator=(const ftkSpiralWoven2DSpacetimeSource&);

private:
  int DW, DH, DT;
  double ScalingFactor;
};

#endif
