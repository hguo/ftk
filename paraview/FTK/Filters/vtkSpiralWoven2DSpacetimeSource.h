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
  
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);
  vtkSetMacro(ScalingFactor, double);

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
  double ScalingFactor;
};

#endif
