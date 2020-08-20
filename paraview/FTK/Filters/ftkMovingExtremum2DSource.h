#ifndef __ftkMovingExtremum2DSource_h
#define __ftkMovingExtremum2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkMovingExtremum2DSource : public vtkImageAlgorithm
{
public:
  static ftkMovingExtremum2DSource *New();
  vtkTypeMacro(ftkMovingExtremum2DSource, vtkImageAlgorithm);

  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);
  vtkSetVector2Macro(X0, double);
  vtkSetVector2Macro(Dir, double);

protected:
  ftkMovingExtremum2DSource();
  ~ftkMovingExtremum2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkMovingExtremum2DSource(const ftkMovingExtremum2DSource&);
  void operator=(const ftkMovingExtremum2DSource&);

private:
  int DW, DH, DT;
  double X0[2], Dir[2];
};

#endif
