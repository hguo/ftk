#ifndef __ftkMovingExtremum3DSource_h
#define __ftkMovingExtremum3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkMovingExtremum3DSource : public vtkImageAlgorithm
{
public:
  static ftkMovingExtremum3DSource *New();
  vtkTypeMacro(ftkMovingExtremum3DSource, vtkImageAlgorithm);

  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DD, int);
  vtkSetMacro(DT, int);
  vtkSetVector3Macro(X0, double);
  vtkSetVector3Macro(Dir, double);

protected:
  ftkMovingExtremum3DSource();
  ~ftkMovingExtremum3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkMovingExtremum3DSource(const ftkMovingExtremum3DSource&);
  void operator=(const ftkMovingExtremum3DSource&);

private:
  int DW, DH, DD, DT;
  double X0[3], Dir[3];
};

#endif
