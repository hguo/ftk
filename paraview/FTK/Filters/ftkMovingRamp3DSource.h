#ifndef __ftkMovingRamp3DSource_h
#define __ftkMovingRamp3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkMovingRamp3DSource : public vtkImageAlgorithm
{
public:
  static ftkMovingRamp3DSource *New();
  vtkTypeMacro(ftkMovingRamp3DSource, vtkImageAlgorithm);

  vtkSetMacro(DW, int);
  vtkGetMacro(DW, int);

  vtkSetMacro(DH, int);
  vtkGetMacro(DH, int);

  vtkSetMacro(DD, int);
  vtkGetMacro(DD, int);

  vtkSetMacro(DT, int);
  vtkGetMacro(DT, int);

  vtkSetMacro(X0, double);
  vtkGetMacro(X0, double);

  vtkSetMacro(Rate, double);
  vtkGetMacro(Rate, double);

protected:
  ftkMovingRamp3DSource();
  ~ftkMovingRamp3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkMovingRamp3DSource(const ftkMovingRamp3DSource&);
  void operator=(const ftkMovingRamp3DSource&);

private:
  int DW, DH, DD, DT;
  double X0, Rate;
};

#endif
