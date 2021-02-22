#ifndef __ftkMovingDualRamp3DSource_h
#define __ftkMovingDualRamp3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkMovingDualRamp3DSource : public vtkImageAlgorithm
{
public:
  static ftkMovingDualRamp3DSource *New();
  vtkTypeMacro(ftkMovingDualRamp3DSource, vtkImageAlgorithm);

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
  
  vtkSetMacro(Offset, double);
  vtkGetMacro(Offset, double);

protected:
  ftkMovingDualRamp3DSource();
  ~ftkMovingDualRamp3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkMovingDualRamp3DSource(const ftkMovingDualRamp3DSource&);
  void operator=(const ftkMovingDualRamp3DSource&);

private:
  int DW, DH, DD, DT;
  double X0, Rate, Offset;
};

#endif
