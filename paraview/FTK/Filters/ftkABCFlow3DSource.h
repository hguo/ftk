#ifndef __ftkABCFlow3DSource_h
#define __ftkABCFlow3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkABCFlow3DSource : public vtkImageAlgorithm
{
public:
  static ftkABCFlow3DSource *New();
  vtkTypeMacro(ftkABCFlow3DSource, vtkImageAlgorithm);

  vtkSetMacro(A, double);
  vtkSetMacro(B, double);
  vtkSetMacro(C, double);
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DD, int);

protected:
  ftkABCFlow3DSource();
  ~ftkABCFlow3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkABCFlow3DSource(const ftkABCFlow3DSource&);
  void operator=(const ftkABCFlow3DSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
