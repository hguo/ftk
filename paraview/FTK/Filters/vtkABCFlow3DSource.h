#ifndef __vtkABCFlow3DSource_h
#define __vtkABCFlow3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class vtkABCFlow3DSource : public vtkImageAlgorithm
{
public:
  static vtkABCFlow3DSource *New();
  vtkTypeMacro(vtkABCFlow3DSource, vtkImageAlgorithm);

  vtkSetMacro(A, double);
  vtkSetMacro(B, double);
  vtkSetMacro(C, double);
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DD, int);

protected:
  vtkABCFlow3DSource();
  ~vtkABCFlow3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  vtkABCFlow3DSource(const vtkABCFlow3DSource&);
  void operator=(const vtkABCFlow3DSource&);

private:
  double A, B, C;
  int DW, DH, DD;
  // double scale;
};

#endif
