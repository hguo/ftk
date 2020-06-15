#ifndef __ftkTornadoFlow3DSource_h
#define __ftkTornadoFlow3DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkTornadoFlow3DSource : public vtkImageAlgorithm
{
public:
  static ftkTornadoFlow3DSource *New();
  vtkTypeMacro(ftkTornadoFlow3DSource, vtkImageAlgorithm);

  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DD, int);

protected:
  ftkTornadoFlow3DSource();
  ~ftkTornadoFlow3DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkTornadoFlow3DSource(const ftkTornadoFlow3DSource&);
  void operator=(const ftkTornadoFlow3DSource&);

private:
  int DW, DH, DD;
  // double scale;
};

#endif
