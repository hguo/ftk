#ifndef __ftkMerger2DSource_h
#define __ftkMerger2DSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

class vtkDataSet;

class ftkMerger2DSource : public vtkImageAlgorithm
{
public:
  static ftkMerger2DSource *New();
  vtkTypeMacro(ftkMerger2DSource, vtkImageAlgorithm);

  vtkSetMacro(TimeScale, double);
  vtkSetMacro(DW, int);
  vtkSetMacro(DH, int);
  vtkSetMacro(DT, int);

protected:
  ftkMerger2DSource();
  ~ftkMerger2DSource();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillOutputPortInformation(int, vtkInformation*) override;

private:
  ftkMerger2DSource(const ftkMerger2DSource&);
  void operator=(const ftkMerger2DSource&);

private:
  double TimeScale;
  int DW, DH, DT;
};

#endif
