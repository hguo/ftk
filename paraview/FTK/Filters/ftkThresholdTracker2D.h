#ifndef __ftkThresholdTracker2D_h
#define __ftkThresholdTracker2D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/threshold_tracker.hh>

class vtkDataSet;

class ftkThresholdTracker2D : public vtkImageAlgorithm
{
public:
  static ftkThresholdTracker2D *New();
  vtkTypeMacro(ftkThresholdTracker2D, vtkImageAlgorithm);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);
  
  vtkSetMacro(Threshold, double);

protected:
  ftkThresholdTracker2D();
  ~ftkThresholdTracker2D();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  ftkThresholdTracker2D(const ftkThresholdTracker2D&);
  void operator=(const ftkThresholdTracker2D&);

private:
  double Threshold;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  ftk::threshold_tracker<> tracker; 
};

#endif
