#ifndef __ftkLevelsetTracker2D_h
#define __ftkLevelsetTracker2D_h

#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <ftk/filters/contour_tracker_2d_regular.hh>

class vtkDataSet;

class ftkLevelsetTracker2D : public vtkImageAlgorithm
{
public:
  static ftkLevelsetTracker2D *New();
  vtkTypeMacro(ftkLevelsetTracker2D, vtkImageAlgorithm);

  vtkSetMacro(InputVariable, std::string);
  vtkGetMacro(InputVariable, std::string);
  
  vtkSetMacro(Threshold, double);
  vtkGetMacro(Threshold, double);

protected:
  ftkLevelsetTracker2D();
  ~ftkLevelsetTracker2D();

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  ftkLevelsetTracker2D(const ftkLevelsetTracker2D&);
  void operator=(const ftkLevelsetTracker2D&);

private:
  double Threshold;
  std::string InputVariable;

  int currentTimestep;
  int inputDataComponents;
  
  ftk::contour_tracker_2d_regular tracker;
};

#endif
