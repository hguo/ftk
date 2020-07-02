#include "ftkLevelsetTracker2D.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkImageData.h"
#include "vtkSphereSource.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(ftkLevelsetTracker2D);

ftkLevelsetTracker2D::ftkLevelsetTracker2D()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  currentTimestep = 0;
  inputDataComponents = 0;
}

ftkLevelsetTracker2D::~ftkLevelsetTracker2D()
{
}

int ftkLevelsetTracker2D::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  return 1;
}

int ftkLevelsetTracker2D::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

int ftkLevelsetTracker2D::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  fprintf(stderr, "requesting data.\n");

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );

  fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
  tracker.set_threshold(Threshold);

  ftk::ndarray<double> field_data;
  field_data.from_vtk_image_data(input, InputVariable);

  tracker.push_scalar_field_data_snapshot(field_data);
  tracker.advance_timestep();

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    tracker.finalize();
    // auto output = tracker.get_traced_levelset_vtk();
  }

  currentTimestep ++;
  return 1;
}
