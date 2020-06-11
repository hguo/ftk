#include "vtkCriticalPointTracker2D.h"
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
#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(vtkCriticalPointTracker2D);

vtkCriticalPointTracker2D::vtkCriticalPointTracker2D()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);

  currentTimestep = 0;
}

vtkCriticalPointTracker2D::~vtkCriticalPointTracker2D()
{
}

void vtkCriticalPointTracker2D::SetUseGPU(bool b)
{
  bUseGPU = b;
}

void vtkCriticalPointTracker2D::SetGaussianKernelSize(double t)
{
  dGaussianKernelSize = t;
}

int vtkCriticalPointTracker2D::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int vtkCriticalPointTracker2D::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() ))
    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  else 
    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );

  auto rtn = TrackCriticalPoints2DSpaceTime(input, output);
  currentTimestep ++;

  return rtn;
}

int vtkCriticalPointTracker2D::TrackCriticalPoints2DSpaceTime(vtkImageData* imageData, vtkPolyData* polyData)
{
  ftk::ndarray<double> scalar;
  scalar.from_vtk_image_data(imageData);
  
  const size_t DW = scalar.shape(0), DH = scalar.shape(1), DT = scalar.shape(2);
  fprintf(stderr, "currentTimestep=%d, DW=%lu, DH=%lu, DT=%lu\n", 
      currentTimestep, DW, DH, DT);

  return 1;

  ftk::critical_point_tracker_2d_regular tracker; 
  tracker.set_domain(ftk::lattice({2, 2}, {DW-3, DH-3}));
  // tracker.set_domain(ftk::lattice({4, 4}, {DW-6, DH-6}));
  tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
  tracker.initialize();

  tracker.push_scalar_field_spacetime(scalar);
  while (tracker.advance_timestep()) {}

  tracker.finalize();
  auto poly = tracker.get_traced_critical_points_vtk();

  polyData->DeepCopy(poly);

  return 1;
}
