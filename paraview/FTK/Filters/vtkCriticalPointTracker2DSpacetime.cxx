#include "vtkCriticalPointTracker2DSpacetime.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkImageData.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(vtkCriticalPointTracker2DSpacetime);

vtkCriticalPointTracker2DSpacetime::vtkCriticalPointTracker2DSpacetime()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  // SetUseGPU(false);
  // SetGaussianKernelSize(2.0);
}

vtkCriticalPointTracker2DSpacetime::~vtkCriticalPointTracker2DSpacetime()
{
}

int vtkCriticalPointTracker2DSpacetime::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int vtkCriticalPointTracker2DSpacetime::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  auto rtn = TrackCriticalPoints2DSpacetime(input, output);
  return rtn;
}

int vtkCriticalPointTracker2DSpacetime::TrackCriticalPoints2DSpacetime(vtkImageData* imageData, vtkPolyData* polyData)
{
  ftk::ndarray<double> scalar;
  scalar.from_vtk_image_data(imageData);
 
  imageData->PrintSelf(std::cerr, vtkIndent(2));

  const size_t DW = scalar.shape(0), DH = scalar.shape(1), DT = scalar.shape(2);
  // fprintf(stderr, "currentTimestep=%d, DW=%lu, DH=%lu, DT=%lu\n", 
  //     currentTimestep, DW, DH, DT);

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

  // transform to match the bounds of the input image data
  const double *bounds = imageData->GetBounds();
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Scale((bounds[1]-bounds[0])/(DW-1), (bounds[3]-bounds[2])/(DH-1), (bounds[5]-bounds[4])/(DT-1));
  transform->Translate(bounds[0], bounds[2], bounds[4]);

  vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetInputData(poly);
  transformFilter->SetTransform(transform);
  transformFilter->Update();

  // polyData->DeepCopy(poly);
  polyData->DeepCopy(transformFilter->GetOutput());

  return 1;
}
