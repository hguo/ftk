#include "ftkCriticalPointTracker3D.h"
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
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(ftkCriticalPointTracker3D);

ftkCriticalPointTracker3D::ftkCriticalPointTracker3D()
  : tracker(diy::mpi::communicator()), 
  UseGPU(false), GaussianKernelSize(1.0)
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);

  currentTimestep = 0;

  inputDataComponents = 0;
}

ftkCriticalPointTracker3D::~ftkCriticalPointTracker3D()
{
}

int ftkCriticalPointTracker3D::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  return 1;
}

int ftkCriticalPointTracker3D::RequestUpdateExtent(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  double *inTimes = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  if (inTimes) 
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), inTimes[currentTimestep]);

  return 1;
}


int ftkCriticalPointTracker3D::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalPointTracker3D::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  fprintf(stderr, "requesting data.\n");

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  const size_t DW = input->GetDimensions()[0], 
               DH = input->GetDimensions()[1],
               DD = input->GetDimensions()[1];
  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  
  if (currentTimestep == 0) { // first timestep
    inputDataComponents = input->GetNumberOfScalarComponents();
    
    fprintf(stderr, "DW=%zu, DH=%zu, DD=%zu, components=%d\n", DW, DH, DD, inputDataComponents);
    tracker.set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    tracker.set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
    tracker.set_input_array_partial(false);
    if (inputDataComponents == 1) { // scalar field
      tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
      tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
      tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
    } else if (inputDataComponents > 1) { // vector field
      tracker.set_scalar_field_source(ftk::SOURCE_NONE);
      tracker.set_vector_field_source(ftk::SOURCE_GIVEN);
      tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
    } else 
      assert(false);
    
    tracker.initialize();
  }
  
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  field_data.from_vtk_image_data(input, InputVariable);

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    if (inputDataComponents == 1) tracker.push_scalar_field_snapshot(field_data);
    else tracker.push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tracker.advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      if (inputDataComponents == 1) tracker.push_scalar_field_snapshot(field_data);
      else tracker.push_vector_field_snapshot(field_data);
      tracker.update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tracker.finalize();
    auto poly = tracker.get_traced_critical_points_vtk();

    // transform to match the bounds of the input image data
    const double *bounds = input->GetBounds();
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale((bounds[1]-bounds[0])/(DW-1), (bounds[3]-bounds[2])/(DH-1), (bounds[5]-bounds[4])/(DD-1));
    transform->Translate(bounds[0], bounds[2], bounds[4]);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    // output->DeepCopy(poly);
    output->DeepCopy(transformFilter->GetOutput());

    tracker.reset();

    return 1;
  }

  currentTimestep ++;
  return 1; 
}
