#include "ftkLevelsetTracker3D.h"
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

vtkStandardNewMacro(ftkLevelsetTracker3D);

ftkLevelsetTracker3D::ftkLevelsetTracker3D()
  : tracker(diy::mpi::communicator()), 
  UseGPU(false)
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  SetUseGPU(false);

  currentTimestep = 0;

  inputDataComponents = 0;
}

ftkLevelsetTracker3D::~ftkLevelsetTracker3D()
{
}

int ftkLevelsetTracker3D::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  return 1;
}

int ftkLevelsetTracker3D::RequestUpdateExtent(
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


int ftkLevelsetTracker3D::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkLevelsetTracker3D::RequestSlicedSurfaces(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  std::vector<double> timesteps(nt);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timesteps[0]);

  double rt; // requested time
  int ts; // actual time step
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    rt = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    ts = std::find_if(timesteps.begin(), timesteps.end(), 
            [rt](double t) { return std::abs(t - rt) < 1e-6; })
          - timesteps.begin();
    if (ts == timesteps.size() - 1) ts = ts - 1; // TODO FIXME: the tracker does not slice the last timestep
  } else {
    rt = timesteps[0];
    ts = 0;
  }
  
  auto poly = tracker.get_isovolume().slice_time(ts).to_vtp();
  vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
  normalGenerator->SetInputData(poly);
  normalGenerator->ConsistencyOn();
  normalGenerator->ComputePointNormalsOff();
  normalGenerator->ComputeCellNormalsOn();
  // normalGenerator->SetFlipNormals(true);
  // normalGenerator->AutoOrientNormalsOn();
  normalGenerator->Update();

  // output->DeepCopy( tracker.get_isovolume().slice_time(ts).to_vtp() );
  output->DeepCopy( normalGenerator->GetOutput() );
  return 1;
}

int ftkLevelsetTracker3D::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  if (std::abs( Threshold - tracker.get_threshold() ) > 1e-6)
    tracker_needs_recompute = true;

  if (tracker_needs_recompute) {
    const size_t DW = input->GetDimensions()[0], 
                 DH = input->GetDimensions()[1],
                 DD = input->GetDimensions()[1];
    const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
   
    fprintf(stderr, "requesting data.\n");
    if (currentTimestep == 0) { // first timestep
      inputDataComponents = input->GetNumberOfScalarComponents();
      
      fprintf(stderr, "DW=%zu, DH=%zu, DD=%zu, components=%d\n", DW, DH, DD, inputDataComponents);
      tracker.reset();
      tracker.set_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD})); 
      tracker.set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
      tracker.set_input_array_partial(false);
      tracker.set_threshold( Threshold );
      tracker.initialize();
    }
    
    ftk::ndarray<double> field_data;
    // std::cerr << "InputVariable: " << InputVariable << std::endl;
    // input->PrintSelf(std::cerr, vtkIndent(2));

    field_data.from_vtk_image_data(input, InputVariable);

    if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
      tracker.push_field_data_snapshot(field_data);
      if (currentTimestep != 0)
        tracker.advance_timestep();

      request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
    } else { // the last timestep
      if (nt == 0) { // the only timestp
        tracker.push_field_data_snapshot(field_data);
        tracker.update_timestep();
      }

      request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
      currentTimestep = 0;
      
      tracker.finalize();
      tracker_needs_recompute = false;
      
      return RequestSlicedSurfaces(request, inputVector, outputVector);
    }

    currentTimestep ++;
    return 1; 
  } else {
    return RequestSlicedSurfaces(request, inputVector, outputVector);
  }
}
