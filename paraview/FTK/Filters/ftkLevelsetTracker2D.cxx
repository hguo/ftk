#include "ftkLevelsetTracker2D.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
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
  : tracker(diy::mpi::communicator()), 
  Threshold(0.0), 
  OutputType(0),
  ZTimeScale(1.0)
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  currentTimestep = 0;
  inputDataComponents = 0;
}

ftkLevelsetTracker2D::~ftkLevelsetTracker2D()
{
}

int ftkLevelsetTracker2D::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkLevelsetTracker2D::RequestUpdateExtent(
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

int ftkLevelsetTracker2D::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  if (OutputType == 1) { // traced outputs; remove timesteps in the outputs
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
  }

  return 1;
}

int ftkLevelsetTracker2D::RequestData(
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
                 DH = input->GetDimensions()[1];
    const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
   
    fprintf(stderr, "requesting data.\n");
    if (currentTimestep == 0) { // first timestep
      inputDataComponents = input->GetNumberOfScalarComponents();
      
      fprintf(stderr, "DW=%zu, DH=%zu, components=%d\n", DW, DH, inputDataComponents);
      tracker.reset();
      tracker.set_domain(ftk::lattice({0, 0}, {DW, DH})); 
      tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
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
    
      if (OutputType == 0)
        return RequestSliced(request, inputVector, outputVector);
      else {
        output->DeepCopy( tracker.get_surfaces().to_vtp(true, true, ZTimeScale) );
      }
    }

    currentTimestep ++;
    return 1; 
  } else {
    if (OutputType == 0)
      return RequestSliced(request, inputVector, outputVector);
    else {
      output->DeepCopy( tracker.get_surfaces().to_vtp(true, true, ZTimeScale) );
      return 1;
    }
  }
}

int ftkLevelsetTracker2D::RequestSliced(
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

  fprintf(stderr, "rt=%f, ts=%d\n", rt, ts);
  output->DeepCopy( tracker.get_surfaces().slice_time(ts).to_vtp({}));
  return 1;
}
