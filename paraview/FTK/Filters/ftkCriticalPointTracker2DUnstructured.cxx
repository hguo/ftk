#include "ftkCriticalPointTracker2DUnstructured.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkAbstractArray.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSphereSource.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(ftkCriticalPointTracker2DUnstructured);

ftkCriticalPointTracker2DUnstructured::ftkCriticalPointTracker2DUnstructured()
//   : tracker(diy::mpi::communicator())
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  currentTimestep = 0;
}

ftkCriticalPointTracker2DUnstructured::~ftkCriticalPointTracker2DUnstructured()
{
}

int ftkCriticalPointTracker2DUnstructured::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  return 1;
}

int ftkCriticalPointTracker2DUnstructured::RequestUpdateExtent(
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


int ftkCriticalPointTracker2DUnstructured::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalPointTracker2DUnstructured::RequestData(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  fprintf(stderr, "requesting data.\n");

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  
  vtkSmartPointer<vtkAbstractArray> da = input->GetPointData()->GetAbstractArray(InputVariable.c_str());
  if (!da) da = input->GetPointData()->GetAbstractArray(0);
 
  if (currentTimestep == 0) { // first timestep
    m2.from_vtu(input);

    tracker.reset(new ftk::critical_point_tracker_2d_unstructured(diy::mpi::communicator(), m2));
    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable.c_str());

#if 0 
    inputDataComponents = da->GetNumberOfComponents();
    const size_t DW = input->GetDimensions()[0], 
                 DH = input->GetDimensions()[1];
    
    fprintf(stderr, "DW=%zu, DH=%zu, components=%d\n", DW, DH, inputDataComponents);
    if (inputDataComponents == 1) { // scalar field
      tracker->set_domain(ftk::lattice({2, 2}, {DW-3, DH-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
      tracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
      tracker->set_input_array_partial(false);
      tracker->set_scalar_field_source(ftk::SOURCE_GIVEN);
      tracker->set_vector_field_source(ftk::SOURCE_DERIVED);
      tracker->set_jacobian_field_source(ftk::SOURCE_DERIVED);
    } else if (inputDataComponents > 1) { // vector field
      tracker->set_domain(ftk::lattice({1, 1}, {DW-2, DH-2})); // the indentation is needed becase the jacoobian field will be automatically derived
      tracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
      tracker->set_input_array_partial(false);
      tracker->set_scalar_field_source(ftk::SOURCE_NONE);
      tracker->set_vector_field_source(ftk::SOURCE_GIVEN);
      tracker->set_jacobian_field_source(ftk::SOURCE_DERIVED);
    } else 
      assert(false);
#endif
    tracker->initialize();
  }
 
#if 1
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  field_data.from_vtk_array(da); // input, InputVariable);

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    tracker->push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tracker->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      tracker->push_vector_field_snapshot(field_data);
      tracker->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tracker->finalize();

#if 1
    auto &trajs = tracker->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });
#endif

    // auto poly = tracker->get_traced_critical_points_vtk();
    auto poly = tracker->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);

    tracker->reset();

    return 1;
  }
#endif

  currentTimestep ++;
  return 1; 
}
