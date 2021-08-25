#include "ftkCriticalPointTracker.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkAbstractArray.h"
#include "vtkCellArray.h"
#include "vtkImageData.h"
#include "vtkSphereSource.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

vtkStandardNewMacro(ftkCriticalPointTracker);

ftkCriticalPointTracker::ftkCriticalPointTracker() :
    UseGPU(false),
    GaussianKernelSize(1.0)
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);

  currentTimestep = 0;

  inputDataComponents = 0;
}

ftkCriticalPointTracker::~ftkCriticalPointTracker()
{
}

#if 0
int ftkCriticalPointTracker::RequestInformation(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  return 1;
}

int ftkCriticalPointTracker::RequestUpdateExtent(
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
#endif

int ftkCriticalPointTracker::FillInputPortInformation(int, vtkInformation *info)
{
  // info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalPointTracker::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalPointTracker::RequestData_vti(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  fprintf(stderr, "requesting vti data.\n");

  // TODO FIXME check if input is 2D or 3D
  
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  
  const int itype = input->GetDataObjectType();
  if (currentTimestep == 0) { // first timestep
    int nd;
    size_t DW, DH, DD;
    
    if (itype == VTK_IMAGE_DATA) {
      vtkImageData *vti = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
      nd = vti->GetDataDimension();
      DW = vti->GetDimensions()[0];
      DH = vti->GetDimensions()[1];
      DD = vti->GetDimensions()[2];
    }
    
    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable.c_str());
    vtkSmartPointer<vtkAbstractArray> da = input->GetPointData()->GetAbstractArray(InputVariable.c_str());
    if (!da) da = input->GetPointData()->GetAbstractArray(0);

    inputDataComponents = da->GetNumberOfComponents();

    if (nd == 2) {
      fprintf(stderr, "DW=%zu, DH=%zu, components=%d\n", DW, DH, inputDataComponents);
      tcpr.reset(new ftk::critical_point_tracker_2d_regular(MPI_COMM_WORLD));
      if (inputDataComponents == 1) { // scalar field
        tcpr->set_domain(ftk::lattice({2, 2}, {DW-3, DH-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
        tcpr->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
        tcpr->set_input_array_partial(false);
        tcpr->set_scalar_field_source(ftk::SOURCE_GIVEN);
        tcpr->set_vector_field_source(ftk::SOURCE_DERIVED);
        tcpr->set_jacobian_field_source(ftk::SOURCE_DERIVED);
      } else if (inputDataComponents > 1) { // vector field
        tcpr->set_domain(ftk::lattice({1, 1}, {DW-2, DH-2})); // the indentation is needed becase the jacoobian field will be automatically derived
        tcpr->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
        tcpr->set_input_array_partial(false);
        tcpr->set_scalar_field_source(ftk::SOURCE_NONE);
        tcpr->set_vector_field_source(ftk::SOURCE_GIVEN);
        tcpr->set_jacobian_field_source(ftk::SOURCE_DERIVED);
      } else 
        assert(false);
    } else if (nd == 3) {
      fprintf(stderr, "DW=%zu, DH=%zu, DD=%zu, components=%d\n", DW, DH, DD, inputDataComponents);
      tcpr.reset(new ftk::critical_point_tracker_3d_regular(MPI_COMM_WORLD));
      if (inputDataComponents == 1) { // scalar field
        tcpr->set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
        tcpr->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
        tcpr->set_input_array_partial(false);
        tcpr->set_scalar_field_source(ftk::SOURCE_GIVEN);
        tcpr->set_vector_field_source(ftk::SOURCE_DERIVED);
        tcpr->set_jacobian_field_source(ftk::SOURCE_DERIVED);
      } else if (inputDataComponents > 1) { // vector field
        tcpr->set_domain(ftk::lattice({1, 1, 1}, {DW-2, DH-2, DD-2})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
        tcpr->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
        tcpr->set_input_array_partial(false);
        tcpr->set_scalar_field_source(ftk::SOURCE_NONE);
        tcpr->set_vector_field_source(ftk::SOURCE_GIVEN);
        tcpr->set_jacobian_field_source(ftk::SOURCE_DERIVED);
      } else 
        assert(false);
    } else 
      assert(false);
    
    tcpr->initialize();
  }
  
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  if (itype == VTK_IMAGE_DATA) {
    vtkImageData *vti = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    field_data.from_vtk_image_data(vti, InputVariable);
  } else 
    assert(false); // TODO FIXME not implemented yet

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    if (inputDataComponents == 1) tcpr->push_scalar_field_snapshot(field_data);
    else tcpr->push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tcpr->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      if (inputDataComponents == 1) tcpr->push_scalar_field_snapshot(field_data);
      else tcpr->push_vector_field_snapshot(field_data);
      tcpr->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tcpr->finalize();
    auto poly = tcpr->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);

    tcpr->reset();

    return 1;
  }

  currentTimestep ++;
  return 1;
}

int ftkCriticalPointTracker::RequestData_vtu(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  fprintf(stderr, "requesting vtu data.\n");
  
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  
  vtkSmartPointer<vtkAbstractArray> da = input->GetPointData()->GetAbstractArray(InputVariable.c_str());
  if (!da) da = input->GetPointData()->GetAbstractArray(0);
 
  if (currentTimestep == 0) { // first timestep
    m2u.from_vtu(input);

    tcp2du.reset(new ftk::critical_point_tracker_2d_unstructured(MPI_COMM_WORLD, m2u));
    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable.c_str());
    
    tcp2du->initialize();
  }
 
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  field_data.from_vtk_array(da); // input, InputVariable);

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    tcp2du->push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tcp2du->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      tcp2du->push_vector_field_snapshot(field_data);
      tcp2du->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tcp2du->finalize();

    auto &trajs = tcp2du->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });

    // auto poly = tracker->get_traced_critical_points_vtk();
    auto poly = tcp2du->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);

    tcp2du->reset();

    return 1;
  }

  currentTimestep ++;
  return 1; 
}

int ftkCriticalPointTracker::ProcessRequest(
    vtkInformation* request, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_INFORMATION())) {
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

    return 1;
  } else if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT())) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    double *inTimes = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
    if (inTimes) 
      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), inTimes[currentTimestep]);

    return 1;
  } else if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_DATA())) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    // vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    // vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
   
    const int itype = input->GetDataObjectType();
    // fprintf(stderr, "itype=%d\n", itype);

    if (itype == VTK_IMAGE_DATA || itype == VTK_STRUCTURED_GRID || itype == VTK_RECTILINEAR_GRID)
      return RequestData_vti(request, inputVector, outputVector);
    else if (itype == VTK_UNSTRUCTURED_GRID)
      return RequestData_vtu(request, inputVector, outputVector);
    else 
      assert(false);

    return 1;
  } else 
    return 1;
}
