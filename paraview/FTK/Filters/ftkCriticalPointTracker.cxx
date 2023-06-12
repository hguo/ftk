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
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh_implicit.hh>

vtkStandardNewMacro(ftkCriticalPointTracker);

ftkCriticalPointTracker::ftkCriticalPointTracker()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
 
  SetComputeDegrees(false);
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);

  currentTimestep = 0;

  inputDataComponents = 0;
}

ftkCriticalPointTracker::~ftkCriticalPointTracker()
{
}

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

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
  vtkImageData *vti = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkRectilinearGrid *vtr = vtkRectilinearGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkStructuredGrid *vts = vtkStructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  const int nt = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  // const double *timesteps = inInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  
  const int itype = input->GetDataObjectType();
  int nd;
  size_t DW, DH, DD;
  ftk::ndarray<double> field_data;
  
  if (itype == VTK_IMAGE_DATA) {
    nd = vti->GetDataDimension();
    DW = vti->GetDimensions()[0];
    DH = vti->GetDimensions()[1];
    DD = vti->GetDimensions()[2];
    field_data.from_vtk_image_data(vti, InputVariable);
  } else if (itype == VTK_RECTILINEAR_GRID) {
    nd = vtr->GetDataDimension();
    DW = vtr->GetDimensions()[0];
    DH = vtr->GetDimensions()[1];
    DD = vtr->GetDimensions()[2];
    field_data.from_vtk_regular_data<vtkRectilinearGrid>(vtr, InputVariable);
  } else if (itype == VTK_STRUCTURED_GRID) {
    nd = vts->GetDataDimension();
    int dims[3];
    vts->GetDimensions(dims);
    DW = dims[0];
    DH = dims[1]; 
    DD = dims[2];
    field_data.from_vtk_regular_data<vtkStructuredGrid>(vts, InputVariable);
  } else 
    assert(false);
  
  if (currentTimestep == 0) { // first timestep
    
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
      
      tcpr->set_enable_computing_degrees( ComputeDegrees );
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

    if (itype == VTK_IMAGE_DATA)
      tcpr->set_coords_bounds(vti);
    else if (itype == VTK_RECTILINEAR_GRID)
      tcpr->set_coords_rectilinear(vtr);
    else if (itype == VTK_STRUCTURED_GRID)
      tcpr->set_coords_explicit(vts);
    else 
      assert(false);
    
    tcpr->initialize();
  }
  
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

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
    // fprintf(stderr, "FINALIZING...\n");
    auto &trajs = tcpr->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });
    
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
  using namespace ftk;
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
    std::set<unsigned char> cellTypesSet;
    vtkSmartPointer<vtkCellTypes> cellTypes = vtkNew<vtkCellTypes>(); 
    input->GetCellTypes(cellTypes);
    for (vtkIdType i = 0; i < cellTypes->GetNumberOfTypes(); i ++) {
      cellTypesSet.insert( cellTypes->GetCellType( i ) );
    }

    if (cellTypesSet.find( VTK_TETRA ) != cellTypesSet.end()) { // the input is regarded as a simplicial 3D mesh
      // TODO
    } else if (cellTypesSet.find( VTK_TRIANGLE) != cellTypesSet.end()) { // the input is regarded as a simplicial 2D mesh
      fprintf(stderr, "treating the input as a 2D grid..\n");

      std::shared_ptr<simplicial_unstructured_2d_mesh<>> m2(new simplicial_unstructured_2d_mesh<>);
      m2->from_vtu(input);

      std::shared_ptr<simplicial_unstructured_extruded_2d_mesh_implicit<>> m3i(new simplicial_unstructured_extruded_2d_mesh_implicit<>(m2));
      auto m3 = std::dynamic_pointer_cast<simplicial_unstructured_extruded_2d_mesh<>>(m3i);

      tcp.reset(new critical_point_tracker_2d_unstructured(MPI_COMM_WORLD, m3));
      tcp->set_enable_computing_degrees( ComputeDegrees );
    } else {
      fprintf(stderr, "treating the input as a 3D grid..\n");
      assert(false); // unsupported mesh type
    }

    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable.c_str());
    
    tcp->initialize();
  }
 
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  field_data.from_vtk_array(da); // input, InputVariable);

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    tcp->push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tcp->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      tcp->push_vector_field_snapshot(field_data);
      tcp->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tcp->finalize();

    auto &trajs = tcp->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });

    // auto poly = tracker->get_traced_critical_points_vtk();
    auto poly = tcp->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);

    tcp->reset();

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
