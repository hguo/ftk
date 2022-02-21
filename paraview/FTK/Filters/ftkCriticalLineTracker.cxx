#include "ftkCriticalLineTracker.h"
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

vtkStandardNewMacro(ftkCriticalLineTracker);

ftkCriticalLineTracker::ftkCriticalLineTracker()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
 
  SetComputeDegrees(false);
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);

  currentTimestep = 0;

  inputDataComponents = 0;
}

ftkCriticalLineTracker::~ftkCriticalLineTracker()
{
}

int ftkCriticalLineTracker::FillInputPortInformation(int, vtkInformation *info)
{
  // info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalLineTracker::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int ftkCriticalLineTracker::RequestData_vti(
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
    field_data.from_vtk_image_data(vti, InputVariable1);
  } else if (itype == VTK_RECTILINEAR_GRID) { // TODO
    nd = vtr->GetDataDimension();
    DW = vtr->GetDimensions()[0];
    DH = vtr->GetDimensions()[1];
    DD = vtr->GetDimensions()[2];
    field_data.from_vtk_regular_data<vtkRectilinearGrid>(vtr, InputVariable1);
  } else if (itype == VTK_STRUCTURED_GRID) { // TODO
    nd = vts->GetDataDimension();
    DW = vts->GetDimensions()[0];
    DH = vts->GetDimensions()[1];
    DD = vts->GetDimensions()[2];
    field_data.from_vtk_regular_data<vtkStructuredGrid>(vts, InputVariable1);
  } else 
    assert(false);
  
  if (currentTimestep == 0) { // first timestep
    
    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable.c_str());
    vtkSmartPointer<vtkAbstractArray> da = input->GetPointData()->GetAbstractArray(InputVariable1.c_str());
    if (!da) da = input->GetPointData()->GetAbstractArray(0);

    inputDataComponents = da->GetNumberOfComponents();

    if (nd == 3) {
      fprintf(stderr, "DW=%zu, DH=%zu, DD=%zu, components=%d\n", DW, DH, DD, inputDataComponents);
      tlpr.reset(new ftk::critical_line_tracker_3d_regular(MPI_COMM_WORLD));
      
      tlpr->set_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD})); 
      tlpr->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
      tlpr->set_input_array_partial(false);
    } else 
      assert(false);

    if (itype == VTK_IMAGE_DATA)
      tlpr->set_coords_bounds(vti);
    else if (itype == VTK_RECTILINEAR_GRID)
      tlpr->set_coords_rectilinear(vtr);
    else if (itype == VTK_STRUCTURED_GRID)
      tlpr->set_coords_explicit(vts);
    else 
      assert(false);
    
    tlpr->initialize();
  }
  
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    tlpr->push_field_data_snapshot(field_data);

    if (currentTimestep != 0)
      tlpr->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      tlpr->push_field_data_snapshot(field_data);
      tlpr->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
   
    tlpr->finalize();
    // fprintf(stderr, "FINALIZING...\n");
    
    auto &surfs = tlpr->get_traced_surfaces();
#if 0 // TODO
    auto &trajs = tlpr->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });
    
    auto poly = tlpr->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);
#endif

    tlpr->reset();

    return 1;
  }

  currentTimestep ++;
  return 1;
}

int ftkCriticalLineTracker::RequestData_vtu(
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
  
  vtkSmartPointer<vtkAbstractArray> da = input->GetPointData()->GetAbstractArray(InputVariable1.c_str());
  if (!da) da = input->GetPointData()->GetAbstractArray(0);

  assert(false); // unsupported for now
#if 0
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
      assert(false); // unsupported mesh type
      // m2u.from_vtu(input);
      // tlp.reset(new ftk::critical_line_tracker_2d_unstructured(MPI_COMM_WORLD, m2u));
      // tlp->set_enable_computing_degrees( ComputeDegrees );
    } else {
      fprintf(stderr, "treating the input as a 3D grid..\n");
      assert(false); // unsupported mesh type
    }

    // vtkSmartPointer<vtkDataArray> da = input->GetPointData()->GetArray(InputVariable1.c_str());
    
    tlp->initialize();
  }
 
  ftk::ndarray<double> field_data;
  // std::cerr << "InputVariable: " << InputVariable << std::endl;
  // input->PrintSelf(std::cerr, vtkIndent(2));

  field_data.from_vtk_array(da); // input, InputVariable1);

  if (currentTimestep < inInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() )) {
    // fprintf(stderr, "currentTimestep=%d\n", currentTimestep);
    tlp->push_vector_field_snapshot(field_data);

    if (currentTimestep != 0)
      tlp->advance_timestep();

    request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
  } else { // the last timestep
    if (nt == 0) { // the only timestp
      tlp->push_vector_field_snapshot(field_data);
      tlp->update_timestep();
    }

    request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );
    currentTimestep = 0;
    
    tlp->finalize();

    auto &trajs = tlp->get_traced_critical_points();
    trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
    });

    // auto poly = tracker->get_traced_critical_points_vtk();
    auto poly = tlp->get_traced_critical_points().to_vtp({});
    output->DeepCopy(poly);

    tlp->reset();

    return 1;
  }
#endif
  currentTimestep ++;
  return 1; 
}

int ftkCriticalLineTracker::ProcessRequest(
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
