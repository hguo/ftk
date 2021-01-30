#include "ftkTornadoFlow3DSource.h"
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

vtkStandardNewMacro(ftkTornadoFlow3DSource);

ftkTornadoFlow3DSource::ftkTornadoFlow3DSource() : 
  DW(32), DH(32), DD(10), DT(20)
{
  SetNumberOfInputPorts(0);
  SetNumberOfOutputPorts(1);
}

ftkTornadoFlow3DSource::~ftkTornadoFlow3DSource()
{
}

int ftkTornadoFlow3DSource::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

int ftkTornadoFlow3DSource::RequestInformation(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector* outVec)
{
  int extent[6] = {0, DW-1, 0, DH-1, 0, DD-1};
  double cell_lengths[3] = {1.0, 1.0, 1.0}, 
         origins[3] = {0.0, 0.0, 0.0};

  vtkInformation *outInfo = outVec->GetInformationObject(0);

  // time varying data
  double TimeRange[2] = {0.0, DT - 1.0};
  std::vector<double> TimeSteps;
  for (int i = 0; i < DT; i ++) 
    TimeSteps.push_back(static_cast<double>(i));

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), TimeRange, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &TimeSteps[0], DT);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
  outInfo->Set(vtkDataObject::SPACING(), cell_lengths, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), origins, 3);

  return 1;
}

int ftkTornadoFlow3DSource::RequestData(
    vtkInformation*, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *imageData = 
    vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  double currentTime;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    currentTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  else 
    currentTime = 0;

  auto vector_field = ftk::synthetic_tornado<float>(DW, DH, DD, static_cast<int>(currentTime));
  auto imageData1 = vector_field.to_vtk_image_data(); // true);
  imageData->DeepCopy(imageData1);
  
  int extent[6] = {0, DW-1, 0, DH-1, 0, DD-1};
  double cell_lengths[3] = {1.0, 1.0, 1.0}, 
         origins[3] = {0.0, 0.0, 0.0};
  
  double TimeRange[2] = {0.0, DT - 1.0};
  std::vector<double> TimeSteps;
  for (int i = 0; i < DT; i ++)
    TimeSteps.push_back(static_cast<double>(i));
  
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), TimeRange, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &TimeSteps[0], DT);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
  outInfo->Set(vtkDataObject::SPACING(), cell_lengths, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), origins, 3);

  return 1;
}
