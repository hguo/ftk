#include "ftkSpiralWoven2DSource.h"
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

vtkStandardNewMacro(ftkSpiralWoven2DSource);

ftkSpiralWoven2DSource::ftkSpiralWoven2DSource() 
  : DW(32), DH(32), DT(10), 
    ScalingFactor(15.0), StartTime(0.0), TimeScale(1.0), 
    NoiseInjection(0.0)
{
  SetNumberOfInputPorts(0);
  SetNumberOfOutputPorts(1);
}

ftkSpiralWoven2DSource::~ftkSpiralWoven2DSource()
{
}

int ftkSpiralWoven2DSource::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

int ftkSpiralWoven2DSource::RequestInformation(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector* outVec)
{
  int extent[6] = {0, DW-1, 0, DH-1, 0, 0};
  double cell_lengths[3] = {1.0, 1.0, 1.0}, 
         origins[3] = {0.0, 0.0, 0.0};

  vtkInformation *outInfo = outVec->GetInformationObject(0);

  // time varying data
  double timeRange[2] = {0.0, DT - 1.0};
  std::vector<double> timeSteps;
  for (int i = 0; i < DT; i ++)
    timeSteps.push_back(static_cast<double>(i) * TimeScale + StartTime);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0], DT);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
  outInfo->Set(vtkDataObject::SPACING(), cell_lengths, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), origins, 3);

  return 1;
}

int ftkSpiralWoven2DSource::RequestData(
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
    currentTime = StartTime;

  auto scalar = ftk::synthetic_woven_2D<float>(DW, DH, currentTime+1e-4, ScalingFactor);
  if (NoiseInjection > 0) {
    scalar.perturb(NoiseInjection);
  }
  auto imageData1 = scalar.to_vtk_image_data();
  imageData->ShallowCopy(imageData1);
  
  int extent[6] = {0, DW-1, 0, DH-1, 0, 0};
  double cell_lengths[3] = {1.0, 1.0, 1.0}, 
         origins[3] = {0.0, 0.0, 0.0};

  double timeRange[2] = {0.0, DT - 1.0};
  std::vector<double> timeSteps;
  for (int i = 0; i < DT; i ++)
    timeSteps.push_back(static_cast<double>(i) * TimeScale + StartTime);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0], DT);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
  outInfo->Set(vtkDataObject::SPACING(), cell_lengths, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), origins, 3);

  return 1;
}
