#include "vtkCriticalPoint2DTracker.h"
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

vtkStandardNewMacro(vtkCriticalPoint2DTracker);

vtkCriticalPoint2DTracker::vtkCriticalPoint2DTracker()
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
  SetUseGPU(false);
  SetGaussianKernelSize(2.0);
}

vtkCriticalPoint2DTracker::~vtkCriticalPoint2DTracker()
{
}

void vtkCriticalPoint2DTracker::SetUseGPU(bool b)
{
  bUseGPU = b;
}

void vtkCriticalPoint2DTracker::SetGaussianKernelSize(double t)
{
  dGaussianKernelSize = t;
}

int vtkCriticalPoint2DTracker::FillOutputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int vtkCriticalPoint2DTracker::RequestData(
    vtkInformation*, 
    vtkInformationVector** inputVector, 
    vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  return TrackCriticalPoints2D(input, output);
}

int vtkCriticalPoint2DTracker::TrackCriticalPoints2D(vtkImageData* imageData, vtkPolyData* polyData)
{
#if 0
  // TODO: check compatability
  vtkSmartPointer<vtkDataArray> dataArrayRho, dataArrayPhi, dataArrayRe, dataArrayIm;
  vtkSmartPointer<vtkDataArray> dataArrayB, dataArrayPBC, dataArrayJxext, dataArrayKx, dataArrayV;
  int index;

  dataArrayRho = imageData->GetPointData()->GetArray("rho", index);
  dataArrayPhi = imageData->GetPointData()->GetArray("phi", index);
  dataArrayRe = imageData->GetPointData()->GetArray("re", index);
  dataArrayIm = imageData->GetPointData()->GetArray("im", index);
  dataArrayB = imageData->GetFieldData()->GetArray("B", index);
  dataArrayPBC = imageData->GetFieldData()->GetArray("pbc", index);
  dataArrayJxext = imageData->GetFieldData()->GetArray("Jxext", index);
  dataArrayKx = imageData->GetFieldData()->GetArray("Kx", index);
  dataArrayV = imageData->GetFieldData()->GetArray("V", index);

  // fprintf(stderr, "dataType=%d\n", dataArrayRho->GetDataType());

  GLHeader h;
  double origins[3], cell_lengths[3];

  imageData->GetDimensions(h.dims);
  if (h.dims[2] == 1) h.ndims = 2;
  else h.ndims = 3;

  imageData->GetOrigin(origins);
  imageData->GetSpacing(cell_lengths);
  for (int i=0; i<h.ndims; i++) {
    h.origins[i] = origins[i];
    h.cell_lengths[i] = cell_lengths[i];
    h.lengths[i] = h.cell_lengths[i] * h.dims[i];
  }

  double B[3], pbc1[3];
  dataArrayB->GetTuple(0, B);
  dataArrayPBC->GetTuple(0, pbc1);
  for (int i=0; i<3; i++) {
    // h.pbc[i] = (pbc1[i]>0);
    h.pbc[i] = 0; 
    h.B[i] = B[i];
  }

  h.Jxext = dataArrayJxext->GetTuple1(0);
  h.Kex = dataArrayKx->GetTuple1(0);
  h.V = dataArrayV->GetTuple1(0);

  // fprintf(stderr, "B={%f, %f, %f}, pbc={%d, %d, %d}, Jxext=%f, Kx=%f, V=%f\n", 
  //     h.B[0], h.B[1], h.B[2], h.pbc[0], h.pbc[1], h.pbc[2], h.Jxext, h.Kex, h.V);

  const int count = h.dims[0]*h.dims[1]*h.dims[2];
  float *rho = (float*)dataArrayRho->GetVoidPointer(0), 
        *phi = (float*)dataArrayPhi->GetVoidPointer(0), 
        *re = (float*)dataArrayRe->GetVoidPointer(0), 
        *im = (float*)dataArrayIm->GetVoidPointer(0);

  // build data
  GLGPUDataset *ds;
  if (h.ndims == 2) ds = new GLGPU2DDataset;
  else ds = new GLGPU3DDataset;
  ds->BuildDataFromArray(h, rho, phi, re, im); // FIXME

  if (h.ndims == 3 && iMeshType == 1) {
    GLGPU3DDataset *ds3 = (GLGPU3DDataset*)ds;
    ds3->SetMeshType(GLGPU3D_MESH_TET);
  }
  ds->BuildMeshGraph();

  VortexExtractor *ex = new VortexExtractor;
  ex->SetDataset(ds);
  ex->SetArchive(false);
#if WITH_CUDA
  ex->SetGPU(bUseGPU); // FIXME: failure fallback
#endif
  ex->SetCond(true); // TODO
  ex->SetGaugeTransformation(true); 
  ex->SetExtentThreshold(dExtentThreshold);
  ex->ExtractFaces(0);
  ex->TraceOverSpace(0);

  std::vector<VortexLine> vlines = ex->GetVortexLines(0);
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();
  vtkSmartPointer<vtkFloatArray> conditionNumberArray = vtkFloatArray::New();
  // conditionNumberArray->SetNumberOfComponents(1);

  std::vector<float> conditionNumbers;

  bool hasCond = false; 
  if (vlines.size() > 0 && vlines[0].cond.size() > 0)
    hasCond = true;

  if (h.ndims == 3) { // 3D poly lines
    std::vector<int> vertCounts;
    for (int i=0; i<vlines.size(); i++) {
      int vertCount = 0;
      const int nv = vlines[i].size()/3;
      // if (nv<2) continue;
      double p0[3];
      for (int j=0; j<nv; j++) {
        double p[3] = {vlines[i][j*3], vlines[i][j*3+1], vlines[i][j*3+2]};
        points->InsertNextPoint(p);

        if (hasCond) {
          // fprintf(stderr, "cond1=%f\n", vlines[i].cond[j]);
          conditionNumbers.push_back(vlines[i].cond[j]);
        }

        double delta[3] = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};
        double dist = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

        if (j>0 && dist>5) { // FIXME
          vertCounts.push_back(vertCount);
          vertCount = 0;
        }
        memcpy(p0, p, sizeof(double)*3);
        vertCount ++;
      }

      if (vertCount > 0) 
        vertCounts.push_back(vertCount);
    }

    int nv = 0;
    for (int i=0; i<vertCounts.size(); i++) {
      // fprintf(stderr, "vertCount=%d\n", vertCounts[i]);
      vtkSmartPointer<vtkPolyLine> polyLine = vtkPolyLine::New();
      polyLine->GetPointIds()->SetNumberOfIds(vertCounts[i]);
      for (int j=0; j<vertCounts[i]; j++)
        polyLine->GetPointIds()->SetId(j, j+nv);

      cells->InsertNextCell(polyLine);
      nv += vertCounts[i];
    }
    
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    if (hasCond) {
      conditionNumberArray->SetNumberOfValues(conditionNumbers.size());
      for (int i=0; i<conditionNumbers.size(); i++) 
        conditionNumberArray->SetValue(i, conditionNumbers[i]);
      polyData->GetPointData()->SetScalars(conditionNumberArray);
    }
  } else { // 2D data
    for (int i=0; i<vlines.size(); i++) {
      double p[3] = {vlines[i][0], vlines[i][1], vlines[i][2]};
      points->InsertNextPoint(p);
      cells->InsertNextCell(1);
      cells->InsertCellPoint(i);
    }
    polyData->SetPoints(points);
    polyData->SetVerts(cells);
  }

  delete ex;
  delete ds;
#endif
  return 1;
}
