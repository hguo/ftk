#ifndef _FTK_WRITE_VTP_HH
#define _FTK_WRITE_VTP_HH

#include <ftk/config.hh>
#include <ftk/utils/string.hh>
#include <ftk/external/diy/mpi.hpp>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>

#ifndef FTK_HAVE_PARAVIEW
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#endif

namespace ftk {

inline void write_polydata(
    const std::string& filename, 
    vtkSmartPointer<vtkPolyData> poly, 
    std::string format="auto", // format can be vtp, ply...
    diy::mpi::communicator comm = MPI_COMM_WORLD) 
{
  if (format == "auto" || format.empty()) { // determine format by extension
    if (ends_with(filename, "ply")) format = "ply";
    else if (ends_with(filename, "stl")) format = "stl";
    else if (ends_with(filename, "vtk")) format = "legacy";
    else if (ends_with(filename, "pvtp")) format = "pvtp";
    else format = "vtp";
  }

  if (format == "pvtp") {
    vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkXMLPPolyDataWriter::New();
    pwriter->EncodeAppendedDataOff();
    pwriter->SetFileName(filename.c_str());
    pwriter->SetNumberOfPieces(comm.size());
    pwriter->SetStartPiece(0);
    pwriter->SetEndPiece(comm.size()-1);
    pwriter->SetInputData(poly);
    pwriter->Write();
    
    comm.barrier();
    const auto prefix = remove_file_extension(filename);
    const auto f = prefix + "_" + std::to_string(comm.rank()) + ".vtp";
    write_polydata(f, poly, "vtp", MPI_COMM_SELF);
  } else if (format == "vtp") {
    if (comm.rank() == 0) {
      fprintf(stderr, "writing vtp to %s\n", filename.c_str());
      vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
      writer->SetFileName(filename.c_str());
      writer->SetInputData(poly);
      writer->Write();
    }
  } else if (format == "ply") {
#ifndef FTK_HAVE_PARAVIEW
    vtkSmartPointer<vtkPLYWriter> writer = vtkPLYWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
#endif
  } else if (format == "stl") {
#ifndef FTK_HAVE_PARAVIEW
    vtkSmartPointer<vtkSTLWriter> writer = vtkSTLWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
#endif
  } else if (format == "legacy") {
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
  }
}

}

#endif // FTK_HAVE_VTK
#endif
