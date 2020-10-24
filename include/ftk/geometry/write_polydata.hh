#ifndef _FTK_WRITE_VTP_HH
#define _FTK_WRITE_VTP_HH

#include <ftk/ftk_config.hh>
#include <ftk/utils/string.hh>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>

namespace ftk {

inline void write_polydata(
    const std::string& filename, 
    vtkSmartPointer<vtkPolyData> poly, 
    std::string format="auto") // format can be vtp, ply
{
  if (format == "auto" || format.empty()) { // determine format by extension
    if (ends_with(filename, "ply")) format = "ply";
    else if (ends_with(filename, "stl")) format = "stl";
    else if (ends_with(filename, "vtk")) format = "legacy";
    else format = "vtp";
  }

  if (format == "vtp") {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkXMLPolyDataWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
  } else if (format == "ply") {
    vtkSmartPointer<vtkPLYWriter> writer = vtkPLYWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
  } else if (format == "stl") {
    vtkSmartPointer<vtkSTLWriter> writer = vtkSTLWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
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
