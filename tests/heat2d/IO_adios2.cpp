/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ADIOS2.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <iostream>
#include <string>

#include <adios2.h>

adios2::ADIOS *ad = nullptr;
adios2::Engine writer;
adios2::Variable<double> varT;
adios2::Variable<unsigned int> varGndx;

IO::IO(const Settings &s, MPI_Comm comm)
{
    ad = new adios2::ADIOS(s.configfile, comm, adios2::DebugON);

    adios2::IO io = ad->DeclareIO("SimulationOutput");
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default writer
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
        io.AddTransport("File", {{"Library", "POSIX"}});
    }

    if (!s.rank)
    {
//        std::cout << "Using " << io.m_EngineType << " engine for output" << std::endl;
    }

    // define T as 2D global array
    varT = io.DefineVariable<double>(
        "T",
        // Global dimensions
        {s.gndx, s.gndy},
        // starting offset of the local array in the global space
        {s.offsx, s.offsy},
        // local size, could be defined later using SetSelection()
        {s.ndx, s.ndy});

    io.DefineAttribute<std::string>("description", 
            "Temperature from simulation", "T");
    io.DefineAttribute<std::string>("unit", 
            "C", varT.Name());

    // homogeneous last coordinate 2D -> 3D
    const std::string extent = "0 " + std::to_string(s.gndx + 1) + " " +
                               "0 " + std::to_string(s.gndy + 1) + " " +
                               "0 1";

    const std::string imageData = R"(
        <?xml version="1.0"?>
        <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
          <ImageData WholeExtent=")" + extent + R"(" Origin="0 0 0" Spacing="1 1 1">
            <Piece Extent=")" + extent + R"(">
              <CellData Scalars="T">
                  <DataArray Name="T" />
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>)";

    io.DefineAttribute<std::string>("vtk.xml", imageData);

    writer = io.Open(s.outputfile, adios2::Mode::Write, comm);

    // Some optimization:
    // we promise here that we don't change the variables over steps
    // (the list of variables, their dimensions, and their selections)
    writer.LockWriterDefinitions();
}

IO::~IO()
{
    writer.Close();
    delete ad;
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
	//reduce memory footprint, adios will provide memory from its buffer
	if(s.span)
	{
		writer.BeginStep();
		// pre-allocate memory in adios2 buffer and provide a span
		adios2::Variable<double>::Span spanT = writer.Put<double>(varT);
		// populate the span
		ht.data_noghost(spanT.data());
		// collect span data and get min/max
		writer.EndStep();
	}
	else
	{
        writer.BeginStep();
		// using Put() you promise the pointer to the data will be intact
		// until the end of the output step.
		// We need to have the vector object here not to destruct here until the end
		// of function.
		std::vector<double> v = ht.data_noghost();
		writer.Put<double>(varT, v.data());
		writer.EndStep();
	}
}
