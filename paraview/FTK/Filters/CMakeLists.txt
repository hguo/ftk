set(classes
  vtkCriticalPointTracker2D
  vtkCriticalPointTracker2DSpaceTime
  vtkSpiralWoven2DSpacetimeSource
  vtkABCFlow3DVectorFieldSource
  vtkDoubleGyreVectorField2DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPointTracker2D.xml
  CriticalPointTracker2DSpaceTime.xml
  SpiralWoven2DSpacetimeSource.xml
  ABCFlow3DVectorFieldSource.xml
  DoubleGyreVectorField2DSource.xml
)
