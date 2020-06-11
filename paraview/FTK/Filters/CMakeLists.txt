set(classes
  vtkCriticalPointTracker2D
  vtkCriticalPointTracker2DSpacetime
  vtkSpiralWoven2DSpacetimeSource
  vtkSpiralWoven2DSource
  vtkABCFlow3DSource
  vtkDoubleGyre2DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPointTracker2D.xml
  CriticalPointTracker2DSpacetime.xml
  SpiralWoven2DSpacetimeSource.xml
  SpiralWoven2DSource.xml
  ABCFlow3DSource.xml
  DoubleGyre2DSource.xml
)
