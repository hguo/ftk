set(classes
  vtkCriticalPointTracker2D
  vtkCriticalPointTracker2DSpacetime
  vtkCriticalPointTracker3D
  vtkSpiralWoven2DSpacetimeSource
  vtkSpiralWoven2DSource
  vtkDoubleGyre2DSource
  vtkABCFlow3DSource
  ftkTornadoFlow3DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPointTracker2D.xml
  CriticalPointTracker2DSpacetime.xml
  CriticalPointTracker3D.xml
  SpiralWoven2DSpacetimeSource.xml
  SpiralWoven2DSource.xml
  ABCFlow3DSource.xml
  TornadoFlow3DSource.xml
  DoubleGyre2DSource.xml
)
