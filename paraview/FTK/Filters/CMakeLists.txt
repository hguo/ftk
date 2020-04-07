set(classes
  vtkCriticalPoint2DTracker
  vtkSpiralWoven2DStackSource
  vtkABCFlow3DVectorFieldSource
  vtkDoubleGyreVectorField2DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPoint2DTracker.xml
  SpiralWoven2DStackSource.xml
  ABCFlow3DVectorFieldSource.xml
  DoubleGyreVectorField2DSource.xml
)
