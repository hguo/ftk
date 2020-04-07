set(classes
  vtkCriticalPoint2DTracker
  vtkSpiralWoven2DStackSource
  vtkABCFlowSource
  vtkDoubleGyreVectorField2DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPoint2DTracker.xml
  SpiralWoven2DStackSource.xml
  ABCFlowSource.xml
  DoubleGyreVectorField2DSource.xml
)
