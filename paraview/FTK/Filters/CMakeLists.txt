set(classes
  vtkCriticalPoint2DTracker
  vtkSpiral2DStackSource
  vtkABCFlowSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPoint2DTracker.xml
  Spiral2DStackSource.xml
  ABCFlowSource.xml
)
