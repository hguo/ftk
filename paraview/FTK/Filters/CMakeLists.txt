set(classes
  vtkCriticalPoint2DTracker
  vtkSpiral2DSource
  vtkABCFlowSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE ftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPoint2DTracker.xml
  Spiral2DSource.xml
  ABCFlowSource.xml
)
