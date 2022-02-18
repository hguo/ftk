set(classes
  ftkCriticalPointTracker
  # ftkCriticalPointTracker2D
  # ftkCriticalPointTracker2DUnstructured
  # ftkCriticalPointTracker3D
  # ftkThresholdTracker2D
  ftkLevelsetTracker2D
  ftkLevelsetTracker3D
  #ftkSpiralWoven2DSpacetimeSource
  ftkSpiralWoven2DSource
  ftkCappedWovenGradient2DSource
  ftkDoubleGyre2DSource
  ftkMerger2DSource
  ftkMovingExtremum2DSource
  ftkMovingExtremum3DSource
  ftkMovingRamp3DSource
  ftkMovingDualRamp3DSource
  ftkABCFlow3DSource
  ftkTornadoFlow3DSource
)

vtk_module_add_module(FTKFilters
  CLASSES ${classes})
vtk_module_link(FTKFilters
  PRIVATE libftk)

paraview_add_server_manager_xmls(XMLS 
  CriticalPointTracker.xml
  # CriticalPointTracker2D.xml
  # CriticalPointTracker2DUnstructured.xml
  # CriticalPointTracker3D.xml
  # ThresholdTracker2D.xml
  LevelsetTracker2D.xml
  LevelsetTracker3D.xml
  #SpiralWoven2DSpacetimeSource.xml
  SpiralWoven2DSource.xml
  CappedWovenGradient2DSource.xml
  Merger2DSource.xml
  MovingExtremum2DSource.xml
  MovingExtremum3DSource.xml
  MovingRamp3DSource.xml
  MovingDualRamp3DSource.xml
  ABCFlow3DSource.xml
  TornadoFlow3DSource.xml
  DoubleGyre2DSource.xml
)
