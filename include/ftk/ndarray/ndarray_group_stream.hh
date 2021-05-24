#ifndef _FTK_NDARRAY_GROUP_STREAM_HH
#define _FTK_NDARRAY_GROUP_STREAM_HH

#include <ftk/ndarray/ndarray_group.hh>

namespace ftk {

// json config
struct ndarray_group_stream {
  // WIP: possible JSON specifications:
  //  - key/value.  The array name (key) and specification of the array (value)
  // required fields for value: 
  //    - type, string.  Must be one of "synthetic" and "file".  This field may be ommited
  //      if `format' is given.
  //    - name (required if type is synthetic), string.  Must be one of the follows: "woven", 
  //      "double_gyre", "merger_2d".
  //   optional fields:
  //    - filenames (required if type is file), string.  The list of filenames will be 
  //      determined by glob(3)
  //    - mesh_filename (required for unstructured mesh, must be in vtu format in the 
  //      current version), string.
  //    - format (required if type is file and format is float32/float64), string.  If not 
  //      given, the format will be determined by the filename extension.  The value of this 
  //      field must be one of the follows: vti, nc, h5, float32, float64.
  //    - variables (required if format is nc/h5, optional for vti), array of strings.
  //      - the number of components is the length of the array.
  //      - if not given, the defaulat value is ["scalar"]
  //    - components (to be determined), array of number of components per variable
  //    - dimensions (required if format is floaot32/float64), array of integers, e.g.  
  //      [width, height, depth]
  //    - n_timesteps, integer.  The default is 32 for synthetic data; the number can be 
  //      automatically derived from file; the number can be automatically derived from files
  //    - perturbation, number.  Add gaussian perturbation to the data
};

}

#endif
