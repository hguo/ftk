#include <hypermesh/regular_simplex_mesh.hh>
#include <ftk/basic/distributed_union_find.hh>

#ifndef DIM
  #define DIM 3
#endif


// ==========================================================

struct intersection_t {
  template <class Archive> void serialize(Archive & ar) {
    ar(eid, x[0], x[1], x[2], val);
  }

  // std::vector<float> x; // the spacetime coordinates of the trajectory
  float x[DIM];

  // std::vector<float> corner; // the spacetime coordinates of the left corner of the element

  float val; // scalar value at the intersection
  
  std::string eid; // element id
  std::set<std::string> related_elements; 
};


// ==========================================================

// template<int DIM>
struct point_t{

  int&  operator[](unsigned i)                          { return corner[i]; }
  int   operator[](unsigned i) const                    { return corner[i]; }

  int corner[DIM]; // the spacetime coordinates of the left corner of the element
  // std::vector<int> corner; 

  // point_t(int DIM): corner(DIM) {}
};


// ==========================================================


struct Message_Feature : public Message_Union_Find {
  Message_Feature() {

  }

  void send_intersection(const std::pair<std::string, intersection_t>& pair) {
    tag = "intersection"; 

    strs.push_back(pair.first); 

    auto& I = pair.second; 
    
    strs.push_back(I.eid); 
    strs.push_back(std::to_string(I.val)); 
    
    // strs.push_back(std::to_string(I.x[0])); 
    // strs.push_back(std::to_string(I.x[1])); 
    // strs.push_back(std::to_string(I.x[2])); 

    // for(int i = 0; i < I.x.size(); ++i) {
    //   strs.push_back(std::to_string(I.x[i])); 
    // }

    for(int i = 0; i < DIM; ++i) {
      strs.push_back(std::to_string(I.x[i])); 
    }
  }

  void receive_intersection(std::pair<std::string, intersection_t>& pair) {
    pair.first = strs[0]; 
    
    intersection_t& I = pair.second; 
    
    I.eid = strs[1];
    I.val = std::stof(strs[2]); 

    // I.x[0] = std::stof(strs[3]); 
    // I.x[1] = std::stof(strs[4]); 
    // I.x[2] = std::stof(strs[5]); 

    for(int i = 3; i < strs.size(); ++i) {
      // I.x.push_back(std::stof(strs[i])); 
      I.x[i-3] = std::stof(strs[i]); 
    }
  }

};



// ==========================================================
// DIY serialization



namespace diy
{
  template<>
  struct Serialization<intersection_t>
  {
      static void save(BinaryBuffer& bb, const intersection_t& msg)
      {
          diy::save(bb, msg.x);
          // diy::save(bb, msg.corner);
          diy::save(bb, msg.val);
          diy::save(bb, msg.eid);
          diy::save(bb, msg.related_elements);
      }

      static void load(BinaryBuffer& bb, intersection_t& msg)
      {
          diy::load(bb, msg.x);
          // diy::load(bb, msg.corner);
          diy::load(bb, msg.val);
          diy::load(bb, msg.eid);
          diy::load(bb, msg.related_elements);
      }
  };
}


namespace diy
{
  template<>
  struct Serialization<point_t>
  {
      static void save(BinaryBuffer& bb, const point_t& msg)
      {
          diy::save(bb, msg.corner);
      }

      static void load(BinaryBuffer& bb, point_t& msg)
      {
          diy::load(bb, msg.corner);
      }
  };
}



namespace diy
{
  template<>
  struct Serialization<Message_Feature>
  {
      static void save(BinaryBuffer& bb, const Message_Feature& msg)
      {
          diy::save(bb, msg.tag);
          // diy::save(bb, msg.str);
          diy::save(bb, msg.strs);
      }

      static void load(BinaryBuffer& bb, Message_Feature& msg)
      {
          diy::load(bb, msg.tag);
          // diy::load(bb, msg.str);
          diy::load(bb, msg.strs);
      }
  };
}
