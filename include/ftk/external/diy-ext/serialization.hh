#ifndef _DIY_EXT
#define _DIY_EXT

#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy/storage.hpp>
#include <cassert>
#include <cstring>

namespace diy {
  struct StringBuffer : public BinaryBuffer {
    std::string &str;
    size_t pos; 

    explicit StringBuffer(std::string& str_, size_t pos_=0) : str(str_), pos(pos_) {}
    void clear() {str.clear(); pos = 0;}
    void reset() {pos = 0;}

    operator bool() const {return pos < str.size();}

    inline void save_binary(const char *x, size_t count) {
      if (pos + count > str.size()) str.resize(pos + count);
      memcpy((char*)(str.data()+pos), x, count);
      pos += count;
    }

    inline void append_binary(const char *x, size_t count) {
      // TODO: check with tom to see what should append_binary do for binary buffer.
    }

    inline void load_binary(char *x, size_t count) {
      memcpy(x, str.data()+pos, count);
      pos += count;
    }

    inline void load_binary_back(char *x, size_t count) {
      memcpy(x, str.data()+str.size()-count, count);
    }
  };

  //////////
  template <typename T> void serialize(const T& obj, std::string& buf)
  {
    buf.clear();
    diy::StringBuffer bb(buf);
    diy::save(bb, obj);
  }

  template <typename T> void serializeToFile(const T& obj, const std::string& filename)
  {
    FILE *fp = fopen(filename.c_str(), "wb");
    assert(fp);
    diy::detail::FileBuffer bb(fp);
    diy::save(bb, obj);
    fclose(fp);
  }

  template <typename T> void unserialize(std::string& buf, T& obj)
  {
    diy::StringBuffer bb(buf);
    diy::load(bb, obj);
    buf.clear();
  }

  template <typename T> void unserialize(const std::string& str, T& obj)
  {
    std::string buf(str);
    unserialize(buf, obj);
  }

  template <typename T> void unserializeFromFile(const std::string& filename, T& obj)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    assert(fp);
    diy::detail::FileBuffer bb(fp);
    diy::load(bb, obj);
    fclose(fp);
  }
}

#endif
