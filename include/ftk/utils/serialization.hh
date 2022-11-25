#ifndef _DIY_EXT
#define _DIY_EXT

#include <ftk/config.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy/storage.hpp>
#include <cassert>
#include <cstring>

#if FTK_HAVE_TBB
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_unordered_set.h>
#endif

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
      str.append(x, count);
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
  template <typename T> void serializeToString(const T& obj, std::string& buf)
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
  
  template <typename T1, typename T2> void serializeToFile(const T1& obj1, const T2& obj2, const std::string& filename)
  {
    FILE *fp = fopen(filename.c_str(), "wb");
    assert(fp);
    diy::detail::FileBuffer bb(fp);
    diy::save(bb, obj1);
    diy::save(bb, obj2);
    fclose(fp);
  }

  template <typename T> void unserializeFromString(const std::string& buf, T& obj)
  {
    std::string buf1(buf);
    diy::StringBuffer bb(buf1);
    diy::load(bb, obj);
  }

  template <typename T> bool unserializeFromFile(const std::string& filename, T& obj)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) return false;
    diy::detail::FileBuffer bb(fp);
    diy::load(bb, obj);
    fclose(fp);
    return true;
  }
  
  template <typename T1, typename T2> bool unserializeFromFile(const std::string& filename, T1& obj1, T2 &obj2)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) return false;
    diy::detail::FileBuffer bb(fp);
    diy::load(bb, obj1);
    diy::load(bb, obj2);
    fclose(fp);
    return true;
  }

#if FTK_HAVE_TBB
  template <typename Key, typename T, typename HashCompare, typename A>
  struct Serialization<tbb::concurrent_hash_map<Key, T, HashCompare, A>> {
    typedef tbb::concurrent_hash_map<Key, T, HashCompare, A> hash_map;
    
    static void save(BinaryBuffer& bb, const hash_map& m) {
      size_t s = m.size();
      diy::save(bb, s);
      for (typename hash_map::const_iterator it = m.begin(); it != m.end(); it ++)
        diy::save(bb, *it);
    }

    static void load(BinaryBuffer& bb, hash_map& m) {
      size_t s;
      diy::load(bb, s);
      for (size_t i = 0; i < s; i ++) {
        Key k;
        T v;
        diy::load(bb, k);
        diy::load(bb, v);
        m.insert(std::make_pair(k, v));
      }
    }
  };

  template <typename Key,
          typename Hasher, // = tbb::tbb_hash<Key>,
          typename Equality, //  = std::equal_to<Key>,
          typename Allocator> //  = tbb::tbb_allocator<Key>
  struct Serialization<tbb::concurrent_unordered_set<Key, Hasher, Equality, Allocator>> {
    typedef tbb::concurrent_unordered_set<Key, Hasher, Equality, Allocator> unordered_set;
    
    static void save(BinaryBuffer& bb, const unordered_set& m) {
      size_t s = m.size();
      diy::save(bb, s);
      for (const auto &k : m)
        diy::save(bb, k);
    }

    static void load(BinaryBuffer& bb, unordered_set& m) {
      size_t s;
      diy::load(bb, s);
      for (size_t i = 0; i < s; i ++) {
        Key k;
        diy::load(bb, k);
        m.insert(k);
      }
    }

  };
#endif
}

#endif
