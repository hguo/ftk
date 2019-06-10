#ifndef _BDATREADER_H
#define _BDATREADER_H

#include <string>

class BDATReader {
  enum {
    BDAT_STATE_HEADER, 
    BDAT_STATE_DATA
  };

public:
  BDATReader(const std::string &filename); 
  ~BDATReader();

  bool Valid() const {return valid;}

  std::string ReadNextRecordInfo();
  void ReadNextRecordData(std::string *buf); //!< returns recType

  unsigned int RecType() const {return recType;}
  unsigned int RedID() const {return recID;}
  
private:  
  bool Read(int typeName, void *val);
  bool ReadString(int length, std::string *str);
  
private:
  FILE *fp;

private: // the state machine
  unsigned short recID;
  std::string recName; 
  unsigned int recType, recNum, recLen; 

  int valid, state;
};

enum {
  BDAT_INT8 = 0, 
  BDAT_UINT8, 
  BDAT_INT16, 
  BDAT_UINT16, 
  BDAT_INT32, 
  BDAT_UINT32, 
  BDAT_INT64, 
  BDAT_UINT64, 
  BDAT_FLOAT, 
  BDAT_DOUBLE, 
  BDAT_CHAR
};

#endif
