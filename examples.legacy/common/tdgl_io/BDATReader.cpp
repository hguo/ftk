#include "BDATReader.h"
#include <cassert>
#include <cstdio>
#include <limits.h>

static const char BDATTypeNames[] = {
  'b', 'B', 'h', 'H', 'i', 'I', 'q', 'Q', 'f', 'd', 's'};
static const int BDATTypeSizes[] = {
  1, 1, 2, 2, 4, 4, 8, 8, 8, 16, 1};
static const std::string BDATTypeNamesLong[] = {
  "INT8", "UINT8", "INT16", "UINT16", "INT32", "UINT32", 
  "INT64", "UINT64", "FLOAT", "DOUBLE", "CHAR"};
static const std::string BDATSignature = "BDAT";

static unsigned int TypeID2RecType(unsigned int id)
{
  switch (id) {
  case 0x100: return BDAT_INT8; 
  case 0x101: return BDAT_UINT8;
  case 0x200: return BDAT_INT16;
  case 0x201: return BDAT_UINT16;
  case 0x400: return BDAT_INT32;
  case 0x401: return BDAT_UINT32;
  case 0x402: return BDAT_FLOAT;
  case 0x800: return BDAT_INT64;
  case 0x801: return BDAT_UINT64;
  case 0x802: return BDAT_DOUBLE;
  default: return BDAT_CHAR;
  }
}

static std::string RecType2String(unsigned int recType)
{
  if (recType<=BDAT_CHAR)
    return BDATTypeNamesLong[recType];
  else 
    return std::string();
}

BDATReader::BDATReader(const std::string& filename) 
  : state(BDAT_STATE_HEADER)
{
  valid = true;
  fp = fopen(filename.c_str(), "rb");

  if (!fp) {
    valid = false;
    return;
  } 
  
  std::string signature; 
  ReadString(4, &signature);
  if (signature != BDATSignature) {
    valid = false;
  }

  unsigned int BOM; 
  Read(BDAT_UINT32, &BOM);
}

BDATReader::~BDATReader()
{
  if (fp) 
    fclose(fp);
}

std::string BDATReader::ReadNextRecordInfo()
{
  assert(state == BDAT_STATE_HEADER);
  if (!Valid()) return std::string();

  unsigned int IDlen;
  unsigned short length;
 
  if (!Read(BDAT_UINT32, &IDlen)) 
    return std::string();
  length = IDlen & 0xff;
  recID = IDlen >> 8;

  ReadString(length, &recName);

  unsigned int typeID;
  Read(BDAT_UINT32, &typeID);
  Read(BDAT_UINT32, &recNum); 
  Read(BDAT_UINT32, &recLen);

  recType = TypeID2RecType(typeID);

  // fprintf(stderr, "recID=%d, recName=%s, recType=%s, recNum=%d, recLen=%d\n", 
  //     recID, recName.c_str(), RecType2String(recType).c_str(), recNum, recLen);

  state = BDAT_STATE_DATA;

  return recName;
}

void BDATReader::ReadNextRecordData(std::string *buf)
{
  assert(state == BDAT_STATE_DATA);
  if (!Valid()) return;

  ReadString(recLen*recNum, buf);
  state = BDAT_STATE_HEADER;
}

bool BDATReader::Read(int typeName, void *val)
{
  if (!Valid()) return false;

  int typeSize = BDATTypeSizes[typeName];
  size_t count = fread(val, typeSize, 1, fp);

  return count>0;
}

bool BDATReader::ReadString(int length, std::string *str)
{
  if (!Valid()) return false;
  
  str->resize(length);
  size_t count = fread((char*)str->data(), 1, length, fp);
  
  return count>0;
}
