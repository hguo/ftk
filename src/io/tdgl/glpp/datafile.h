//---------------------------------------------------------------------------
// datafile
//---------------------------------------------------------------------------
#if !defined(datafile_H)
#define datafile_H
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#include "fileutils.h"

using namespace std;

//data types (identified by 32bit number)
//second byte is the typically the size in bytes of primitive data types, but not used
//lowest bit = 1 means unsigned
//second lowest = 1 means floating point
//third bit = 1 means complex
//(bits 4-7 length in bits)
#define DT_bit     0x0010  //single bits, the length is in bytes though, if number of relvant bits in an array are needed, use another datablock
#define DT_byte    0x0100  //char
#define DT_ubyte   0x0101  //unsigned char (also identifies an ascii string or UTF-8)

#define DT_word    0x0200 //short int
#define DT_uword   0x0201 //unsigned short int (also UTF-16)

#define DT_int     0x0400 //int
#define DT_uint    0x0401 //unsigned int (also UTF-32)

#define DT_long    0x0800 //long int, long long
#define DT_ulong   0x0801 //unsigned long int

#define DT_single  0x0402 //float, 4 byte
#define DT_double  0x0802 //double, 8 byte
#define DT_quad    0x1002 //quad precision floating point "long double", __float128

//if complex numbers are represented as re+i*im or as amplitude and phase, needs to be determined by user
#define DT_csingle 0x0806 //complex float
#define DT_cdouble 0x1006 //complex double
#define DT_cquad   0x2006 //complex "long double"

#define DT_user    0xFFFFFF00 //any record size is possible


#define TMPSIZE    67108864 //64MB


//---------------------------------------------------------------------------
/*
BDAT (binary data) file format:
- 4 bytes: ID = "BDAT"
- 4 bytes: endian BOM
- <data blocks>
 
 
data block:
- DT_uint: form 0xIIIIIINN, where NN is length of name (<=255) and 0xIIIIII = ID (should be used if NN=0)
- NN bytes: data name (should be used if ID=0)
- DT_uint: type (see above)
- DT_uint: records = N
- DT_uint: record size in bytes=sb [should be 1 for bit type as well, but the actual number of bits needs to be specified in another datablock]
- N*sb bytes: data
 
 */
//---------------------------------------------------------------------------

typedef struct {
    string name;
    unsigned int ID; //only lowest 24 bit are valid
    unsigned int type,records,recordsize; //byte size=records*recordsize
} datablockinfo;


class BDATreader {
  public:
    BDATreader(string fn="") {filename=fn;pos=0;state=0;endian=0;hasTOC=false;TOCN=0;};
    ~BDATreader() {close();};

    
    int open(string fn="",bool makeTOC=false);
    void close() {if(state>0) FileClose(f);state=0;deleteTOC();};
    int createTOC();
    void deleteTOC();
    int getTOClength() {if(hasTOC) return TOCN;return 0;};
    
    //output file data content/structure, restores current file state
    int printlist();
    
    //state needs to be 1 or 3, set state to 2
    int getnextblockinfo(datablockinfo &db);
    //after blockinfo, state needs to be 2, info in lastinfo
    int skipblock(); //moves only the file pointer and sets state to 3
    int readblock(unsigned char *buffer,unsigned int maxlength=0); //buffer is assumed to be records*recordsize bytes, otherwise use maxlength
    
    //same special data type implmentations (which might take care of endian conversion in future)
    //maxlength is the number of records, but byte length
    int readblock(char *buffer,unsigned int maxlength=0);
    int readblock(int *buffer,unsigned int maxlength=0);
    int readblock(unsigned int *buffer,unsigned int maxlength=0);
    int readblock(float *buffer,unsigned int maxlength=0);
    int readblock(double *buffer,unsigned int maxlength=0);
    
    //read and convert
    int readconvert(int *buffer,unsigned int maxlength=0);    //converts read data into type int
    int readconvert(float *buffer,unsigned int maxlength=0);  //converts read data into type float
    int readconvert(double *buffer,unsigned int maxlength=0); //converts read data into type double
    
    //print warnings if system does not support "standard"/expected primitive sizes
    void primitives_warning();
    
    
    //seek functions needing a TOC (will be initialized if used the first time), state and pos will be changed, also executes getnextblockinfo
    int seek(string name,datablockinfo &db);
    int seek(unsigned int ID,datablockinfo &db);
    
  private:
    long pos; //position in file
    int state; //current reading state, 0-file closed, 1-file opened and header recognized, ready to read data, 2 - block info read, 3 - data read, 4 - eof, <0 read error
    int endian;
    long size; //filesize
    string filename;
    fHandle f;
    datablockinfo lastinfo;
    
    //TOC
    bool hasTOC;
    int TOCN;
    long *TOCpos;
    unsigned int *TOCID;
    string *TOCnames;
};

class BDATwriter {
public:
    BDATwriter(string fn="") {filename=fn;state=0;};
    ~BDATwriter() {close();};
    
    
    int open(string fn="",bool append=false);
    void close() {if(state>0) FileClose(f);state=0;};
    
    //state needs to be 1
    int writeblock(string name,unsigned int ID,unsigned int type,unsigned int records,unsigned int recordsize,char *buffer); //this is the main function, all other write functions use it
    int writeblock(datablockinfo &db,char *buffer);
    
    //specific formats
    int write(string name,unsigned int ID,string str);
    int write(string name,unsigned int ID,int *buffer,unsigned int records=1);
    int write(string name,unsigned int ID,unsigned int *buffer,unsigned int records=1);
    int write(string name,unsigned int ID,long *buffer,unsigned int records=1);
    int write(string name,unsigned int ID,unsigned long *buffer,unsigned int records=1);
    int write(string name,unsigned int ID,float *buffer,unsigned int records=1);
    int write(string name,unsigned int ID,double *buffer,unsigned int records=1);
    
    //special complex writers, notice one needs to type cast to float or double, but records is interpreted as complex number count, i.e., half the actual FP numbers
    int splitwrite(string namere,string nameim,unsigned int IDre,unsigned int IDim,float *buffer,unsigned int records=1);
    int splitwrite(string namere,string nameim,unsigned int IDre,unsigned int IDim,double *buffer,unsigned int records=1);

private:
    int state; //current file state, 0-file closed, 1-file opened for writing
    long size; //filesize
    string filename;
    fHandle f;
};

#endif // datafile_H
