/*
 FileUtilites for Non-VCL applications
 wraps file related functions in VCL style to use on any C/C++ compiler

 (c) Andreas Glatz, V2.0
*/
#ifndef fileutilsH
#define fileutilsH

#include <string>

using namespace std;

#define fmOpenRead 0x0
#define fmOpenWrite 0x1
#define fmOpenReadWrite 0x2
#define fmOpenAppend 0x3
#define fmCreate 0xFF

//pointer and length/position size
#define FileUtilsBits 8*sizeof(long)

typedef long fHandle;

#ifdef __cplusplus          

extern "C" {

#endif
             
fHandle FileOpen(string FileName, int Mode);

void FileClose(fHandle Handle);

fHandle FileCreate(string FileName);

bool FileExists(string FileName);

long FileRead(fHandle Handle, void *Buffer, long Count);

long FileWrite(fHandle Handle, const void *Buffer, long Count);
 
long FileSeek(fHandle Handle, long Offset, int Origin);
   
long FileSize(fHandle Handle);

void FileFlush(fHandle Handle);

void DeleteFile(string FileName);

#ifdef __cplusplus          

}

#endif

/* Origin
0	Der Dateizeiger wird Offset Bytes nach dem Dateianfang positioniert. 
1	Der Dateizeiger wird Offset Bytes nach der aktuellen Position positioniert. 
2	Der Dateizeiger wird Offset Bytes vor dem Dateiende positioniert.

Bei Erfolg gibt FileSeek die neue Position des Dateizeigers zurueck. Schlaegt die Funktion fehl, wird -1 zurueckgegeben.
*/
//-- end unit ----------------------------------------------------------------
#endif
