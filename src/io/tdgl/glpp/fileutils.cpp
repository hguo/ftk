/*
 FileUtilites for Non-VCL applications
 wraps file related functions in VCL style to use on any C/C++ compiler

 (c) Andreas Glatz, V2.0
*/
#include <stdio.h>
#include <string>
#include "fileutils.h"

using namespace std;

fHandle FileOpen(string FileName, int Mode)
{fHandle r=-1;
 FILE *f;
 switch(Mode) {

  case fmOpenRead:f=fopen(FileName.c_str(),"rb");break;
  case fmOpenWrite:f=fopen(FileName.c_str(),"r+b");break;
  case fmOpenReadWrite:f=fopen(FileName.c_str(),"r+b");break;
  case fmOpenAppend:f=fopen(FileName.c_str(),"a+b");break;
  case fmCreate:f=fopen(FileName.c_str(),"wb");break;
 }
 if(f!=NULL) r=(fHandle) f;
 //printf("Debug: fopen: f=%d, r=%d\n",(int) f,r);
 return r;
}

void FileClose(fHandle Handle)
{
 if(Handle!=-1) fclose((FILE *) Handle);
}

fHandle FileCreate(string FileName)
{
 return FileOpen(FileName,fmCreate);

}

bool FileExists(string FileName)
{
 fHandle f=FileOpen(FileName,fmOpenRead);
 if(f!=-1)
  {FileClose(f);
   return true;
  }
 return false;
} 

long FileRead(fHandle Handle, void *Buffer, long Count)
{
 return fread(Buffer,Count,1,(FILE *) Handle);
}

long FileWrite(fHandle Handle, const void *Buffer, long Count)
{
 return fwrite(Buffer,Count,1,(FILE *) Handle);
}
 
long FileSeek(fHandle Handle, long Offset, int Origin)
{
 if(Handle!=-1)
  {
   fseek((FILE *) Handle,Offset,Origin);
   return ftell((FILE *) Handle);
  }
 return -1;
}

long FileSize(fHandle Handle)
{long s=-1;
 if(Handle!=-1)
   {
    s = FileSeek(Handle,0,2);
    FileSeek(Handle,0,0);
   }
 return s;
}
void FileFlush(fHandle Handle)
{
 if(Handle!=-1)
   {
    fflush((FILE *) Handle);
   }
}

void DeleteFile(string FileName)
{
  remove(FileName.c_str());
}
