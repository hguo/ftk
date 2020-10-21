//---------------------------------------------------------------------------
// paramfile.cpp
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"

#include "paramfile.h"

//---------------------------------------------------------------------------

int paramfilereader::openread(string fn,bool addparams)
{
 int a,b,vs,i,bsize;
 fHandle f;
 char c;
 char *buffer;
 string str,vname;
 if(fn!="") filename=fn;
 if(filename=="") return -1;
 f=FileOpen(filename,fmOpenRead);
 if(f!=-1)
     {
	  if(!addparams) done(); //remove all existing values first
      bsize = FileSize(f);
	  buffer=new char[bsize+1];
      FileRead(f,buffer,bsize);
      FileClose(f);
      buffer[bsize]=0;
      vs=0;
      for(i=0;i<bsize;i++) {c=buffer[i];if((c==0x0D) || (c==0x0A)) buffer[i]=0;if(c=='=') vs++;}
	  a=0;
	  if(addparams)
       {
		   string *tmp=new string[vnum+vnum+vs+vs];
           int *tmpvrc=new int[vnum+vs];
		   for(i=0;i<vnum+vnum;i++) tmp[i]=vlist[i];
           for(i=0;i<vnum;i++) tmpvrc[i]=valreadcount[i];
           for(i=vnum;i<vnum+vs;i++) tmpvrc[i]=0;
		   delete[] vlist;
           delete[] valreadcount;
		   vlist=tmp;
           valreadcount=tmpvrc;
	   }
	  else {vlist=new string[vs+vs];valreadcount=new int[vs];for(i=0;i<vs;i++) valreadcount[i]=0;}

      while(a<bsize)
	   {
        if(buffer[a]!=0)
         {str=&buffer[a];
          while((buffer[a]!=0) && (a<bsize)) a++;
          if((b=str.find('#'))!=string::npos) str.erase(b,str.length()-b); //remove comments
          if((b=str.find('='))!=string::npos) //paramter split
		   {
			vname=TrimStr(str.substr(0,b));
			i=getvalidx(vname);
			if(i>=0) vlist[i+1]=TrimStr(str.substr(b+1,str.length()-b-1));
			else
			 {
			  vlist[vnum+vnum]=vname;
			  vlist[vnum+vnum+1]=TrimStr(str.substr(b+1,str.length()-b-1));
			  vnum++;
			 }
		   }
         }
        else a++;
       }
      delete[] buffer;
     }
 else return -1;

 return vnum;
}
//---------------------------------------------------------------------------

int paramfilereader::cmdlineopenread(int n,char *cmdl[])
{
    if(n<2) return -2;
    if((n==2) && FileExists(cmdl[1])) return openread(cmdl[1]);
    
    int i,a,b,vs,ps;
    string str,vname;
    
    vs=0;ps=0;
    for(i=1;i<n;i++) {str=cmdl[i];if(str.find('=')!=string::npos) vs++;else ps++;}
    done();
    if(vs>0)
    {
        vlist=new string[vs+vs];
        valreadcount=new int[vs];
        for(i=1;i<n;i++)
        {
            str=cmdl[i];
            if((b=str.find('='))!=string::npos)
            {
                vname=TrimStr(str.substr(0,b));
                a=getvalidx(vname);
                if(a>=0) vlist[a+1]=TrimStr(str.substr(b+1,str.length()-b-1));
                else
                {
                    vlist[vnum+vnum]=vname;
                    vlist[vnum+vnum+1]=TrimStr(str.substr(b+1,str.length()-b-1));
                    valreadcount[vnum]=0;
                    vnum++;
                }
            }
        }
    }
    if(ps>0)
    {
        plist=new string[ps];
        ptype=new int[ps];
        for(i=1;i<n;i++)   //cmdl[0] is disregarded
        {
            str=cmdl[i];
            if((b=str.find('='))==string::npos)
            {
                switch (str[0]) {
                    case '-':
                        if(str[1]=='-') {ptype[pnum]=PTminm;plist[pnum]=str.substr(2,str.length()-2);}
                        else {ptype[pnum]=PTminus;plist[pnum]=str.substr(1,str.length()-1);}
                        break;
                    case '+':ptype[pnum]=PTplus;plist[pnum]=str.substr(1,str.length()-1);break;
                    case '\\':ptype[pnum]=PTBS;plist[pnum]=str.substr(1,str.length()-1);break;
                    case '/':ptype[pnum]=PTslash;plist[pnum]=str.substr(1,str.length()-1);break;
                    default:
                        ptype[pnum]=PTunk;plist[pnum]=str;
                }
                pnum++;
            }
        }
    }
    
    str=getstring("--readparams"); //with this value-option a parameter file is read additionally
    if((str!="") && FileExists(str)) openread(str,true);
    
    return vnum;
};
//---------------------------------------------------------------------------

string paramfilereader::getstring(int n)
{
  if((n<0) || (n>=vnum)) return "";
  valreadcount[n]++;
  return vlist[n+n+1];
};

string paramfilereader::getvname(int n)
{
  if((n<0) || (n>=vnum)) return "";
  return vlist[n+n];
};

int paramfilereader::getreadcount(int n) {
    if((n<0) || (n>=vnum)) return 0;
    return valreadcount[n];
};

//---------------------------------------------------------------------------

int paramfilereader::renameval(string vname,string newvname) {
    int idx=getvalidx(vname);
    if(idx>=0) vlist[idx]=newvname;
    return idx;
};

int paramfilereader::changeval(string vname,string newval) {
    int idx=getvalidx(vname);
    if(idx>=0) {vlist[idx+1]=newval;valreadcount[idx/2]=0;}
    return idx;
};

//---------------------------------------------------------------------------
//private !!!
int paramfilereader::getvalidx(string vname)
{
    for(int i=0;i<vnum+vnum;i+=2)
    {
        if(vname==vlist[i]) return i;
    }
    return -2;
};

string paramfilereader::getstring(string vname,string def)
{
    int idx=getvalidx(vname);
    if(idx>=0) {valreadcount[idx/2]++;return vlist[idx+1];}
    return def;
};

int paramfilereader::getint(string vname,int def)
{
    int idx=getvalidx(vname);
    if(idx>=0)  {valreadcount[idx/2]++;return StrToInt(vlist[idx+1]);}
    return def;
};

double paramfilereader::getdouble(string vname,double def)
{
    int idx=getvalidx(vname);
    if(idx>=0)  {valreadcount[idx/2]++;return StrToFloat(vlist[idx+1]);}
    return def;
};

//---------------------------------------------------------------------------

//returns the number of read doubles
int paramfilereader::getarraycount(string vname,char sep)
{
 int n,p;
 string s=getstring(vname);
 if ((s.length()==0) || (sep=='.')) return 0;
 n=1;
 for(p=0;p<s.length();p++) if(s[p]==sep) n++;
 return n;
};

//returns the number of read doubles
int paramfilereader::getstringarray(string vname,string *data,int maxl,char sep)
{
 int n,p;
 string s=getstring(vname);
 if ((s=="") || (sep=='.')) return 0;
 n=0;
 while(((p=s.find(sep,0))!=string::npos) && (n<maxl))
 {
  data[n]=s.substr(0,p);
  s.erase(0,p+1);
  n++;
 }
 if(n<maxl) {data[n]=s;n++;}
 return n;
};

//returns the number of read integers
int paramfilereader::getintarray(string vname,int *data,int maxl,char sep)
{
 int n,p;
 string s=getstring(vname);
 if (s=="") return 0;
 n=0;
 while(((p=s.find(sep,0))!=string::npos) && (n<maxl))
 {
  data[n]=StrToInt(s.substr(0,p));
  s.erase(0,p+1);
  n++;
 }
 if((s!="") && (n<maxl)) {data[n]=StrToInt(s);n++;}
 return n;
};

//returns the number of read doubles
int paramfilereader::getdoublearray(string vname,double *data,int maxl,char sep)
{
 int n,p;
 string s=getstring(vname);
 if ((s=="") || (sep=='.')) return 0;
 n=0;
 while(((p=s.find(sep,0))!=string::npos) && (n<maxl))
 {
  data[n]=StrToFloat(s.substr(0,p));
  s.erase(0,p+1);
  n++;
 }
 if((s!="") && (n<maxl)) {data[n]=StrToFloat(s);n++;}
 return n;
};

//returns the number of read floats
int paramfilereader::getfloatarray(string vname,float *data,int maxl,char sep)
{
    int n,p;
    string s=getstring(vname);
    if ((s=="") || (sep=='.')) return 0;
    n=0;
    while(((p=s.find(sep,0))!=string::npos) && (n<maxl))
    {
        data[n]=StrToFloat(s.substr(0,p));
        s.erase(0,p+1);
        n++;
    }
    if((s!="") && (n<maxl)) {data[n]=StrToFloat(s);n++;}
    return n;
};

//---------------------------------------------------------------------------

int paramfilereader::getparamidx(string pname)
{
 for(int i=0;i<pnum;i++)
 {
  if(pname==plist[i]) return i;
 }
 return -1;
};

string paramfilereader::getparam(int n,int &type)
{
 type=PTnone;
 if((n<0) || (n>=pnum)) return "";
 type=ptype[n];
 return plist[n];
};

int paramfilereader::checkparam(string pname)
{
 int idx=getparamidx(pname);
 if(idx==-1) return PTnone;
 return ptype[idx];
};
//---------------------------------------------------------------------------

