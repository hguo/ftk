#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "stringutils.h"

using namespace std;


int StrToInt(string s)
{
    if((s.length()>1) && (s[1]=='x')) return strtol(s.c_str(),NULL,16);
    return strtol(s.c_str(),NULL,10);
}

string IntToStr(int value)
{char s[16];
 snprintf(s,16,"%d",value);
 return string(s);
}

string IntToStr(long value)
{char s[32];
 snprintf(s,32,"%ld",value);
 return string(s);
}

string IntToStrF(int value,int l)
{char s[32];
 string fmt="%0"+IntToStr(l)+"d";
 snprintf(s,32,fmt.c_str(),value);
 return string(s);
}

string FloatToStr(double value)
{char s[32];
 snprintf(s,32,"%le",value);
 return string(s);
}

string DoubleToStr(double value)
{char s[32];
 snprintf(s,32,"%.16le",value);
 return string(s);
}

double StrToFloat(string s)
{
 return atof(s.c_str());	
}

string TrimStr(string s)
{int a,b,l;
 l=s.length();
 a=0;while((a<l) && (s[a]==' ')) a++;	
 b=l-1;while((b>=0) && (s[b]==' ')) b--;
 if(a<=b) return s.substr(a,b-a+1);
 return "";
}

