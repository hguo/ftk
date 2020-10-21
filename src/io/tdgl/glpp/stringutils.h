#ifndef stringutilsH
#define stringutilsH

#include <string>

using namespace std;

#ifdef __cplusplus          

// extern "C" {

#endif

int StrToInt(string s);

string IntToStr(int value);

string IntToStrF(int value,int l);

double StrToFloat(string s);

string FloatToStr(double value);

string DoubleToStr(double value);

string TrimStr(string s);

#ifdef __cplusplus          

// }

#endif


//-- end unit ----------------------------------------------------------------
#endif
