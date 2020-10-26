//---------------------------------------------------------------------------
// paramfile
//---------------------------------------------------------------------------
#if !defined(paramfile_H)
#define paramfile_H
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

using namespace std;

//none value parameters or switches
#define PTnone  0  //return of checkparam if pname does not exist
#define PTunk   1 //just a string
#define PTminus 2 // '-'
#define PTplus  3 // '+'
#define PTBS    4 // '\'
#define PTslash 5 // '/'
#define PTminm  6 // '--'
//---------------------------------------------------------------------------
class paramfilereader {
    
public:
    paramfilereader(string fn="") {filename=fn;vnum=0;vlist=NULL;plist=NULL;ptype=NULL;valreadcount=NULL;pnum=0;};
    ~paramfilereader() {done();};
    
    int openread(string fn="",bool addparams=false);
    
    int cmdlineopenread(int n,char *cmdl[]);
    
    int getvalnum(string vname) {int idx=getvalidx(vname);if(idx>=0) return (idx/2); return idx;};
    string getstring(int n); //n is not an array index of vlist, but numer of the name/value pair
    string getvname(int n); //n is not an array index of vlist
    int getreadcount(int n); //gets the read count of parameter n
    int vcount() {return vnum;} //values
    int pcount() {return pnum;} //other parameters/options
    
    //direct manipulation of existing parameters names&values
    int renameval(string vname,string newvname);
    int changeval(string vname,string newval); //only direct string value can be changed
    
    string getstring(string vname,string def="");
    int getint(string vname,int def=0);
    double getdouble(string vname,double def=0.0);
    
    int getarraycount(string vname,char sep=',');
    int getstringarray(string vname,string *data,int maxl,char sep=',');
    int getintarray(string vname,int *data,int maxl,char sep=',');
    int getdoublearray(string vname,double *data,int maxl,char sep=',');
    int getfloatarray(string vname,float *data,int maxl,char sep=',');
    
    int getparamnum(string pname) {return getparamidx(pname);}; //to be safe..
    string getparam(int n,int &type);
    int checkparam(string pname);
    
    void done() {
        if(vlist!=NULL) {delete[] vlist;};
        if(valreadcount!=NULL) {delete[] valreadcount;};
        if(plist!=NULL) {delete[] plist;delete[] ptype;};
        vlist=NULL;vnum=0;
        valreadcount=NULL;
        plist=NULL;ptype=NULL;pnum=0;};
    
private:
    string filename;
    string *vlist, *plist;
    int *valreadcount; //keeps track of read parameter values (to identify unused parameters mostly)
    int *ptype;
    int vnum,pnum;
    
    int getvalidx(string vname); //this returns real array indices (not a user function), not the number of the name/value pair
    int getparamidx(string pname);
};



#endif // paramfile_H
