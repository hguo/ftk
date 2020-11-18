//---------------------------------------------------------------------------
// graphicstools
// (c) Andreas Glatz, 2011
//---------------------------------------------------------------------------
#if !defined(graphicstools_H)
#define graphicstools_H
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#include "fileutils.h"

#define PXM_P0 0 //universal, see description below
#define PXM_P1 1 //Portable bitmap ASCII
#define PXM_P2 2 //Portable graymap ASCII
#define PXM_P3 3 //Portable pixmap ASCII
#define PXM_P4 4 //Portable bitmap binary
#define PXM_P5 5 //Portable graymap binary
#define PXM_P6 6 //Portable pixmap binary


#define MAXGRADS 1024
#define GRAD_RAINBOW 0 //default built-in gradient

#define GRAD_USER_LOADED 0x10000 //ID flag to identify gradient as user loaded (needs deallocayion)
#define GRAD_USER_SET 0x20000 //ID flag to identify gradient as user loaded (no deallocation, array needs to stay valid)




using namespace std;

typedef struct { //RGB or HSV where H is also normalized to 1 (and not the usual 360 deg or 2pi)
 union { double R,H; };
 union { double G,S; };
 union { double B,V; };
 double A;  //alpha channel
} dcolor;

typedef struct { //RGB or HSV where H is also normalized to 1 (and not the usual 360 deg or 2pi)
    union { double R,H; };
    union { double G,S; };
    union { double B,V; };
    double A;  //alpha channel
    double X1; //(float) eXtra info 1, like position within a gradient
    int X2;    //(int) eXtra info 2, like (gnu) interpolation function 
} Xcolor;

typedef struct {
    unsigned int ID;
    int N;
    int type; //1: RGBA, 2: HSVA, 3: CMYK
    Xcolor *xcols;
} gradient;


static const char * ctypes[]={"BW","RGB","GRAY","RGBA","HSV","HSVA","CMYK"};
static const char * PXMext[]={".pxm",".pbm",".pgm",".ppm",".pbm",".pgm",".ppm"};
//---------------------------------------------------------------------------
//built-in RGB gradients 
static Xcolor CGrainbow[] ={{0.471412,0.108766,0.527016,0,0,3},{0.31106,0.11758,0.664469,0,0.0625,3},{0.250728,0.225386,0.769152,0,0.125,3},{0.24408,0.361242,0.816084,0,0.1875,3},{0.266122,0.486664,0.802529,0,0.25,3},{0.305919,0.585575,0.739666,0,0.3125,3},{0.36048,0.655759,0.645692,0,0.375,3},{0.429842,0.701849,0.540321,0,0.4375,3},{0.513417,0.72992,0.440682,0,0.5,3},{0.607651,0.743718,0.358588,0,0.5625,3},{0.705038,0.742591,0.299167,0,0.625,3},{0.794549,0.721158,0.260829,0,0.6875,3},{0.863512,0.670771,0.236564,0,0.75,3},{0.901014,0.582826,0.216542,0,0.8125,3},{0.902853,0.453964,0.192014,0,0.875,3},{0.878107,0.293208,0.160481,0,0.9375,3},{0.857359,0.131106,0.132128,0,1.,3}}; //len=17
static Xcolor CGrainbow2[] ={{0,0,0,0,0,3},{0.75,0,0.5,0,0.1,3},{0,0,0.5,0,0.4,3},{0,0.625,0,0,0.5,3},{1,1,0,0,0.7,3},{1,0.125,0,0,0.9,3},{1,1,1,0,1.0,3}}; //7
static Xcolor CGrainbow3[] ={{1,0,0,0,0,3},{1,1,0,0,0.16667,3},{0,1,0,0,0.33333,3},{0,1,1,0,0.5,3},{0,0,1,0,0.66667,3},{1,0,1,0,0.833333,3},{1,0,0,0,1.0,3}};//7
static Xcolor CGdarkrainbow[] ={{0.237736,0.340215,0.575113,0,0,3},{0.253651,0.344893,0.558151,0,0.1,3},{0.264425,0.423024,0.3849,0,0.2,3},{0.291469,0.47717,0.271411,0,0.3,3},{0.416394,0.555345,0.24182,0,0.4,3},{0.624866,0.673302,0.264296,0,0.5,3},{0.813033,0.766292,0.303458,0,0.6,3},{0.877875,0.731045,0.326896,0,0.7,3},{0.812807,0.518694,0.303459,0,0.8,3},{0.72987,0.239399,0.230961,0,0.9,3},{0.72987,0.239399,0.230961,0,1.,3}};//11
static Xcolor CGtemperature[] ={{0.164,0.043,0.85,0,0,3},{0.15,0.306,1,0,0.09090909,3},{0.25,0.63,1,0,0.18181818,3},{0.45,0.853,1,0,0.27272727,3},{0.67,0.973,1,0,0.36363636,3},{0.88,1,1,0,0.45454545,3},{1,1,0.75,0,0.54545454,3},{1,0.88,0.6,0,0.63636363,3},{1,0.679,0.45,0,0.72727272,3},{0.97,0.43,0.37,0,0.81818181,3},{0.85,0.15,0.196,0,0.9090909,3},{0.65,0,0.13,0,1,3}};//12
static Xcolor CGtemperature2[] ={{0.178927,0.305394,0.933501,0,0,3},{0.308746,0.441842,0.940894,0,0.08333,3},{0.453318,0.567063,0.950106,0,0.1667,3},{0.642359,0.720535,0.964988,0,0.25,3},{0.819984,0.859297,0.982692,0,0.3333,3},{0.935699,0.951565,0.993729,0,0.4167,3},{0.984192,0.987731,0.911643,0,0.5,3},{0.995282,0.992317,0.727853,0,0.5833,3},{0.992503,0.986373,0.425376,0,0.6667,3},{0.955963,0.863115,0.283425,0,0.75,3},{0.904227,0.657999,0.241797,0,0.8333,3},{0.858405,0.449932,0.203562,0,0.9167,3},{0.817319,0.134127,0.164218,0,1.,3}};//13
static Xcolor CGthermo[]={{0.163302,0.119982,0.79353,0,0,3},{0.254221,0.313173,0.892833,0,0.09091,3},{0.407119,0.543513,0.938275,0,0.1818,3},{0.572715,0.73338,0.95065,0,0.2727,3},{0.720374,0.855234,0.928635,0,0.3636,3},{0.831017,0.903518,0.868326,0,0.4545,3},{0.894452,0.880139,0.77279,0,0.5455,3},{0.907999,0.789417,0.652903,0,0.6364,3},{0.874505,0.639254,0.522424,0,0.7273,3},{0.79915,0.446142,0.391971,0,0.8182,3},{0.685695,0.242449,0.268261,0,0.9091,3},{0.534081,0.0853132,0.16669,0,1.,3}};//12
static Xcolor CGsolar[]={{0.468742,0.,0.0158236,0,0,3},{0.822129,0.122225,0.0039559,0,0.25,3},{0.969963,0.376081,0.0322881,0,0.5,3},{1.,0.646929,0.0801709,0,0.75,3},{1.,0.820127,0.126955,0,1.,3}};//5
static Xcolor CGsunset[]={{0.,0.,0.,0,0,3},{0.372793,0.1358,0.506503,0,0.1667,3},{0.788287,0.259816,0.270778,0,0.3333,3},{0.979377,0.451467,0.0511329,0,0.5,3},{1.,0.682688,0.129771,0,0.6667,3},{1.,0.882236,0.491094,0,0.8333,3},{1.,1.,1.,0,1.,3}};//7
static Xcolor CGneon[]={{0.720287,0.923781,0.297597,0,0,3},{0.741831,0.759184,0.2915,0,0.1111,3},{0.763375,0.594586,0.285403,0,0.2222,3},{0.784919,0.429989,0.279306,0,0.3333,3},{0.824907,0.318689,0.25001,0,0.4444,3},{0.858159,0.314389,0.27497,0,0.5556,3},{0.857705,0.331729,0.385417,0,0.6667,3},{0.840408,0.290678,0.51353,0,0.7778,3},{0.823111,0.249627,0.641643,0,0.8889,3},{0.805814,0.208576,0.769757,0,1.,3}};//10

//---------------------------------------------------------------------------


//color functions for RGBA (HSVA/HSBA) using double color channel information (type dcolor)
class colorfunction {
public:
    
    colorfunction();
    ~colorfunction();
    
    //color conversion
    dcolor RGBtoHSV(dcolor RGB);
    dcolor HSVtoRGB(dcolor HSV);
    unsigned int get32bitcolor(dcolor RGBA) {
        unsigned int r,g,b,a;
        r=(unsigned int) (256.0*RGBA.R);if(r>0xFF) r=0xFF;
        g=(unsigned int) (256.0*RGBA.G);if(g>0xFF) g=0xFF;
        b=(unsigned int) (256.0*RGBA.B);if(b>0xFF) b=0xFF;
        a=(unsigned int) (256.0*RGBA.A);if(a>0xFF) a=0xFF;
        return (b+(g<<8)+(r<<16)+(a<<24));};
    dcolor set32bitcolor(unsigned int rgba) {
        dcolor col;
        col.R=0.003921568627451*((rgba & 0x00FF0000) >> 16);
        col.G=0.003921568627451*((rgba & 0x0000FF00) >> 8);
        col.B=0.003921568627451*((rgba & 0x000000FF));
        col.A=0.003921568627451*((rgba & 0xFF000000) >> 24);
        return col;
    };
    
    unsigned int getgrayscale(dcolor RGB,unsigned int maxval=0xFF) //by default 8-bit, can be changed to 16 or 32 bit, by setting maxval
    {//uses: (0.299*r + 0.587*g + 0.114*b) [alpha channel not used - mix with background color before]
        double gray=0.299*RGB.R + 0.587*RGB.G + 0.114*RGB.B;
        unsigned int g=(unsigned int) ((maxval+1)*gray);if(g>maxval) g=maxval;
        return g;
    };
    unsigned int getgrayscale(unsigned int rgb,unsigned int maxval=0xFF) //by default 8-bit, can be changed to 16 or (32) bit, by setting maxval
    {//uses: (0.299*r/255 + 0.587*g/255 + 0.114*b/255) [alpha channel not used - mix with background color before]
        double gray=0.0039216*(0.299*((rgb>>16) & 0xFF) + 0.587*((rgb>>8) & 0xFF) + 0.114*(rgb & 0xFF));
        unsigned int g=(unsigned int) ((maxval+1)*gray);if(g>maxval) g=maxval;
        return g;
    };
    
    //color functions
    double inline gnufunc(double val, int ID); //gnu(plot) color functions for single components
    
    
    //gradient functions
    dcolor getgnugradientcolor(double val,int rf,int gf=-1,int bf=-1,int af=0);
    
    unsigned int loadgradient(string fn); //returns the ID, does not select
    unsigned int definegradient(int N,Xcolor *cols); //sets a defined gradient (array cols is used and not copied -> deallocation by user), returns ID, no select
    int selectgradient(unsigned int ID); //built-in or user defined, return=0 success
    int selectgradientidx(unsigned int idx); //built-in or user defined, return=0 success
    
    int getgradientID() {return gradID;};
    int getgradientidx() {return gradidx;};
    int getgradientlen() {return grads[gradidx].N;};
    int getgradienttype() {return grads[gradidx].type;};
    int getgradientN() {return Ngrads;};
    
    dcolor getgradientcolor(double val,bool periodic=false); //uses the currently selected gradient
    
    
private:
    int Ngrads, gradidx,gradID;
    gradient *grads;
    
    
    void cvalidate(dcolor &col)
    {if(col.R<0.0) col.R=0.0;else if(col.R>1.0) col.R=1.0;
        if(col.G<0.0) col.G=0.0;else if(col.G>1.0) col.G=1.0;
        if(col.B<0.0) col.B=0.0;else if(col.B>1.0) col.B=1.0;
        if(col.A<0.0) col.A=0.0;else if(col.A>1.0) col.A=1.0;
    };
    

};

//simple netpbm graphics files with PXM extension
//supports binary & ascii: PPM (P3, P6), PGM (P2, P5), PBM (P1, P4), and PXM (P0)
//PAM (P7) is not supported
class PXMfile {
public:
    PXMfile(string fn, int type=PXM_P6);
    ~PXMfile();
    
    void setname(string fn) {name=fn;};
    void setsize(unsigned int dx,unsigned int dy,unsigned int dz=1) {Lx=dx;Ly=dy;Lz=dz;};
    void setascii(bool ascii=true) {if(ptype==PXM_P0) binary=!ascii;};
    void settype(int type) {ptype=type;binary=true;if((ptype>=PXM_P1) && (ptype<=PXM_P3)) binary=false;};
    void setmaxgray(unsigned int mg,double BWth=0.5) {maxgray=mg;BWthres=BWth;};
    void setbits(unsigned int b) {bits=b;};
    
    int writefile(dcolor *col,int N);
    int writefile(double *gray,int N);
    int writefile(unsigned int *rgba,int N,bool isrgb=true); //if isrgb=false, assume gray (no rescaling, just cut to "bits")
    int writefile(unsigned char *gray,int N); //no rescaling, needs to be consistent with maxgray (bits are reduced to 8 if larger)
    
    //for sequential (animated/movie) output (P0 only), Lz is number of frames - needs to be known upfront
    void setfrate(double fr) {frate=fr;};
    int startanifile();
    int writeframe(dcolor *col,int N);
    int endanifile() {if(fop) {FileClose(f);fop=false; if(fcount<Lz) return -2;} else return -1;return 0;};
    
private:
    int ptype,cmodel;
    unsigned int Lx,Ly,Lz,maxgray,channels,bits,fcount,ctype;
    double BWthres;
    bool binary,fop;
    double frate;
    fHandle f;
    string name;
    colorfunction *colf;
    
    int checkendian();
    
    string composehead();
    string getext() {return PXMext[ptype];};
    
    int convertbitarray(unsigned char *input,int N,unsigned char *output);
    
    void writeascii(fHandle af,unsigned char* buf,int N);
  
};

/* PXM files
 - supports alpha channel, or CMYK etc. (number of components/channels is variable)
 - supports 16 & 32 colors components (i.e. up to 128bit color per pixel: RGBA) [16&32 bit needs endian info for binary format]
 - full 3D images: height, width and depths (can be used to store movies if interpreted as sequence of frames, use keyword framerate)
 - can encode also P1 through P6 formats
 - 
 
Format:

--------------
P0
SIZE <width> <height> [<depth>]
CH <number of channels>
BITS <1,8,16, or 32>
[MAX <max gray value>] //only for gray images
DATA <A or B>
[ENDIAN <B or L>] //for binary, bits=16 or 32
[MODEL <RGB, GRAY, BW, RGBA, CMYL, HSV, ...>]
[FRAMERATE <rate, float number as frames per second>]
[USER_DEFINED_PARAMETERS]
ENDHDR
<data>
---------------
The whitespace charaters after P0 determine the end of line encoding: UNIX, DOS, MAC (LF,CR+LF, CR).
After "ENDHDR" one EOL marker is expected before the (binary) data starts
"DATA A" is for ASCII, "DATA B" for binary
There should be not data after the <data> part
The header info order (lines between "P0" and "ENDHDR") is arbitrary, but SIZE, CH, BITS, and DATA (and the keyword "ENDHDR") are required

*/
 


#endif // graphicstools_H
