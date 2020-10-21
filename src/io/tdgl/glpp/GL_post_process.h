//
//  GL_post_process.h
//  
//
//  Created by Glatz, Andreas on 2013-11-06.
//
//

#ifndef ____GL_post_process__
#define ____GL_post_process__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;


//constants&macros, mostly from GLtools.h

#define EPS 1E-10

#define PI	3.1415926535897932384626433832795
#define TWOPI 6.28318530717958647692528676655901
#define INVTWOPI 0.15915494309189533576888376337251
#define SQRT2 1.4142135623730950488016887242097
#define PI180 0.017453292519943295769236907684886 //PI/180

#define ABS(a) ((a)<0.0?(-(a)):(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

//---------------------------------------------------------------------------
#define Bcomp_X 1 //0b001
#define Bcomp_Y 2 //0b010
#define Bcomp_Z 4 //0b100
#define Bcomp_XY 3 //0b011
#define Bcomp_XZ 5 //0b101
#define Bcomp_YZ 6 //0b110
//---------------------------------------------------------------------------


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

//make a b-periodic
#define PBC(AA,BB) {while(AA>=(BB)) AA-=(BB);while(AA<0) AA+=(BB);}
#define k_INDEX(i,L) ((i)<=((L)/2)?(i):((i)-(L)))
//---------------------------------------------------------------------------



//OP - order parameter (complex scalar field, or (rho,phi) tupel)
//VP - vector potential (vector field)
//SP - scalar potential (scalar field)
//SC - supercurrent (vector field)
//NC - normal current (vector field)
//SNC - total current = SC+NC (vector field)
//TC  - epsilon (Tc modulation) (scalar field)

#define READ_single 4
#define READ_double 8

#define READ_OP_reim  0x0100
#define READ_OP_ap    0x0200
#define READ_OP       0x0300 //mask to check the OP type
#define READ_VP       0x0400
#define READ_SC       0x0800
#define READ_SP       0x1000
#define READ_TC       0x8000

#define WRITE_OP_reim 0x0100
#define WRITE_OP_amp2 0x0200 //amplitude^2 only
#define WRITE_VP      0x0400
#define WRITE_SC      0x0800
#define WRITE_SP      0x1000
#define WRITE_OP_ph   0x2000 //phase only
#define WRITE_TC      0x8000


typedef struct {
	double re;
	double im;
} COMPLEX;

typedef struct {
	float re;
	float im;
} fCOMPLEX;


typedef struct {
    double time,zaniso,Bx,By,Bz;
	double vfrac,psisqfrac,gradfrac;
    double Js_av,Js_min,Js_max,Js_stddev; //these are actually for |J_s|^2
} Adata;




//reads single and double precision data files
//writes single and double precision data files
//internally double precision is used

class GLPP {
public:
    
    GLPP();
    ~GLPP();
    
    void resetdata();
    
    
    
    void setSV(int ox,int oy,int oz,int nx,int ny,int nz) {SVOx=ox;SVOy=oy;SVOz=oz;SVNx=nx;SVNy=ny;SVNz=nz;};
    void loadparams(int n, char* argv[]);
    int run();
    
    
    //load functions
    int loadarray(string fn);
    int loadBDAT(string fn);
    
    int loadASCIImatrix(string fn,double* &data,int &cols,int &rows,int skiprows=0);

    
    //data processing
    int calc_vector_pot_kappa_inf(); //allocate and calculate the vector potential
    
    int calc_current(double *gradsq=NULL,bool calcnormal=false); //allocate and calculate the (super)currents using the vector potential if allocated or magnetic field in kappa=inf limit
    
    
    //data analysis functions
    int analysis(Adata &adat,double *gradsq=NULL,double psi2thres=0.1); //uses subvolume information if defined

    
    
    //helper
    void print_info();
    int size() {return NN;};
    
    
    //output in BDAT format
    int save_compressed_data(int strideX,int strideY,int strideZ); //saves a "compressed" data set (e.g., remove every 2nd grid point by strideX=2), boundary conditions might be undefined after compression, used subvolume information if defined, uses readdatatype information
    int save_user(unsigned int user_readdatatype); //uses a user defined readdatatype instead of original (can generate A & Js, but not OP if not read)
   
    
    void savearray(string fn,bool binary=false,bool singleprec=false);
    
    //XMDF output
    int write_XMDF(int frame, int totalframes,unsigned int writedata);

    //void output(COMPLEX *Gz,string fn,int omode=0x010101,simdata *sd=NULL);
    //void writeBM(string fn,unsigned int *rgb,int nx,int ny,int deco);
    //void writeBM_real(string fn,REAL *a,REAL m,REAL M,int nx,int ny);
    
    
    void delaunay_analysis();

    
// private:
    int numinputfiles;  //number of input files
    string *inputfiles; //name of input files
    string outputfileprefix; //used only for common output from all input files, e.g. XDMF files
    
    int action; //type of post processing
    
    //system parameters
	int dim,Nx,Ny,Nz,NN,btype;
    double Lx,Ly,Lz,zaniso,time;
    double dx,dy,dz;
    unsigned int readdatatype; //lowest byte 4 or 8 for single or double prec, bit 8&9: 0 no order parameter, 1 re&im order parameter, 2 amph/phase OP, bit 10: 1 if vector pot read; bit 11: 1 if supercurrent read
    
    //for subvolume calculations, used when saving data as well
    int SVOx,SVOy,SVOz; //subvolume origin
    int SVNx,SVNy,SVNz; //size of subvolume
    
    //control parameter
    //megnetic field onlt valid for kappa->inf, otherwise we need to read the vector potential
    double t,Bx,By,Bz,Jex,fT,kappa;
    
    //calculated values
    double KExdot,KEx;
    
    //data structures, working format: double precision
    
    //order paramter
    //always in re & im format (might need conversion while reading/writing)
    COMPLEX *psi;
    //vector potential & (super)current, sclar pot, Tc modulation
    double *Ax,*Ay,*Az,*Jx,*Jy,*Jz,*mu,*epsilon;
    
};


#endif /* defined(____GL_post_process__) */
