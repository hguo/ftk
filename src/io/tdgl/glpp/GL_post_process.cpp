#define GLPPversion 0x00010200

//
//  GL_post_process.cpp
//  
//
//  Created by Glatz, Andreas on 2013-11-06.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "paramfile.h"
#include "fileutils.h"
#include "datafile.h"
#include "stringutils.h"
#include "graphicstools.h"


#include "GL_post_process.h"
//---------------------------------------------------------------------------

/*

actions (* not implemented, + limited):

//all 1xx & 2xx use subvolume information if present (origins might be moved)
//input binary format is preserved, but converted into double internally
100*: read each file (format detected) and write to BDAT format
101+: read each file, calculate supercurrent and save to BDAT (keeps all existing data)
102*: read each file, calculate vector potential (only infinite kappa limit) and save to BDAT
103*: combined 101&102

200: convert all input files into one XMDF, all binary arrays as input, filenames get new structure
201: convert all input files into one XMDF, all binary arrays as input plus supercurrent
 
210: read each file and save to single XMDF, amplitude^2 only
211: read each file, calculate supercurrent and save to single XMDF, supercurrent only saved

//1xxx function use subvolume info
1002: analysis of supercurrent and total energy
 
*/


//---------------------------------------------------------------------------
// data input/output
//---------------------------------------------------------------------------

GLPP::GLPP() {
    dim=2;
    Nx=16;
    Ny=16;
    Nz=1;
    NN=Nx*Ny*Nz;
    btype=0x000000;
    Lx=8.0;
    Ly=8.0;
    Lz=0.0;
    zaniso=1.0;
    kappa=1e6;
    
    dx=Lx/Nx;
    dy=Ly/Ny;
    dz=Lz/Nz;
    
    t=0.0;
    Bx=0.0;
    By=0.0;
    Bz=0.0;
    Jex=0.0;
    fT=0.0;
    
    //subvolume
    SVOx=SVOy=SVOz=0;
    SVNx=SVNy=SVNz=0;
    
    
    readdatatype=0; //nothing read
    psi=NULL;
    Ax=Ay=Az=NULL;
    Jx=Jy=Jz=NULL;
    mu=NULL;
    epsilon=NULL;
    
    inputfiles=NULL;
};


void GLPP::resetdata() {
    if(psi!=NULL) delete[] psi;
    if(Ax!=NULL) delete[] Ax;
    if(Ay!=NULL) delete[] Ay;
    if(Az!=NULL) delete[] Az;
    if(Jx!=NULL) delete[] Jx;
    if(Jy!=NULL) delete[] Jy;
    if(Jz!=NULL) delete[] Jz;
    if(mu!=NULL) delete[] mu;
    if(epsilon!=NULL) delete[] epsilon;
    
    
    psi=NULL;
    Ax=Ay=Az=NULL;
    Jx=Jy=Jz=NULL;
    mu=NULL;
    epsilon=NULL;
    readdatatype=0;
}


GLPP::~GLPP() {
    resetdata();
    if(inputfiles!=NULL) delete[] inputfiles;
};

//---------------------------------------------------------------------------
// data processing
//---------------------------------------------------------------------------


int GLPP::calc_vector_pot_kappa_inf() {
    int i,j,k,m,M;
    double x,y;
    if(Ax!=NULL) return -1; //if a vector potential is already present it is not recalculated
    if((readdatatype&0xFF)==0) return -2; //without valid input data parameter, we cannot calculate a vector potential
    
    Ax=new double[NN];
    Ay=new double[NN];
    if(dim==3) Az=new double[NN];
    
    //we need the type of gauge used for the link variables
    if((dim==3) && ABS(By)>EPS) {
        M=Ny*Nz;
        for(i=0;i<Nx;i++) {
            x=dx*(i-0.5*Nx);
            for(j=0;j<M;j++) {
                m=i+j*Nx;
                Ax[m]=0.0;
                Ay[m]=Bz*x;
                Az[m]=-By*x;}
        }
    }
    else if(dim==3) {
        for(j=0;j<Ny;j++) {
            y=dy*(j-0.5*Ny);
            for(i=0;i<Nx;i++) for(k=0;k<Nz;k++)
            {m=i+Nx*(j+k*Ny);
                Ax[m]=-Bz*y;
                Ay[m]=0;
                Az[m]=Bx*y;}
        }
        
    }
    else if(dim==2) {
        for(j=0;j<Ny;j++) {
            y=dy*(j-0.5*Ny);
            for(i=0;i<Nx;i++)
            {m=i+Nx*j;
                Ax[m]=-Bz*y;
                Ay[m]=0;
            }
        }
    }
    return 0;
};

//---------------------------------------------------------------------------

int GLPP::calc_current(double *gradsq,bool calcnormal) { // uses the vector potential if present, Note: the term $\partial_t {\tilde A}$ is not calculated for the normal part
    int i,j,k,idx,m,p,lp,lm,bc;
    double x,y,v,dx2i,dy2i,dz2i;
    COMPLEX z,zm,zp,U,UK,zd,QP;
    unsigned int gauge,Bcomp;
    
    if(Jx!=NULL) return -1; //if a supercurrent is already present it is not recalculated
    if(psi==NULL) return -2; //without order parameter, we cannot calculate the supercurrent
    
    
    Bcomp=0;
    gauge=0;
    if(ABS(Bx)>EPS) Bcomp+=Bcomp_X;
    if(ABS(By)>EPS) Bcomp+=Bcomp_Y;
    if(ABS(Bz)>EPS) Bcomp+=Bcomp_Z;
    
    if(Ax==NULL) {
        if(Bcomp&Bcomp_Y) gauge=1; //x-dependent gauge with By,Bz
        else if(Bcomp>0) gauge=2; //y-dependent, even when only Bz!=0
    }
    else gauge=16; //use the vector potential
    
    if((gauge<16) && (kappa<5e5)) return -3; //for finite kappa, we need the vector potential, other wise Js cannot be calculated
    
    Jx=new double[NN];
    Jy=new double[NN];
    if(dim==3) Jz=new double[NN];
    
    
    dx2i=1/(2*dx);
    dy2i=1/(2*dy);
    dz2i=1/(2*dz*zaniso);
    
    if(calcnormal && (mu==NULL)) calcnormal=false; //we cannot calculate the normal part w/o the vector potential
    
    //"K-LV"
    x=dx*KEx;
    UK.re=cos(x);UK.im=sin(x);
    
    if(dim==3) {
        for(k=0;k<Nz;k++) {
            for(j=0;j<Ny;j++) {
                for(i=0;i<Nx;i++) {
                    idx=i+Nx*(j+Ny*k);
                    z=psi[idx];
                    
                    
                    //---- x-direction
                    bc=(btype&0xFF);  //x-boundary condition (0 - no current \partial_x\psi=0 & \partial_x\mu=0, 1 - periodic)
                    
                    if((bc==0) && ((i==0) || (i==(Nx-1)))) Jx[idx]=0;  //definition of the no current condition
                    else {
                        lp=i+1;if(lp==Nx) lp=0; //the no-current case is already handled
                        lm=i-1;if(lm<0) lm=Nx-1;
                        m=lm+Nx*(j+k*Ny);
                        p=lp+Nx*(j+k*Ny);
                        zp=psi[p];zm=psi[m];
                        
                        if((gauge==0) || (gauge==1)) { //no LV, except the "UK" factor needed
                            //nothing to do
                        } else if(gauge==2) { //here we use QP-BC
                            //modify the zp & zm terms
                            y=(j-0.5*Ny)*dy*Bz*dx;
                            U.re=cos(y);
                            U.im=sin(y);
                            x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                            x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                            
                            if(bc==1) {//quasi boundary conditions, we assume that the fields are multiples of the flux quantum
                                if(i==0) { //modify zm by QP
                                    x=(k*dz*By-j*dy*By)*Lx;
                                    QP.re=cos(x);
                                    QP.im=sin(x);
                                    x=zm.re;zm.re=zm.re*QP.re-zm.im*QP.im;zm.im=zm.im*QP.re+x*QP.im;
                                }
                                else if(i==(Nx-1)) { //modify zp by QP*
                                    x=(k*dz*By-j*dy*By)*Lx;
                                    QP.re=cos(x);
                                    QP.im=sin(x);
                                    x=zp.re;zp.re=zp.re*QP.re+zp.im*QP.im;zp.im=zp.im*QP.re-x*QP.im;
                                }
                            }
                        } else {//general case, here we use a periodic A, no QP-BC
                            y=-0.5*dx*(Ax[p]+Ax[m]);
                            U.re=cos(y);
                            U.im=sin(y);
                            x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                            x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                        }
                        //apply UK
                        x=UK.re*(zp.re-zm.re)-UK.im*(zp.im+zm.im);
                        y=UK.re*(zp.im-zm.im)+UK.im*(zp.re+zm.re);
                        //calc Im(psi^(x+iy))
                        v=z.re*y-z.im*x;
                        if(calcnormal) v+=(mu[m]-mu[p]);
                        if(gradsq!=NULL) gradsq[idx]=x*x+y*y; //energy value of the gradient term
                        Jx[idx]=dx2i*v;
                    }
                    
                    
                    
                    //---- y-direction
                    bc=((btype>>8)&0xFF); //y-boundary condition (0 - no current \partial_y\psi=0 & \partial_y\mu=0, 1 - quasi-periodic)
                    
                    if((bc==0) && ((j==0) || (j==(Ny-1)))) Jy[idx]=0;  //definition of the no current condition
                    else {
                        lp=j+1;if(lp==Ny) lp=0;
                        lm=j-1;if(lm<0) lm=Ny-1;
                        m=i+Nx*(lm+k*Ny);
                        p=i+Nx*(lp+k*Ny);
                        zp=psi[p];
                        zm=psi[m];
                        
                        if((gauge==0) || (gauge==1)) {//LV
                            x=-(i-0.5*Nx)*dx*Bz*dy;
                            U.re=cos(x);
                            U.im=sin(x);
                            x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                            x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                        } else if(gauge==2) {//QPBC
                            if((bc==1) && (j==0)) { //affects zm by QP
                                y=(i*dx*Bz-k*dz*Bx)*Ly;
                                QP.re=cos(y);
                                QP.im=sin(y);
                                x=zm.re;zm.re=zm.re*QP.re-zm.im*QP.im;zm.im=zm.im*QP.re+x*QP.im;
                            }
                            else if((bc==1) && (j==(Ny-1))) { //affects zp by QP*
                                y=(i*dx*Bz-k*dz*Bx)*Ly;
                                QP.re=cos(y);
                                QP.im=sin(y);
                                x=zp.re;zp.re=zp.re*QP.re+zp.im*QP.im;zp.im=zp.im*QP.re-x*QP.im;
                            }
                        } else {//general case, here we use a periodic A, no QP-BC
                            y=-0.5*dy*(Ay[p]+Ay[m]);
                            U.re=cos(y);
                            U.im=sin(y);
                            x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                            x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                        }
                        x=zp.re-zm.re;
                        y=zp.im-zm.im;
                        v=z.re*y-z.im*x;
                        if(calcnormal) v+=(mu[m]-mu[p]);
                        if(gradsq!=NULL) gradsq[idx]+=(x*x+y*y);
                        Jy[idx]=dy2i*v;
                        
                    }
                    
                    //---- z-direction
                    if(dim>2) {
                        bc=((btype>>16)&0xFF); //z-boundary condition (0 - no current, 1 - periodic)
                        
                        if((bc==0) && ((k==0) || (k==(Nz-1)))) Jz[idx]=0;  //definition of the no current condition
                        else {
                            lp=k+1;if(lp==Nz)  lp=0;
                            lm=k-1;if(lm<0) lm=Nz-1;
                            m=i+Nx*(j+lm*Ny);
                            p=i+Nx*(j+lp*Ny);
                            zp=psi[p];
                            zm=psi[m];
                            
                            if((gauge==0) || (gauge==1)) {//LV
                                x=(i-0.5*Nx)*dx*By*dz;
                                U.re=cos(x);
                                U.im=sin(x);
                            } else if(gauge==2) {//LV
                                y=-(j-0.5*Ny)*dy*Bx*dz;
                                U.re=cos(y);
                                U.im=sin(y);
                            } else {//general case, here we use a periodic A, no QP-BC
                                y=-0.5*dz*(Az[p]+Az[m]);
                                U.re=cos(y);
                                U.im=sin(y);
                            }
                            x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                            x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                            x=zp.re-zm.re;
                            y=zp.im-zm.im;
                            v=z.re*y-z.im*x;
                            if(calcnormal) v+=(mu[m]-mu[p]);
                            if(gradsq!=NULL) gradsq[idx]+=(x*x+y*y);
                            Jz[idx]=dz2i*v;
                        }
                    }
                }
            }
        }
    }
    else if(dim==2) {
        
        for(j=0;j<Ny;j++) {
            for(i=0;i<Nx;i++) {
                idx=i+Nx*j;
                z=psi[idx];
                
                //---- x-direction
                bc=(btype&0xFF);  //x-boundary condition (0 - no current \partial_x\psi=0 & \partial_x\mu=0, 1 - periodic)
                
                if((bc==0) && ((i==0) || (i==(Nx-1)))) Jx[idx]=0.0;  //definition of the no current condition
                else {
                    lp=i+1;if(lp==Nx) lp=0; //the no-current case is already handled
                    lm=i-1;if(lm<0) lm=Nx-1;
                    m=lm+Nx*(j+k*Ny);
                    p=lp+Nx*(j+k*Ny);
                    zp=psi[p];zm=psi[m];
                    
                    if((gauge==0) || (gauge==1)) { //no LV, except the "UK" factor needed
                        //nothing to do
                    } else if(gauge==2) { //here we use QP-BC
                        //modify the zp & zm terms
                        y=(j-0.5*Ny)*dy*Bz*dx;
                        U.re=cos(y);
                        U.im=sin(y);
                        x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                        x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                        //boundary conditions
                        if(bc==1) {//quasi boundary conditions, we assume that the fields are multiples of the flux quantum
                            if(i==0) { //modify zm by QP
                                x=(k*dz*By-j*dy*By)*Lx;
                                QP.re=cos(x);
                                QP.im=sin(x);
                                x=zm.re;zm.re=zm.re*QP.re-zm.im*QP.im;zm.im=zm.im*QP.re+x*QP.im;
                            }
                            else if(i==(Nx-1)) { //modify zp by QP*
                                x=(k*dz*By-j*dy*By)*Lx;
                                QP.re=cos(x);
                                QP.im=sin(x);
                                x=zp.re;zp.re=zp.re*QP.re+zp.im*QP.im;zp.im=zp.im*QP.re-x*QP.im;
                            }
                        }
                    } else {//general case, here we use a periodic A, no QP-BC, not valid for 2D really
                        y=-0.5*dx*(Ax[p]+Ax[m]);
                        U.re=cos(y);
                        U.im=sin(y);
                        x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                        x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                    }
                    //apply UK
                    x=UK.re*(zp.re-zm.re)-UK.im*(zp.im+zm.im);
                    y=UK.re*(zp.im-zm.im)+UK.im*(zp.re+zm.re);
                    //calc Im(psi^(x+iy))
                    v=z.re*y-z.im*x;
                    if(calcnormal) v+=(mu[m]-mu[p]);
                    if(gradsq!=NULL) gradsq[idx]=(x*x+y*y);
                    Jx[idx]=dx2i*v;
                }
                
                
                
                //---- y-direction
                bc=((btype>>8)&0xFF); //y-boundary condition (0 - no current \partial_y\psi=0 & \partial_y\mu=0, 1 - quasi-periodic)
                
                if((bc==0) && ((j==0) || (j==(Ny-1)))) Jy[idx]=0;  //definition of the no current condition
                else {
                    lp=j+1;if(lp==Ny) lp=0;
                    lm=j-1;if(lm<0) lm=Ny-1;
                    m=i+Nx*(lm+k*Ny);
                    p=i+Nx*(lp+k*Ny);
                    zp=psi[p];
                    zm=psi[m];
                    
                    if((gauge==0) || (gauge==1)) {//LV
                        x=-(i-0.5*Nx)*dx*Bz*dy;
                        U.re=cos(x);
                        U.im=sin(x);
                        x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                        x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                    } else if(gauge==2) {//QPBC
                        if((bc==1) && (j==0)) { //affects zm by QP
                            y=(i*dx*Bz-k*dz*Bx)*Ly;
                            QP.re=cos(y);
                            QP.im=sin(y);
                            x=zm.re;zm.re=zm.re*QP.re-zm.im*QP.im;zm.im=zm.im*QP.re+x*QP.im;
                        }
                        else if((bc==1) && (j==(Ny-1))) { //affects zp by QP*
                            y=(i*dx*Bz-k*dz*Bx)*Ly;
                            QP.re=cos(y);
                            QP.im=sin(y);
                            x=zp.re;zp.re=zp.re*QP.re+zp.im*QP.im;zp.im=zp.im*QP.re-x*QP.im;
                        }
                    } else {//general case, here we use a periodic A, no QP-BC
                        y=-0.5*dy*(Ay[p]+Ay[m]);
                        U.re=cos(y);
                        U.im=sin(y);
                        x=zp.re;zp.re=U.re*x-U.im*zp.im;zp.im=U.re*zp.im+U.im*x;
                        x=zm.re;zm.re=U.re*x+U.im*zm.im;zm.im=U.re*zm.im-U.im*x;
                    }
                    x=zp.re-zm.re;
                    y=zp.im-zm.im;
                    v=z.re*y-z.im*x;
                    if(calcnormal) v+=(mu[m]-mu[p]);
                    if(gradsq!=NULL) gradsq[idx]+=(x*x+y*y);
                    Jy[idx]=dy2i*v;
                }
            }
        }
    }
    return 0;
}
//---------------------------------------------------------------------------

int GLPP::analysis(Adata &adat,double *gradsq,double psi2thres) {
    int ox,oy,oz,ex,ey,ez;
    int i,j,k,idx,vcount;
    double vfrac,psi2av,grad2av,Js2av,Js2min,Js2max,Js4av,Jabs2,x;
    COMPLEX z;
    
    ox=oy=oz=0;
    ex=Nz;ey=Ny;ez=Nz;
    if(SVNx>0) {
        ox=SVOx;oy=SVOy;oz=SVOz;
        ex=ox+SVNx;ey=oy+SVNy;ez=oz+SVNz;
        if(ex>Nx) ex=Nx;
        if(ey>Ny) ey=Ny;
        if(ez>Nz) ez=Nz;
    }
    
    if(dim==2) {oz=0;ez=1;}
    
    vcount=0;
    psi2av=0.0;
    grad2av=0.0;
    Js2av=0.0;
    Js2min=1e6;
    Js2max=0.0;
    Js4av=0.0;
    
    
    for(k=oz;k<ez;k++) {
        for(j=oy;j<ey;j++) {
            for(i=ox;i<ex;i++) {
                idx=i+Nx*(j+k*Nz);
                if(Jx!=NULL) {
                    x=Jx[idx];Jabs2=x*x;
                    x=Jy[idx];Jabs2+=x*x;
                    if(dim==3) {x=Jz[idx];Jabs2+=x*x;}
                    if(Jabs2<Js2min) Js2min=Jabs2;
                    if(Jabs2>Js2max) Js2max=Jabs2;
                    Js2av+=Jabs2;
                    Js4av+=Jabs2*Jabs2;
                }
                if (gradsq!=NULL) {x=gradsq[idx];grad2av+=x;}
                if(psi!=NULL) {
                    z=psi[idx];
                    x=z.re*z.re+z.im*z.im;
                    if(x<psi2thres) vcount++;
                    psi2av+=x;
                }
            }
        }
    }
    x=(ex-ox)*(ey-oy)*(ez-oz);
    x=1/x;
    vfrac=x*vcount;
    psi2av*=x;
    Js2av*=x;
    Js4av*=x;
    grad2av*=x;
    
    adat.time=time;
    adat.zaniso=zaniso;
    adat.Bx=Bx;adat.By=By;adat.Bz=Bz;
    adat.vfrac=vfrac;
    adat.psisqfrac=psi2av;
    adat.gradfrac=grad2av;
    adat.Js_av=Js2av;
    adat.Js_min=Js2min;
    adat.Js_max=Js2max;
    adat.Js_stddev=sqrt(Js4av-Js2av*Js2av);
    
    //printf("%le\t%le\t%le\n",time,energy,vfrac);
    return 0;
}


//---------------------------------------------------------------------------
// data input/output
//---------------------------------------------------------------------------


//old data format
void GLPP::savearray(string fn,bool binary,bool singleprec)
{
    fHandle f;
    string s;
    int i,M,FPs,OPtype;
    float ftmp;
    fCOMPLEX *fctmp;
    
    s="CA01\n";
    if(binary) s="CA02";
    f=FileCreate(fn);
    FileWrite(f,s.c_str(),s.length());
    
    if(singleprec) FPs=sizeof(float);
    else FPs=sizeof(double);
    
    OPtype=1;
    
    M=Nx;
    if(dim>1) M*=Ny;
    if(dim>2) M*=Nz;
    if(binary)
    {
        i=0xFEFF; //BOM
        FileWrite(f,&i,sizeof(int));
        FileWrite(f,&dim,sizeof(int));
        
        FileWrite(f,&FPs,sizeof(int));
        
        FileWrite(f,&Nx,sizeof(int));
        if(singleprec) {
            ftmp=Lx; //conversion
            FileWrite(f,&ftmp,FPs);
        }
        else FileWrite(f,&Lx,FPs);
        
        if(dim>1) {
            FileWrite(f,&Ny,sizeof(int));
            if(singleprec) {ftmp=Ly;FileWrite(f,&ftmp,FPs);}
            else FileWrite(f,&Ly,FPs);}
        if(dim>2) {
            FileWrite(f,&Nz,sizeof(int));
            if(singleprec) {ftmp=Lz;FileWrite(f,&ftmp,FPs);}
            else FileWrite(f,&Lz,FPs);}
        
        //more header information
        i=8*FPs+2*sizeof(int); //length
        FileWrite(f,&i,sizeof(int));
        
        if(singleprec) {
            ftmp=t;FileWrite(f,&ftmp,FPs); //time
            ftmp=fT;FileWrite(f,&ftmp,FPs); //fluctuation amplitude
            ftmp=Bx;FileWrite(f,&ftmp,FPs); //magnetic field in x-direction
            ftmp=By;FileWrite(f,&ftmp,FPs); //magnetic field in y-direction
            ftmp=Bz;FileWrite(f,&ftmp,FPs); //magnetic field in z-direction
            ftmp=Jex;FileWrite(f,&ftmp,FPs); //current in x-direction
            FileWrite(f,&btype,sizeof(int)); //boundary conditions
            FileWrite(f,&OPtype,sizeof(int)); //Re&Im or Ampl&Phase
            ftmp=KEx;FileWrite(f,&ftmp,FPs); //integral of the voltage
            ftmp=KExdot;FileWrite(f,&ftmp,FPs); //voltage
        }
        else  {
            FileWrite(f,&t,FPs); //time
            FileWrite(f,&fT,FPs); //fluctuation amplitude
            FileWrite(f,&Bx,FPs); //magnetic field in x-direction
            FileWrite(f,&By,FPs); //magnetic field in y-direction
            FileWrite(f,&Bz,FPs); //magnetic field in z-direction
            FileWrite(f,&Jex,FPs); //current in x-direction
            FileWrite(f,&btype,sizeof(int)); //boundary conditions
            FileWrite(f,&OPtype,sizeof(int)); //Re&Im or Ampl&Phase
            FileWrite(f,&KEx,FPs); //integral of the voltage
            FileWrite(f,&KExdot,FPs); //voltage
        }

        if(singleprec) {
            fctmp=new fCOMPLEX[M];
            for(i=0;i<M;i++) {fctmp[i].re=psi[i].re;fctmp[i].im=psi[i].im;}
            FileWrite(f,fctmp,M*sizeof(fCOMPLEX));
            delete[] fctmp;
        }
        else
            FileWrite(f,psi,M*sizeof(COMPLEX));
    }
    else //does not correspond to the binary format
    {
        s=IntToStr(dim)+" "+IntToStr(Nx);
        if(dim>1) s=s+" "+IntToStr(Ny);
        if(dim>2) s=s+" "+IntToStr(Nz);
        s=s+"\n";
        FileWrite(f,s.c_str(),s.length());
        for(i=0;i<M;i++)
        {s=FloatToStr(psi[i].re)+" "+FloatToStr(psi[i].im)+"\n";
            FileWrite(f,s.c_str(),s.length());
        }
    }
    FileClose(f);
}
//---------------------------------------------------------------------------

//not compatible with the above
int GLPP::loadarray(string fn) {
    fHandle f;
    string s;
    char *buf;
    int i,M,FPs,len,OPtype;
    bool binary,singleprec;
    float ftmp;
    fCOMPLEX *fctmp;
    
    f=FileOpen(fn,fmOpenRead);
    if(f==-1) return -1;
    i=FileSize(f);
    buf=new char[i];
    FileRead(f,buf,i);
    FileClose(f);
    
    s="xxxx";for(i=0;i<4;i++) s[i]=buf[i];
    if(s=="CA01") binary=false;
    else if(s=="CA02") binary=true;
    else {delete[] buf;return -2;}
    
    if(psi!=NULL) delete[] psi;
    psi=NULL;
    
    if(binary)
    {
        memmove(&i,&buf[4],4);
        if(i!=0xFEFF) {delete[] buf;return -3;} //check endian ...
        i=8; //offset
        memmove(&dim,&buf[i],4);i+=4; //dimension
        memmove(&FPs,&buf[i],4);i+=4; //floating point size
        singleprec=(FPs==sizeof(float));
        memmove(&Nx,&buf[i],4);i+=4;
        if(singleprec) {memmove(&ftmp,&buf[i],FPs);Lx=ftmp;}
        else memmove(&Lx,&buf[i],FPs);
        i+=FPs;
        M=Nx;
        if(dim>1) {
            memmove(&Ny,&buf[i],4);i+=4;
            if(singleprec) {memmove(&ftmp,&buf[i],FPs);Ly=ftmp;}
            else memmove(&Ly,&buf[i],FPs);
            i+=FPs;
            M*=Ny;
        }
        if(dim>2) {
            memmove(&Nz,&buf[i],4);i+=4;
            if(singleprec) {memmove(&ftmp,&buf[i],FPs);Lz=ftmp;}
            else memmove(&Lz,&buf[i],FPs);
            i+=FPs;
            M*=Nz;
        }
        memmove(&len,&buf[i],4);i+=4; //length of additional information
        if(len==8*FPs+2*sizeof(int)) {
            //corresponds to format of the current writer:
            //t,fT,Bx,By,Bz,Jx,btype,OPtype,Kex,Kexdot
            if(singleprec) {
                memmove(&ftmp,&buf[i],FPs);i+=FPs;t=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;fT=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;Bx=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;By=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;Bz=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;Jex=ftmp;
                memmove(&btype,&buf[i],4);i+=4;
                memmove(&OPtype,&buf[i],4);i+=4;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;KEx=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;KExdot=ftmp;
            }
            else  {
                memmove(&t,&buf[i],FPs);i+=FPs;
                memmove(&fT,&buf[i],FPs);i+=FPs;
                memmove(&Bx,&buf[i],FPs);i+=FPs;
                memmove(&By,&buf[i],FPs);i+=FPs;
                memmove(&Bz,&buf[i],FPs);i+=FPs;
                memmove(&Jex,&buf[i],FPs);i+=FPs;
                memmove(&btype,&buf[i],4);i+=4;
                memmove(&OPtype,&buf[i],4);i+=4;
                memmove(&KEx,&buf[i],FPs);i+=FPs;
                memmove(&KExdot,&buf[i],FPs);i+=FPs;
            }
        }
        else if(len==4*FPs) { //older format with t,fT,Bz,Jx
            if(singleprec) {
                memmove(&ftmp,&buf[i],FPs);i+=FPs;t=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;fT=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;Bz=ftmp;
                memmove(&ftmp,&buf[i],FPs);i+=FPs;Jex=ftmp;
            }
            else  {
                memmove(&t,&buf[i],FPs);i+=FPs;
                memmove(&fT,&buf[i],FPs);i+=FPs;
                memmove(&Bz,&buf[i],FPs);i+=FPs;
                memmove(&Jex,&buf[i],FPs);i+=FPs;
            }
        }
        psi=new COMPLEX[M];
        if(singleprec) {
            fctmp=(fCOMPLEX *) &buf[i];
            for(i=0;i<M;i++) {psi[i].re=fctmp[i].re;psi[i].im=fctmp[i].im;}
        }
        else
            memmove(psi,&buf[i],M*sizeof(COMPLEX));
    }
    else
    {
        //ASCII reader to be implemented
    } 
    
    delete[] buf;
    return 0;
}


int GLPP::loadASCIImatrix(string fn,double* &data,int &cols,int &rows,int skiprows) {
    fHandle f;
    string s;
    char *buf;
    char c;
    double *tmp;
    int i,len,l,pos,end,n,m;
    bool del;
    
    f=FileOpen(fn,fmOpenRead);
    if(f==-1) return -1;
    len=FileSize(f);
    buf=new char[len+1];
    FileRead(f,buf,len);
    buf[len]=0;
    FileClose(f);
    
    
    
    //prepare buffer  (remove dos/unix line endings, remove comments)
    del=false;
    for(i=0;i<len;i++)
    {c=buf[i];
        if((c==0x0D) || (c==0x0A)) {buf[i]=0;del=false;}
        else if(c=='#') {buf[i]=0;del=true;}
        else if (del) buf[i]=0;
    }
    
    
    //count lines (lines without data are also counted!!
    l=0;
    pos=1;
    while(pos<len)
    {
        while((buf[pos]==0) && (pos<len)) pos++;  //skip 0's
        if(buf[pos]!=0) l++;
        while((buf[pos]!=0) && (pos<len)) pos++; //skip the contents
    }
    
    pos=0;while((buf[pos]==0) && (pos<len)) pos++;
    for(i=0;i<skiprows;i++) {
        while((buf[pos]!=0) && (pos<len)) pos++;
        l--;
        while((buf[pos]==0) && (pos<len)) pos++;
    }
    if(l>0) {
        rows=l;
        n=0;
        //read lines
        while(pos<len)
        {
            //find line
            end=pos;
            while((buf[end]!=0) && (end<len)) end++;
            
            //analyze line
            m=0;
            for(i=pos;i<end;i++)
            {c=buf[i];
                if((c==' ') || (c==0x9) || (c==',') || (c==';') || (c=='"')) buf[i]=0;
                else m++;
            }
            
            if(m>0)
            {
                //trimming
                while(buf[pos]==0) pos++;
                while(buf[end-1]==0) end--;
                
                if(n==0) //initialization of the data structure
                {
                    i=pos;
                    cols=0;
                    //count the word in the first line
                    while(i<end)
                    {
                        while((buf[i]==0) && (i<end)) i++;  //skip 0's
                        if(buf[i]!=0) cols++;
                        while((buf[i]!=0) && (i<end)) i++; //skip the word
                    }
                    data=new double[rows*cols];
                }
                i=pos;
                m=0;
                while(i<end)
                {
                    while((buf[i]==0) && (i<end)) i++;  //skip 0's
                    data[n*cols+m]=StrToFloat(&buf[i]);m++;
                    while((buf[i]!=0) && (i<end)) i++; //skip the word
                }
                n++;
            }
            pos=end+1;while((buf[pos]==0) && (pos<len)) pos++;
        }
        if(n<rows) {
            rows=n;
            tmp=new double[cols*rows];
            memmove(tmp,data,cols*rows*sizeof(double));
            delete[] data;
            data=tmp;
            tmp=NULL;
        }
        
    }
    delete[] buf;
    return 0;
};



//"try" to load the order parameter, works only if systems sizes are compatible
int GLPP::loadBDAT(string fn) {
    double zabs,zph;
    int k;
    datablockinfo db;
    BDATreader *dfr;
    
    readdatatype=0;
    dfr=new BDATreader(fn);
    
    dfr->open();
    
    if(dfr->seek("dim",db)!=0) {delete dfr;return -1;};
    dfr->readblock(&dim,1);

    
    if(dfr->seek("Nx",db)!=0) {delete dfr;return -1;};
    dfr->readblock(&Nx,1);

    
    if(dfr->seek("Ny",db)!=0) {delete dfr;return -1;};
    dfr->readblock(&Ny,1);

    
    if(dim>2) {
        if(dfr->seek("Nz",db)!=0) {delete dfr;return -1;};
        dfr->readblock(&Nz,1);
    } else Nz=1;
    
    NN=Nx*Ny;
    if(dim==3) NN=NN*Nz;
    
    if(dfr->seek("BC",db)==0) dfr->readblock(&btype,1);
    else btype=0;
    
    if(dfr->seek("Lx",db)!=0) {delete dfr;return -1;};
    dfr->readconvert(&Lx,1);
    if(dfr->seek("Ly",db)!=0) {delete dfr;return -1;};
    dfr->readconvert(&Ly,1);
    if(dim>2) {
        if(dfr->seek("Lz",db)!=0) {delete dfr;return -1;};
        dfr->readconvert(&Lz,1);
    } else Lz=1.0;
    
    if(dfr->seek("Bx",db)==0) dfr->readconvert(&Bx,1);
    else Bx=0.0;
    if(dfr->seek("By",db)==0) dfr->readconvert(&By,1);
    else By=0.0;
    if(dfr->seek("Bz",db)==0) dfr->readconvert(&Bz,1);
    else Bz=0.0;

    
    if(dfr->seek("K",db)==0) { //read the voltage integral value
        dfr->readconvert(&KEx,1);
    } else KEx=0.0;
    
    if(dfr->seek("t",db)==0) { //read the voltage integral value
        dfr->readconvert(&time,1);
    } else time=0.0;
    
    if(dfr->seek("V",db)==0) { //read the voltage 
        dfr->readconvert(&KExdot,1);
    } else KExdot=0.0;
    
    
    //finally read the order parameter
    if(dfr->seek("psi",db)!=0) {delete dfr;return -2;};
    
    
    if(db.type==DT_single) readdatatype+=4;
    else readdatatype+=8;
    
    psi=new COMPLEX[NN];
    
    dfr->readconvert((double *)psi);
    
    if(db.ID==2001) { //means OP contains amp^2&phase data -> convert to re&im
        for(k=0;k<NN;k++) {
            zabs=sqrt(psi[k].re);
            zph=psi[k].im;
            psi[k].re=zabs*cos(zph);
            psi[k].im=zabs*sin(zph);
        }
        readdatatype|=READ_OP_ap; //means that the original format was amp2,phase, but internally we always use re,im
    } else readdatatype|=READ_OP_reim;
    
    delete dfr;
    
    
    
    dx=Lx/Nx;if((btype&0x0000FF)==0) dx=Lx/(Nx-1.0);
    dy=Ly/Ny;if((btype&0x00FF00)==0) dy=Ly/(Ny-1.0);
    dz=Lz/Nz;if((dim>2) && ((btype&0xFF0000)==0)) dz=Lz/(Nz-1.0);
    
    
    return 0;
}

//*********** XDMF output - no reader, since not an output format ********
//3D only right now !!!!!
//frame starts at 1! if totalframes>1 frames with frame>1 are appended and file footer is only written for frame=totalframe
//outputfileprefix needs to be valid
int GLPP::write_XMDF(int frame, int totalframes,unsigned int writedata) {
    fHandle f;
    int i;
    string s,topo,geo,origin,fn;
    float *floattmp;
    double *doubletmp;
    fCOMPLEX *fCtmp;
    COMPLEX z;
    bool wfloat=((writedata&0xFF)==4);
    
    s="";
    //Topology and Geometry
    topo=IntToStr(Nz)+" "+IntToStr(Ny)+" "+IntToStr(Nx);
    geo=FloatToStr(dz)+" "+FloatToStr(dy)+" "+FloatToStr(dx);
    origin=FloatToStr(-0.5*Lz)+" "+FloatToStr(-0.5*Ly)+" "+FloatToStr(-0.5*Lx);
    
    if(frame==1) {
        s="<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
        s=s+"<Domain>\n<Topology name=\"system\" TopologyType=\"3DCoRectMesh\" Dimensions=\""+topo+ "\">\n</Topology>\n<Geometry name=\"geo\" Type=\"ORIGIN_DXDYDZ\">\n<!-- Origin -->\n<DataItem Format=\"XML\" Dimensions=\"3\">\n"+origin+"\n</DataItem>\n<!-- DxDyDz -->\n<DataItem Format=\"XML\" Dimensions=\"3\">\n"+geo+"\n</DataItem>\n</Geometry>\n";
        //Grid information, with possibility to add temporal data
        s=s+"<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n\n";
    }
    
    //Frame data
    fn=outputfileprefix+"_"+IntToStrF(frame,4);
    s=s+"<Grid Name=\"T"+IntToStr(frame)+"\" GridType=\"Uniform\">\n<Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n<Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n";
    s=s+"<Time Value=\""+FloatToStr(t)+"\"/>\n\n";
    
    s=s+"<Attribute Center=\"Grid\" Name=\"B\" Type=\"Vector\">\n<DataItem DataType=\"Float\" Dimensions=\"1 3\" Format=\"XML\" Precision=\"8\">\n";
    s=s+FloatToStr(Bx)+" "+FloatToStr(By)+" "+FloatToStr(Bz)+"\n";
    s=s+"</DataItem>\n</Attribute>\n";
    
    s=s+"<Attribute Center=\"Grid\" Name=\"pbc\" Type=\"Vector\">\n<DataItem DataType=\"UChar\" Dimensions=\"1 3\" Format=\"XML\" Precision=\"1\">\n";
    s=s+IntToStr(btype&0xFF)+" "+IntToStr((btype&0xFF00)>>8)+" "+IntToStr((btype&0xFF0000)>>16)+"\n";
    s=s+"</DataItem>\n</Attribute>\n";
    
    s=s+"<Attribute Center=\"Grid\" Name=\"Jxext\" Type=\"Scalar\">\n<DataItem DataType=\"Float\" Dimensions=\"1\" Format=\"XML\" Precision=\"8\">\n";
    s=s+FloatToStr(Jex)+"\n";
    s=s+"</DataItem>\n</Attribute>\n";
    
    s=s+"<Attribute Center=\"Grid\" Name=\"Kx\" Type=\"Scalar\">\n<DataItem DataType=\"Float\" Dimensions=\"1\" Format=\"XML\" Precision=\"8\">\n";
    s=s+FloatToStr(KEx)+"\n";
    s=s+"</DataItem>\n</Attribute>\n";
    
    s=s+"<Attribute Center=\"Grid\" Name=\"V\" Type=\"Scalar\">\n<DataItem DataType=\"Float\" Dimensions=\"1\" Format=\"XML\" Precision=\"8\">\n";
    s=s+FloatToStr(KExdot)+"\n";
    s=s+"</DataItem>\n</Attribute>\n";
    
    s=s+"\n";
    
    if(((writedata&WRITE_OP_amp2)==WRITE_OP_amp2) && (psi!=NULL)) {
        if(wfloat) {
            floattmp=new float[NN];
            for(i=0;i<NN;i++) {z=psi[i];floattmp[i]=(float) (z.re*z.re+z.im*z.im);}
            f=FileCreate(fn+"_amp2.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            delete[] floattmp;
        } else {
            doubletmp=new double[NN];
            for(i=0;i<NN;i++) {z=psi[i];doubletmp[i]=(z.re*z.re+z.im*z.im);}
            f=FileCreate(fn+"_amp2.bin");
            FileWrite(f,doubletmp,sizeof(double)*NN);
            FileClose(f);
            delete[] doubletmp;
        }
        s=s+"<Attribute Name=\"rho\" Center=\"Node\">\n<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+"\">\n"+fn+"_amp2.bin\n</DataItem>\n</Attribute>\n";
    }
    
    if(((writedata&WRITE_OP_reim)==WRITE_OP_reim) && (psi!=NULL)) {
        if(wfloat) {
            floattmp=new float[NN];
            for(i=0;i<NN;i++) {z=psi[i];floattmp[i]=(float) z.re;}
            f=FileCreate(fn+"_re.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            for(i=0;i<NN;i++) {z=psi[i];floattmp[i]=(float) z.im;}
            f=FileCreate(fn+"_im.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            delete[] floattmp;
        } else {
            doubletmp=new double[NN];
            for(i=0;i<NN;i++) {z=psi[i];doubletmp[i]=z.re;}
            f=FileCreate(fn+"_re.bin");
            FileWrite(f,doubletmp,sizeof(double)*NN);
            FileClose(f);
            for(i=0;i<NN;i++) {z=psi[i];doubletmp[i]=z.im;}
            f=FileCreate(fn+"_im.bin");
            FileWrite(f,doubletmp,sizeof(double)*NN);
            FileClose(f);
            delete[] doubletmp;
        }
        s=s+"<Attribute Name=\"re\" Center=\"Node\">\n<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+"\">\n"+fn+"_re.bin\n</DataItem>\n</Attribute>\n";
        s=s+"<Attribute Name=\"im\" Center=\"Node\">\n<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+"\">\n"+fn+"_im.bin\n</DataItem>\n</Attribute>\n";
    }
    if(((writedata&WRITE_OP_ph)==WRITE_OP_ph) && (psi!=NULL)) {
        if(wfloat) {
            floattmp=new float[NN];
            for(i=0;i<NN;i++) {z=psi[i];floattmp[i]=(float) (atan2(z.im,z.re));}
            f=FileCreate(fn+"_ph.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            delete[] floattmp;
        } else {
            doubletmp=new double[NN];
            for(i=0;i<NN;i++) {z=psi[i];doubletmp[i]=atan2(z.im,z.re);}
            f=FileCreate(fn+"_ph.bin");
            FileWrite(f,doubletmp,sizeof(double)*NN);
            FileClose(f);
            delete[] doubletmp;
        }
        s=s+"<Attribute Name=\"phi\" Center=\"Node\">\n<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+"\">\n"+fn+"_ph.bin\n</DataItem>\n</Attribute>\n";
    }
    if(((writedata&WRITE_SC)==WRITE_SC) && (Jx!=NULL)) {
        if(wfloat) {
            floattmp=new float[NN];
            for(i=0;i<NN;i++) floattmp[i]=(float) (Jx[i]);
            f=FileCreate(fn+"_Jx.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            for(i=0;i<NN;i++) floattmp[i]=(float) (Jy[i]);
            f=FileCreate(fn+"_Jy.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            for(i=0;i<NN;i++) floattmp[i]=(float) (Jz[i]);
            f=FileCreate(fn+"_Jz.bin");
            FileWrite(f,floattmp,sizeof(float)*NN);
            FileClose(f);
            delete[] floattmp;
        } else {
            f=FileCreate(fn+"_Jx.bin");
            FileWrite(f,Jx,sizeof(double)*NN);
            FileClose(f);
            f=FileCreate(fn+"_Jy.bin");
            FileWrite(f,Jy,sizeof(double)*NN);
            FileClose(f);
            f=FileCreate(fn+"_Jz.bin");
            FileWrite(f,Jz,sizeof(double)*NN);
            FileClose(f);
        }
        s=s+"<Attribute Name=\"supercurrent\" Center=\"Node\" AttributeType=\"Vector\">\n";
        s=s+"<DataItem Dimensions=\""+topo+" 3\" Function=\"JOIN($0, $1, $2)\" ItemType=\"Function\">\n";
        s=s+"<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+" 1\">\n"+fn+"_Jx.bin\n</DataItem>\n";
        s=s+"<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+" 1\">\n"+fn+"_Jy.bin\n</DataItem>\n";
        s=s+"<DataItem Format=\"Binary\" DataType=\"Float\" Precision=\""+IntToStr(wfloat?4:8)+"\" Endian=\"Little\" Dimensions=\""+topo+" 1\">\n"+fn+"_Jz.bin\n</DataItem>\n";
        s=s+"</DataItem>\n</Attribute>\n";
    }
    s=s+"</Grid>\n\n";
    
    
    //Footer
    if(frame==totalframes)
        s=s+"</Grid>\n</Domain>\n</Xdmf>\n";
    
    
    //output s to xmf file
    if(frame==1) f=FileCreate(outputfileprefix+".xmf");
    else f=FileOpen(outputfileprefix+".xmf",fmOpenAppend);
    if(f==-1) return -1;
    FileWrite(f,s.c_str(),s.length());
    FileClose(f);
    
    return 0;
}



//---------------------------------------------------------------------------

void GLPP::delaunay_analysis() {
    int n,vortices,i,j,k,m;
    fHandle f;
    unsigned int *pixels,*ipts;
    double *matrix;
    int *hist;
    int cols,rows,tc,nh;
    unsigned int p1,p2,p3,pa,pb,pc,pd,px,ic;
    string s;
    double deff;
    bool found;
    
    f=FileCreate(outputfileprefix+".txt");
    for(n=0;n<numinputfiles;n++) {
        vortices=0;
        deff=0.0;
        
        matrix=NULL;
        printf("loading %s: ",inputfiles[n].c_str());
        loadASCIImatrix(inputfiles[n],matrix,cols,rows,1); //expect the imagej delaunary edge syntax: ID x1 y1 x2 y2 in pixels
        printf("%d rows, %d cols\n",rows,cols);
        
        if((cols==5) && (rows>0)) {
            pixels=new unsigned int[rows+rows];
            for(i=0;i<rows;i++) { //transform the x,y corrdinates to single integer
                p1=((unsigned int) (matrix[5*i+1]+0.5));
                p2=((unsigned int) (matrix[5*i+2]+0.5));
                pa=((unsigned int) (matrix[5*i+3]+0.5));
                pb=((unsigned int) (matrix[5*i+4]+0.5));
                //if(i<10) printf("%d %d %d %d\n",p1,p2,pa,pb);
                pixels[i+i]=  p1+p2*65536;
                pixels[i+i+1]=pa+pb*65536;
            }
            delete[] matrix; matrix=NULL;
            ipts=new unsigned int[rows];ic=0; //invalid points
            
            // now remove edges belonging only to one triangle
            i=0;
            m=0;
            while(i<rows) {
                p1=pixels[i+i];
                p2=pixels[i+i+1];
                if((p1!=0) && (p2!=0)) {//valid edge,should always be the case
                    tc=0;
                    for(j=0;j<rows;j++) {
                        if(j==i) j++;
                        pa=pixels[j+j];
                        pb=pixels[j+j+1];
                        found=false;
                        if(pa==p1)      {p3=pb;px=p2;found=true;}
                        else if(pb==p1) {p3=pa;px=p2;found=true;}
                        else if(pa==p2) {p3=pb;px=p1;found=true;}
                        else if(pb==p2) {p3=pa;px=p1;found=true;}
                        if(found) { //one connecting edge found
                            for(k=0;k<rows;k++) {
                                if((k==j) || (k==i)) k++;
                                if((k==j) || (k==i)) k++;
                                pc=pixels[k+k];
                                pd=pixels[k+k+1];
                                if(((pc==p3) && (pd==px)) || ((pc==px) && (pd==p3))) tc++;
                            }
                        }
                    }
                    if(tc!=4) {//invalidate edge and save points for later invalidation in all edges
                        pixels[i+i]=0;pixels[i+i+1]=0;m++;
                        found=false;k=0;while((k<ic)&& !found) {if(ipts[k]==p1) found=true;k++;}
                        if(!found) {ipts[ic]=p1;ic++;}
                        found=false;k=0;while((k<ic)&& !found) {if(ipts[k]==p2) found=true;k++;}
                        if(!found) {ipts[ic]=p2;ic++;}
                    }
                    if(tc>4) printf("not a valid triangulation...\n");
                    //printf("%d ",tc);
                }
                i++;
            }
            printf("%d boundary points and %d boundary edges found\n",ic,m);
            //invalidate all bounadry points in all edges
            for(k=0;k<ic;k++) {
                p1=ipts[k];
                for(i=0;i<rows+rows;i++) if(pixels[i]==p1) pixels[i]=0;
            }
            delete[] ipts;
            
            nh=100;
            hist=new int[nh];
            for(i=0;i<nh;i++) hist[i]=0;
            i=0;
            while(i<rows+rows) {
                p1=pixels[i];
                if(p1!=0) {
                    tc=1; //p1
                    for(j=i+1;j<rows+rows;j++) {
                        p2=pixels[j];
                        if(p1==p2) {tc++;pixels[j]=0;}
                    }
                    //printf("%d ",tc);
                    vortices++;
                    if(tc<nh) hist[tc]++;
                    else hist[0]++;
                }
              i++;
            }
            if(vortices>0) deff=1.0-1.0*hist[6]/(1.0*vortices);
            delete[] hist;
            delete[] pixels;
        }
        if(matrix!=NULL) delete[] matrix;
        
        s=IntToStr(n)+"\t"+IntToStr(vortices)+"\t"+FloatToStr(deff)+"\t# "+inputfiles[n].c_str()+"\n";
        printf("%s",s.c_str());
        FileWrite(f,s.c_str(),s.length());
    }
    FileClose(f);
}



//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
void GLPP::print_info() {
    if(readdatatype!=0) {
        
    } else {
        printf("no data read.\n");
    }
}

void GLPP::loadparams(int n, char* argv[]) {
    int i;
    /*
     if(argv[1][0]=='-') {
     s=argv[1];
     if(s=="-i") option=1;
     else if(s=="-analysis") option=2;
     else if(s=="-analysish") option=3; //same a 2, but write a header line
     else option=9999;
     }
     */
    
    action=0;
    if(n>1) action=atoi(argv[1]);
    if(n>2) outputfileprefix=argv[2];
    numinputfiles=n-3;
    if(numinputfiles>0) {
        inputfiles=new string[numinputfiles];
        for(i=0;i<numinputfiles;i++) inputfiles[i]=argv[3+i];
    }
};

int GLPP::run() {
    Adata adat;
    string s;
    int n;
    int res=0;
    double *gradsq;
    
    
    
    if(action==1002) {
        if(numinputfiles>0) {
            printf("#time\tBx\tBy\tBz\t<|psi|^2>\t<|grad psi|^2>\tvfrac\tJs2_av\tJs2_min\tJs2_max\tJs2_stddev\n");
            for(n=0;n<numinputfiles;n++) {
                if(loadBDAT(inputfiles[n])==0) {
                    gradsq=new double[NN];
                    res=calc_current(gradsq);
                    if(res!=0) printf("##  current calc fail %d\n",res);
                    //setSV(90,80,15,30,30,30);
                    analysis(adat,gradsq,0.2);
                    delete[] gradsq;
                    printf("%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",adat.time,adat.Bx,adat.By,adat.Bz,adat.psisqfrac,adat.gradfrac,adat.vfrac,adat.Js_av,adat.Js_min,adat.Js_max,adat.Js_stddev);
                }
                resetdata();
            }
        }
    } else if(action==211) {
        if(numinputfiles>0) {
           for(n=0;n<numinputfiles;n++) {
               if(loadBDAT(inputfiles[n])==0) {
                   res=calc_current();
                   write_XMDF(n+1,numinputfiles,(readdatatype&0xFF)|WRITE_SC|WRITE_OP_amp2|WRITE_OP_ph|WRITE_OP_reim);
               } else
                   printf("problem loading BDAT file: %s\n",inputfiles[n].c_str());
               resetdata();
           }
        }
        } else if(action==3201) {
            delaunay_analysis();
    }
    else
        print_info();
    
    return res;
};



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if 0
int main(int argc, char* argv[])
{
    GLPP *gpp;
    int res;
    if(argc==1) printf("GL data post processor, version %d.%d\n\n",((GLPPversion>>16)&0xFF),((GLPPversion>>8)&0xFF));
    
    if(argc>1)
    {
        gpp=new GLPP();
        gpp->loadparams(argc,argv);
        res=gpp->run();
        if(res!=0) printf("#There was an error: %d\n",res);
        delete gpp;
  	}
    
    return 0;
};
#endif
