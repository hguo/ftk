//---------------------------------------------------------------------------
// data file IO
//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string>

#pragma hdrstop

#include "fileutils.h"
#include "stringutils.h"

#include "datafile.h"

//---------------------------------------------------------------------------

int BDATreader::open(string fn,bool makeTOC) {
    if(state==0) {
        char FID[8];
        unsigned int BOM;
        if(fn.length()>0) filename=fn;
        f=FileOpen(filename,fmOpenRead);
        if(f<0) return -2;
        size=FileSize(f);
        pos=0;
        FileRead(f,FID,4);FID[4]=0;pos+=4;
        FileRead(f,&BOM,4);pos+=4;
        if(string(FID)!="BDAT") return -3;
        
        //check for endian encoding of file data - not implemented/used
        switch (BOM) {
            case 0xFEFF:
                endian=0; //match
                break;
            default:
                endian=-1;
                
                break;
        }
        if(makeTOC) createTOC();
        state=1;
        lastinfo.records=0;
        
    }
    else return -1;
    return 0;
};

//---------------------------------------------------------------------------

int BDATreader::createTOC() {
    int res,laststate,n;
    long lastpos,p;
    datablockinfo db,saveddb;

    if(hasTOC) deleteTOC();
    if(state>0) {
        laststate=state;
        lastpos=pos;
        saveddb=lastinfo;
        state=1;pos=8;
        FileSeek(f,pos,SEEK_SET);
        res=getnextblockinfo(db);
        while(res==0) {
            TOCN++;
            res=skipblock();
            if(res==0) res=getnextblockinfo(db);
        }
        if(TOCN>0) {
            hasTOC=true;
            TOCpos=new long[TOCN];
            TOCID=new unsigned int[TOCN];
            TOCnames=new string[TOCN];
            state=1;p=pos=8;
            FileSeek(f,pos,SEEK_SET);
            res=getnextblockinfo(db);
            n=0;
            while(res==0) {
                res=skipblock();
                TOCpos[n]=p;
                TOCID[n]=lastinfo.ID;
                TOCnames[n]=lastinfo.name;
                p=pos;
                if(res==0) res=getnextblockinfo(db);
                n++;
            }
        }
        state=laststate;
        pos=lastpos;
        lastinfo=saveddb;
        FileSeek(f,pos,SEEK_SET);
    }
    else return -1;
    return 0;
};

void BDATreader::deleteTOC() {
    if(hasTOC) {
        delete[] TOCpos;
        delete[] TOCID;
        delete[] TOCnames;
        hasTOC=false;
        TOCN=0;
    }
};

int BDATreader::seek(string name,datablockinfo &db) {
    int n=0;
    if(!hasTOC) createTOC();
    if(!hasTOC) return -1; //no TOC available (empty file?)
    
    while((n<TOCN) && (TOCnames[n]!=name)) n++;
    if(n==TOCN) return -2; //entry not found
    
    pos=TOCpos[n];state=1;
    FileSeek(f,pos,SEEK_SET);
    
    return getnextblockinfo(db);
};

int BDATreader::seek(unsigned int ID,datablockinfo &db) {
    int n=0;
    if(!hasTOC) createTOC();
    if(!hasTOC) return -1; //no TOC available (empty file?)
    
    while((n<TOCN) && (TOCID[n]!=ID)) n++;
    if(n==TOCN) return -2; //entry not found
    
    pos=TOCpos[n];state=1;
    FileSeek(f,pos,SEEK_SET);
    
    return getnextblockinfo(db);
};

//---------------------------------------------------------------------------

int BDATreader::printlist() {
    int res,laststate,n;
    long lastpos;
    string tname;
    datablockinfo db,saveddb;
    if(state>0) {
        laststate=state;
        lastpos=pos;
        saveddb=lastinfo;
        state=1;pos=8;
        FileSeek(f,pos,SEEK_SET);
        res=getnextblockinfo(db);
        printf("BDAT binary data file information:\n - name %s\n - size: %ld\n",filename.c_str(),size);
        n=1;
        while(res==0) {
            printf(" - data block %d\n   - name: %s\n   - ID: %d\n   - type: (%08X) ",n,lastinfo.name.c_str(),lastinfo.ID,lastinfo.type);
            switch (lastinfo.type) {
                case DT_bit:tname="bits";break;
                case DT_byte:tname="byte";break;
                case DT_ubyte:tname="unsigned byte";break;
                case DT_word:tname="word/2-byte integer";break;
                case DT_uword:tname="unsigned word/2-byte integer";break;
                case DT_int:tname="int/4-byte integer";break;
                case DT_uint:tname="unsigned int/4-byte integer";break;
                case DT_long:tname="long int/8-byte integer";break;
                case DT_ulong:tname="unsigned long int/8-byte integer";break;
                case DT_single:tname="float/4-byte floating point number";break;
                case DT_double:tname="double/8-byte floating point number";break;
                case DT_quad:tname="quad/16-byte floating point number";break;
                case DT_csingle:tname="complex float/two 4-byte floating point numbers";break;
                case DT_cdouble:tname="complex double/two 8-byte floating point numbers";break;
                case DT_cquad:tname="complex quad/two 16-byte floating point number";break;
                case DT_user:tname="user defined format";break;
                default:tname="unknown";break;
            }
            printf("%s\n   - records: %d\n   - record size: %d bytes\n   - data size: %d bytes\n",tname.c_str(),lastinfo.records,lastinfo.recordsize,lastinfo.records*lastinfo.recordsize);
            res=skipblock();
            if(res==0) res=getnextblockinfo(db);
            n++;
        }
        state=laststate;
        pos=lastpos;
        lastinfo=saveddb;
        FileSeek(f,pos,SEEK_SET);
    }
    else return -1; //file needs to be open
    return 0;
};

//---------------------------------------------------------------------------

int BDATreader::getnextblockinfo(datablockinfo &db) {
    unsigned int a,l;
    char buf[256];
    if((state==1) || (state==3)) {
        if(pos>=size) {
            state=4;
            db.name="end of file";
            return -4;
        }
        FileRead(f,&a,4);pos+=4;
        l=a&0xFF;
        FileRead(f,buf,l);buf[l]=0;pos+=l;
        db.name=buf;
        db.ID=a>>8;
        FileRead(f,&db.type,4);pos+=4;
        FileRead(f,&db.records,4);pos+=4;
        FileRead(f,&db.recordsize,4);pos+=4;
        state=2;
        lastinfo=db;
    }
    else return -1;
    return 0;
};

//---------------------------------------------------------------------------

int BDATreader::skipblock() {
    unsigned int l;
    if(state==2) {
        l=lastinfo.records*lastinfo.recordsize;
        if(pos+l>size) return -4;
        FileSeek(f,l,SEEK_CUR);pos+=l;
        state=3;
    }
    else return -1;
    return 0;
};

int BDATreader::readblock(unsigned char *buffer,unsigned int maxlength) {
    unsigned int tl,l;
    if(state==2) {
        tl=lastinfo.records*lastinfo.recordsize;
        if(pos+tl>size) return -4;
        if((maxlength==0) || (maxlength>tl)) l=tl;
        else l=maxlength;
        FileRead(f,buffer,l);
        if(tl>l) FileSeek(f,tl-l,SEEK_CUR);
        pos+=tl;
        state=3;
    }
    else return -1;
    return 0;
};

int BDATreader::readblock(char *buffer,unsigned int maxlength) {
    return readblock((unsigned char *) &buffer[0],maxlength);
};

int BDATreader::readblock(int *buffer,unsigned int maxlength) {
    return readblock((unsigned char *) &buffer[0],maxlength*sizeof(int));
};

int BDATreader::readblock(unsigned int *buffer,unsigned int maxlength) {
    return readblock((unsigned char *) &buffer[0],maxlength*sizeof(unsigned int));
};

int BDATreader::readblock(float *buffer,unsigned int maxlength) {
    return readblock((unsigned char *) &buffer[0],maxlength*sizeof(float));
};

int BDATreader::readblock(double *buffer,unsigned int maxlength) {
    return readblock((unsigned char *) &buffer[0],maxlength*sizeof(double));
};

//---------------------------------------------------------------------------

//read and convert
//use these also for complex types with double length
//unfortunately cannot be done as template
int BDATreader::readconvert(int *buffer,unsigned int maxlength) {//converts read data into type int
    unsigned int tl,read,n,k,i,N,l,blocksize;
    unsigned char *tmp;
    if(state==2) {
        if((lastinfo.type==DT_int) || (lastinfo.type==DT_uint)) return readblock(buffer,maxlength); //no need for conversion
        tl=lastinfo.records*lastinfo.recordsize; //total size
        if(pos+tl>size) return -4;
        
        //we want to limited the temporary data needed for reading (max ~64MB) -> block-wise reading
        if(tl>TMPSIZE) {
            blocksize=TMPSIZE/lastinfo.recordsize;
            blocksize*=lastinfo.recordsize;
        }
        else blocksize=tl;
        tmp=new unsigned char[blocksize];
        n=0;
        N=lastinfo.records;if((maxlength>0) && (maxlength<N)) N=maxlength;
        while(tl>0) {
            read=(blocksize<tl?blocksize:tl);
            k=read/lastinfo.recordsize;
            FileRead(f,tmp,read);pos+=read;
            tl-=read;
            i=0;
            while((n<N) && (k>0)) {
                switch (lastinfo.type) {
                    case DT_byte:
                        buffer[n]=(int) tmp[i];i++;
                        break;
                    case DT_ubyte:
                        buffer[n]=(int) ((unsigned char *) tmp)[i];i++;
                        break;
                    case DT_word:
                        buffer[n]=(int) ((short *) tmp)[i];i++;
                        break;
                    case DT_uword:
                        buffer[n]=(int) ((unsigned short *) tmp)[i];i++;
                        break;
                    case DT_long:
                        buffer[n]=(int) ((long *) tmp)[i];i++;
                        break;
                    case DT_ulong:
                        buffer[n]=(int) ((unsigned long *) tmp)[i];i++;
                        break;
                    case DT_single:
                        buffer[n]=(int) ((float *) tmp)[i];i++; //floor to integer
                        break;
                    case DT_double:
                        buffer[n]=(int) ((double *) tmp)[i];i++; //floor to integer
                        break;
                    case DT_quad:
                        buffer[n]=(int) ((long double *) tmp)[i];i++; //floor to integer
                        break;
                    default: buffer[n]=0; //do not know how to convert any other type
                        break;
                }
                k--;
                n++;
            }
        }
        delete[] tmp;
        state=3;
    }
    else return -1;
    return 0;
};    

int BDATreader::readconvert(float *buffer,unsigned int maxlength) {//converts read data into type float
    unsigned int tl,read,n,k,i,N,l,blocksize;
    unsigned char *tmp;
    if(state==2) {
        if(lastinfo.type==DT_single) return readblock(buffer,maxlength); //no need for conversion
        tl=lastinfo.records*lastinfo.recordsize; //total size
        //we want to limited the temporary data needed for reading (max ~64MB) -> block-wise reading
        if(tl>TMPSIZE) {
            blocksize=TMPSIZE/lastinfo.recordsize;
            blocksize*=lastinfo.recordsize;
        }
        else blocksize=tl;
        tmp=new unsigned char[blocksize];
        n=0;
        N=lastinfo.records;if((maxlength>0) && (maxlength<N)) N=maxlength;
        while(tl>0) {
            read=(blocksize<tl?blocksize:tl);
            k=read/lastinfo.recordsize;
            FileRead(f,tmp,read);pos+=read;
            tl-=read;
            i=0;
            while((n<N) && (k>0)) {
                switch (lastinfo.type) {
                    case DT_byte:case DT_ubyte:
                        buffer[n]=(float) tmp[i];i++;
                        break;
                    case DT_word:case DT_uword:
                        buffer[n]=(float) ((short *) tmp)[i];i++;
                        break;
                    case DT_int:case DT_uint:
                        buffer[n]=(float) ((int *) tmp)[i];i++;
                        break;
                    case DT_long:case DT_ulong:
                        buffer[n]=(float) ((long *) tmp)[i];i++;
                        break;
                    case DT_double:
                        buffer[n]=(float) ((double *) tmp)[i];i++;
                        break;
                    case DT_quad:
                        buffer[n]=(float) ((long double *) tmp)[i];i++;
                        break;
                    default: buffer[n]=0.0; //do not know how to convert any other type
                        break;
                }
                k--;
                n++;
            }
        }
        delete[] tmp;
        state=3;
    }
    else return -1;
    return 0;
};  

int BDATreader::readconvert(double *buffer,unsigned int maxlength) {//converts read data into type double
    unsigned int tl,read,n,k,i,N,l,blocksize;
    unsigned char *tmp;
    if(state==2) {
        if(lastinfo.type==DT_double) return readblock(buffer,maxlength); //no need for conversion
        tl=lastinfo.records*lastinfo.recordsize; //total size
        //we want to limited the temporary data needed for reading (max ~64MB) -> block-wise reading
        if(tl>TMPSIZE) {
            blocksize=TMPSIZE/lastinfo.recordsize;
            blocksize*=lastinfo.recordsize;
        }
        else blocksize=tl;
        tmp=new unsigned char[blocksize];
        n=0;
        N=lastinfo.records;if((maxlength>0) && (maxlength<N)) N=maxlength;
        while(tl>0) {
            read=(blocksize<tl?blocksize:tl);
            k=read/lastinfo.recordsize;
            FileRead(f,tmp,read);pos+=read;
            tl-=read;
            i=0;
            while((n<N) && (k>0)) {
                switch (lastinfo.type) {
                    case DT_byte:case DT_ubyte:
                        buffer[n]=(double) tmp[i];i++;
                        break;
                    case DT_word:case DT_uword:
                        buffer[n]=(double) ((short *) tmp)[i];i++;
                        break;
                    case DT_int:case DT_uint:
                        buffer[n]=(double) ((int *) tmp)[i];i++;
                        break;
                    case DT_long:case DT_ulong:
                        buffer[n]=(double) ((long *) tmp)[i];i++;
                        break;
                    case DT_single:
                        buffer[n]=(double) ((float *) tmp)[i];i++;
                        break;
                    case DT_quad:
                        buffer[n]=(double) ((long double *) tmp)[i];i++;
                        break;
                    default: buffer[n]=0.0; //do not know how to convert any other type
                        break;
                }
                k--;
                n++;
            }
        }
        delete[] tmp;
        state=3;
    }
    else return -1;
    return 0;
};

//---------------------------------------------------------------------------

void BDATreader::primitives_warning() {
    if(sizeof(short)!=2) printf("*** BDATreader warning: 2-byte integer (short) conversion will fail!\n");
    if(sizeof(int)!=4) printf("*** BDATreader warning: 4-byte integer (int) conversion will fail!\n");
    if(sizeof(long)!=8) printf("*** BDATreader warning: 8-byte integer (long) conversion will fail!\n");
    if(sizeof(float)!=4) printf("*** BDATreader warning: single precision floating point (foat) conversion will fail!\n");
    if(sizeof(double)!=8) printf("*** BDATreader warning: double precision floating point (double) conversion will fail!\n");
    //this could fail anyway... (depending on the actual precision: 80 bit,96 bit, 106bit, 128bit)
    if(sizeof(long double)!=16) printf("*** BDATreader warning: quad precision floating point (long double, __float128) conversion will fail!\n");
    
    if(endian==-1) printf("*** BDATreader warning: endian mismatch detected - read data might have wrong order!\n");
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int BDATwriter::open(string fn,bool append) {
    char FID[8];
    unsigned int BOM;
    if(state==0) {
        if(fn.length()>0) filename=fn;
        if(!append) {
            size=0;
            f=FileCreate(filename);
            if(f<0) return -2;
            FID[0]='B';FID[1]='D';FID[2]='A';FID[3]='T';
            FileWrite(f,&FID,4);size+=4;
            BOM=0xFEFF;
            FileWrite(f,&BOM,4);size+=4;
            
            state=1;
        }
        else {
            f=FileOpen(filename,fmOpenRead);
            if(f<0) return -2;
            size=FileSize(f);
            FileRead(f,FID,4);FID[4]=0;
            FileRead(f,&BOM,4);
            if(string(FID)!="BDAT") return -3;
            FileClose(f);
            
            f=FileOpen(fn,fmOpenAppend);
            state=1;
        }
    }
    else return -1;
    return 0;
};

//---------------------------------------------------------------------------

int BDATwriter::writeblock(string name,unsigned int ID,unsigned int type,unsigned int records,unsigned int recordsize,char *buffer) {
    unsigned int a,l;
    if(state==1) {
        l=name.length();
        if(l>255) l=255;
        a=l+(ID<<8);
        FileWrite(f,&a,4);
        FileWrite(f,name.c_str(),l);size+=4+l;
        FileWrite(f,&type,4);size+=4;
        FileWrite(f,&records,4);size+=4;
        FileWrite(f,&recordsize,4);size+=4;
        l=records*recordsize;
        FileWrite(f,buffer,l);size+=l;
    }
    else return -1;
    return 0;
};

int BDATwriter::writeblock(datablockinfo &db,char *buffer) {
    return writeblock(db.name,db.ID,db.type,db.records,db.recordsize,buffer);
};

int BDATwriter::write(string name,unsigned int ID,string str) {
    return writeblock(name,ID,DT_byte,str.length(),1,(char *) str.c_str());
};

int BDATwriter::write(string name,unsigned int ID,int *buffer,unsigned int records) {
    return writeblock(name,ID,DT_int,records,4,(char *) buffer);
};

int BDATwriter::write(string name,unsigned int ID,unsigned int *buffer,unsigned int records) {
    return writeblock(name,ID,DT_uint,records,4,(char *) buffer);
};

int BDATwriter::write(string name,unsigned int ID,long *buffer,unsigned int records) {
    return writeblock(name,ID,DT_long,records,8,(char *) buffer);
};

int BDATwriter::write(string name,unsigned int ID,unsigned long *buffer,unsigned int records) {
    return writeblock(name,ID,DT_ulong,records,8,(char *) buffer);
};

int BDATwriter::write(string name,unsigned int ID,float *buffer,unsigned int records) {
    return writeblock(name,ID,DT_single,records,4,(char *) buffer);
};

int BDATwriter::write(string name,unsigned int ID,double *buffer,unsigned int records) {
    return writeblock(name,ID,DT_double,records,8,(char *) buffer);
};

//---------------------------------------------------------------------------

int BDATwriter::splitwrite(string namere,string nameim,unsigned int IDre,unsigned int IDim,float *buffer,unsigned int records) {
    float *tmp;
    unsigned int a,l,n,k,N,type,recordsize,blocklen;
    type=DT_single;
    recordsize=4;
    if(state==1) {
        l=namere.length();
        if(l>255) l=255;
        a=l+(IDre<<8);
        FileWrite(f,&a,4);
        FileWrite(f,namere.c_str(),l);size+=4+l;
        FileWrite(f,&type,4);size+=4;
        FileWrite(f,&records,4);size+=4;
        FileWrite(f,&recordsize,4);size+=4;
        
        l=records*recordsize;
        if(l>TMPSIZE) blocklen=TMPSIZE/recordsize;
        else blocklen=records;
        tmp=new float[blocklen];
        N=records;
        k=0;
        while(N>0) {
            n=0;
            while((N>0) && (n<blocklen)) {tmp[n]=buffer[k];n++;k+=2;N--;}
            l=4*n;
            FileWrite(f,tmp,l);size+=l;
        }
        
        //second block, the imaginary numbers
        l=nameim.length();
        if(l>255) l=255;
        a=l+(IDim<<8);
        FileWrite(f,&a,4);
        FileWrite(f,nameim.c_str(),l);size+=4+l;
        FileWrite(f,&type,4);size+=4;
        FileWrite(f,&records,4);size+=4;
        FileWrite(f,&recordsize,4);size+=4;
        N=records;
        k=1;
        while(N>0) {
            n=0;
            while((N>0) && (n<blocklen)) {tmp[n]=buffer[k];n++;k+=2;N--;}
            l=4*n;
            FileWrite(f,tmp,l);size+=l;
        }
        delete[] tmp;
    }
    else return -1;
    return 0;
};

int BDATwriter::splitwrite(string namere,string nameim,unsigned int IDre,unsigned int IDim,double *buffer,unsigned int records) {
    double *tmp;
    unsigned int a,l,n,k,N,type,recordsize,blocklen;
    type=DT_double;
    recordsize=8;
    if(state==1) {
        l=namere.length();
        if(l>255) l=255;
        a=l+(IDre<<8);
        FileWrite(f,&a,4);
        FileWrite(f,namere.c_str(),l);size+=4+l;
        FileWrite(f,&type,4);size+=4;
        FileWrite(f,&records,4);size+=4;
        FileWrite(f,&recordsize,4);size+=4;
        
        l=records*recordsize;
        if(l>TMPSIZE) blocklen=TMPSIZE/recordsize;
        else blocklen=records;
        tmp=new double[blocklen];
        N=records;
        k=0;
        while(N>0) {
            n=0;
            while((N>0) && (n<blocklen)) {tmp[n]=buffer[k];n++;k+=2;N--;}
            l=4*n;
            FileWrite(f,tmp,l);size+=l;
        }
        
        //second block, the imaginary numbers
        l=nameim.length();
        if(l>255) l=255;
        a=l+(IDim<<8);
        FileWrite(f,&a,4);
        FileWrite(f,nameim.c_str(),l);size+=4+l;
        FileWrite(f,&type,4);size+=4;
        FileWrite(f,&records,4);size+=4;
        FileWrite(f,&recordsize,4);size+=4;
        N=records;
        k=1;
        while(N>0) {
            n=0;
            while((N>0) && (n<blocklen)) {tmp[n]=buffer[k];n++;k+=2;N--;}
            l=4*n;
            FileWrite(f,tmp,l);size+=l;
        }
        delete[] tmp;
    }
    else return -1;
    return 0;
};


//---------------------------------------------------------------------------

