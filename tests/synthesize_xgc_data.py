#!/usr/bin/env python

# modified script based on https://github.com/rmchurch/synthetic_blobs
print('executing xgc blob synthesizer...')

import os
import h5py
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation,LinearTriInterpolator

def interpolate_fieldLineFollow(Lstart,phiEnd,Binterp,dphi = None):

    #define RHS of ODE system of equations, dy/dt = f(y,t)
    def f(L,phi,Binterp):
        R=L[0]
        Z=L[1]
        B=Binterp(R,Z)
        BR=B[0]
        BZ=B[1]
        Bphi=B[2]
        #model equations
        f0 = R*BR/Bphi
        f1 = R*BZ/Bphi
        #f2 = 1.
        return [f0,f1]#,f2]

    #create an array of phi coordinates for which the particle position
    #will be calculated, in between the initial and end phi poisitions
    Npts = 1000
    if dphi is not None:
        Npts = int(np.abs((phiEnd-Lstart[2])/dphi)) + 1
        phi = Lstart[2] + np.sign(phiEnd-Lstart[2])*np.arange(Npts)*dphi
    else:
        phi = np.linspace(Lstart[2],phiEnd,Npts)

    soln = odeint(f,Lstart[0:2],phi,args=(Binterp,))
    Lout = np.hstack((soln,phi[:,np.newaxis]))
    return Lout


class syntheticBlobs():
    def __init__(self,RZ,psin,tri,Bgrid,sml_nphi):
        self.RZ = RZ
        self.R0,self.Z0 = RZ[0,:]
        self.psin = psin
        self.tri = tri
        self.triObj = Triangulation(RZ[:,0],RZ[:,1],tri)
        self.theta = self.calc_theta(RZ[:,0], RZ[:,1])
        ##find x-point, to exclude from interpolator
        #Bpol = np.sqrt(np.sum(Bgrid[:,0:2]**2.,axis=1))
        #ind = np.argmin(Bpol[10:])+10
        #eq_x_r,eq_x_z = RZ[ind,:]
        #goodinds = ~((psin>=1) | ((psin<=1) & (RZ[:,1]>eq_x_z)))
        #mask = np.all(goodinds[tri],axis=1)
        #self.triObj_psintheta = Triangulation(self.psin,self.theta,self.tri,mask=mask)
        #self.fRZ2psin = LinearTriInterpolator(self.triObj_psintheta,self.RZ[:,0]) 
        #self.fpsintheta2Z = LinearTriInterpolator(self.triObj_psintheta,self.RZ[:,1]) 
        self.fRZ2psin = LinearTriInterpolator(self.triObj,self.psin)
        self.Binterp = LinearNDInterpolator(RZ, Bgrid, fill_value = np.inf)
        self.sml_nphi = sml_nphi
        
    def psintheta2RZ(self,psin,theta):
        return self.fpsinitheta2R(psin,theta),self.fpsintheta2Z(psin,theta)
    
    def RZ2psin(self,R,Z):
        return self.fRZ2psin(R,Z)
             
    def calc_theta(self,R,Z):
        """Calculate poloidal angle, with 0 deg at LFS midplane"""
        return np.arctan2(Z-self.Z0, R-self.R0)
    
    def generate(self,xcenter,ntor,Lpol,Lrad,dnOvern,use_RZ=True):
        """ Generate a blob from the input characteristics
        xcenter [3]: Blob center coordinates (psin, theta, phi)
        ntor [1]: Toroidal mode number of the blob
        Lpol [1]: Blob diameter in poloidal direction
        Lrad [1]: Blob diameter in radial direction
        dnOvern [1]: Scalar of the magnitude of the blob at center, dn/n
        use_RZ [bool]: input xcenter with (R,Z,phi) or (psin,theta,phi)
        """
        
        if use_RZ:
            R1 = xcenter[0]; Z1 = xcenter[1]; phi1 = xcenter[2]
            psin1 = self.fRZ2psin(R1,Z1)
            theta1 = self.calc_theta(R1,Z1)
        else:
            raise ValueError("use_RZ==False not implemented yet")
            psin1 = xcenter[0]
            theta1 = xcenter[1]
            phi1 = xcenter[2]
            R1,Z1 = self.psintheta2RZ(psin1,theta1)
        
        #force quantized phi1
        dphi = 2*np.pi/self.sml_nphi
        phis = np.arange(self.sml_nphi)*dphi
        phiInd = int(np.round(phi1 / dphi) % self.sml_nphi)
        phi1 = phis[phiInd]

        #assume toridal mode number ntor = 2*pi/lambda_tor, lambda_tor the toroidal wavelength
        dphiEnd = 2*np.pi/ntor*R1/self.R0/2.
        Lstart = np.array([R1,Z1,phi1])
        
        #generate field-line path 
        LoutFwd = interpolate_fieldLineFollow(Lstart, phi1+dphiEnd,Binterp,dphi=dphi) 
        LoutBwd = interpolate_fieldLineFollow(Lstart, phi1-dphiEnd,Binterp,dphi=dphi) 
        Lout = np.concatenate((LoutBwd[1:,:][::-1,:],LoutFwd) ) #remove duplicate point, concatenate
        phioutInds = (np.round(Lout[:,2] / dphi) % self.sml_nphi).astype(int)
        
        tmp = np.sin(np.arange(Lout.shape[0])/(Lout.shape[0]-1)*np.pi) 
        dn_par = np.zeros((self.sml_nphi,))
        dn_par[phioutInds] = tmp
        #loop through toroidal planes, interpolate onto XGC R,Z mesh, witha  cutoff of 3*sigma
        #interpolate onto the phi XGC grid
        #(wont be needed if using dphi input to interpolate_fieldLineFollow)
        
        #loop through toroidal planes, interpolate onto XGC R,Z mesh, witha  cutoff of 3*sigma
        Bfield1 = self.Binterp(R1,Z1)
        Bpol1 = np.sqrt(np.sum(Bfield1[0:2]**2.))
        B1 = np.sqrt(np.sum(Bfield1**2.))
        alpha1 = np.arccos(Bfield1[0]/Bpol1) 
        
        dnXGC = np.zeros((self.sml_nphi,RZ.shape[0]))
        R = self.RZ[:,0]; Z = self.RZ[:,1]
        for p,phip in enumerate(phioutInds):
            Rp = Lout[p,0]; Zp = Lout[p,1]
            #first, adjust blob size in radial and poloidal directions, based on flux expansion
            Bfieldp = self.Binterp(Lout[p,0],Lout[p,1])
            Bpolp = np.sqrt(np.sum(Bfieldp[0:2]**2.))
            Bp = np.sqrt(np.sum(Bfieldp**2.))
            Lradp = Lrad*(Rp*Bpolp)/(R1*Bpol1)
            Lpolp = Lpol*B1/Bp*Lrad/Lradp 
            #adjust the angle
            alphap = np.arccos(Bfieldp[0]/Bpolp)
            alpha = alpha1 - alphap
            dnXGC[phip,:] = dnOvern*dn_par[phip]*np.exp(-(((R-Rp)*np.cos(alpha) + (Z-Zp)*np.sin(alpha))/Lrad)**2 + -(((R-Rp)*np.sin(alpha) - (Z-Zp)*np.cos(alpha))/Lpol)**2)

        return dnXGC

fileDir = os.getenv('FTK_XGC_TEST_DATA_PATH') + "/"
print('xgc_file_dir=', fileDir)
fileBfield = fileDir + 'xgc.bfield.h5'
fb = h5py.File(fileBfield,'r')
Bgrid = fb['node_data[0]/values'][:]

fileMesh = fileDir + 'xgc.mesh.h5'
fm = h5py.File(fileMesh,'r')
RZ = fm['coordinates/values'][:]
#you may have to replace this with a hardcoded value from units.m
try:
    fileEq = fileDir + 'xgc.equil.h5'
    feq = h5py.File(fileEq,'r')
    psi_x = feq['eq_psi_x'][...]
except:
    feq = open(fileDir+'units.m')
    for line in feq:
        if 'psi_x' in line:
            psi_x = float(line.split()[1])

psin = fm['psi'][:]/psi_x
tri= fm['cell_set[0]/node_connect_list'][...]
triObj = Triangulation(RZ[:,0],RZ[:,1],tri)
file3d = fileDir + 'xgc.3d.00001.h5'
f3d = h5py.File(file3d,'r')
sml_nphi = f3d['nphi'][0]
sml_iphi = f3d['iphi'][0]

#setup Bfield interpolator (could use higher order interpolation scheme)
Binterp = LinearNDInterpolator(RZ, Bgrid, fill_value = np.inf)

blob_generator = syntheticBlobs(RZ,psin,tri,Bgrid,sml_nphi)




#now generate some blobs
#xcenter = np.array([0.95,0,0]) #(psin,theta,phi)
#xcenter = np.array([2.26,0,0]) #(R,Z,phi)

for timestep in range (0, 5):
    xcenter = np.array([2.26, timestep * 0.01, 0]) #(R,Z,phi)
    ntor = 5
    Lpol = 0.01
    Lrad = Lpol/1.5
    dnOvernMag = 0.1
    dnOvernXGC = blob_generator.generate(xcenter,ntor,Lpol,Lrad,dnOvernMag)

    print("synthesizing xgc timestep", timestep) #dnOvernXGC.shape)

    file_output = 'xgc.synthetic.%04d.h5' % timestep
    fo = h5py.File(file_output,'w')
    fo['/dnOvernXGC'] = dnOvernXGC.transpose()
    fo['/nphi'] = sml_nphi
    fo['/iphi'] = sml_iphi
    fo.close()
