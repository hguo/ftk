#ifndef _GLGPU_IO_HELPER_H
#define _GLGPU_IO_HELPER_H

#include "GLHeader.h"
#include "BDATReader.h"

bool GLGPU_IO_Helper_ReadBDAT(
    const std::string& filename, 
    GLHeader &hdr,
    float **rho, float **phi, float **re, float **im, float **Jx, float **Jy, float **Jz,
    bool header_only=false, bool supercurrent=false);

bool GLGPU_IO_Helper_ReadLegacy(
    const std::string& filename, 
    GLHeader &hdr, 
    float **rho, float **phi, float **re, float **im, float **Jx, float **Jy, float **Jz,
    bool header_only=false, bool supercurrent=false);

void GLGPU_IO_Helper_ComputeSupercurrent(
    GLHeader &h, const float *re, const float *im, float **Jx, float **Jy, float **Jz);

bool GLGPU_IO_Helper_ReadNetCDF(
    const std::string& filename, 
    GLHeader &hdr, 
    float **psi); 

bool GLGPU_IO_Helper_WriteNetCDF(
    const std::string& filename, 
    GLHeader &hdr, 
    const float *rho, const float *phi, 
    const float *re, const float *im, const float *Jx, const float *Jy, const float *Jz);

bool GLGPU_IO_Helper_WriteRaw_rho_phi(
    const std::string& filename_rho, 
    const std::string& filename_phi, 
    GLHeader& hdr, 
    const float *rho, const float *phi);

#endif
