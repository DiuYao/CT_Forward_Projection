/*
* GPU程序所用的栅数据和探测器响应数据都是计算后的系数
* 透过系数和吸收吸收系数
*/

#include "device_launch_parameters.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <cmath>

#include "..\DataType.h"


#define CUDAFREE(varP)\
 if(varP != nullptr) \
{ \
cudaFree(varP); \
 varP = nullptr;\
}

#define PI acosf(-1)
#define BLOCKSIZEX 8
#define BLOCKSIZEY 8
#define BLOCKSIZEZ 8



struct PolyForwardProj; 
struct Coordinate;
struct CTScanSystemInfo;
struct CTScanParas;

void forwardProjGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void forwardSinMatProjGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void forwardSinMatNoResponseProjGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);


void forwardProjNoGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void forwardSinMatProjNoGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void forwardSinMatNoResponseProjNoGridGPU(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);


/* ************ 初始化 ************ */
void initDeviceGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoResponseGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);

void initDeviceNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoResponseNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);


/* GPU内存申请 */
// 探测器响应
void mallocDetResponse(PolyForwardProj & d_mPolyForwardProj, CTScanParas mCTScanParas);


// GPU内存释放
void freeDeviceMemory(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate);


// ----------------------------------------------------------------------------------------
// ---------------------------------Kernel function----------------------------------------
__global__ void computeIntPointCoordinatesKernel(Coordinate d_Coordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

__global__ void transformKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, float theta);

// 计算单材质密度投影的核函数
__global__ void computeSinMatIndensityProjKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);

/* *************** Grid ********************* */ 
// 带栅的多材质多色有探元响应的正投核函数
__global__ void forwardProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);


// 带栅的单材质多色有探元响应的正投核函数
__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);
__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);


// 带栅的单材质多色无探元响应的正投核函数
__global__ void forwardSinMatNoResponseProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

/* ************** No Grid *************/
// 无栅的多材质多色有探元响应的正投核函数
__global__ void forwardProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);

// 带栅的单材质多色有探元响应的正投核函数
__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);
__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

// 无栅的单材质多色无探元响应的正投核函数
__global__ void forwardSinMatNoResponseProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);


// -----------------------------------------------------------------------------------

// 初始化纹理
void initTexture3D(cudaTextureObject_t& texObj, cudaArray_t &d_cuArray3D, float* h_data, cudaExtent volumeSize);

// 更新纹理, 更新cudaArray
void updateTex(cudaArray_t &d_cuArray3D, float* h_data, cudaExtent volumeSize);




