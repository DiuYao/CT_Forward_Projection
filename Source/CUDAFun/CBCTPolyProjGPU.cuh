/*
* GPU�������õ�դ���ݺ�̽������Ӧ���ݶ��Ǽ�����ϵ��
* ͸��ϵ������������ϵ��
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


/* ************ ��ʼ�� ************ */
void initDeviceGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoResponseGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);

void initDeviceNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);
void initDeviceSinMatNoResponseNoGrid(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj & h_mPolyForwardProj);


/* GPU�ڴ����� */
// ̽������Ӧ
void mallocDetResponse(PolyForwardProj & d_mPolyForwardProj, CTScanParas mCTScanParas);


// GPU�ڴ��ͷ�
void freeDeviceMemory(PolyForwardProj & d_mPolyForwardProj, Coordinate & d_mCoordinate);


// ----------------------------------------------------------------------------------------
// ---------------------------------Kernel function----------------------------------------
__global__ void computeIntPointCoordinatesKernel(Coordinate d_Coordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

__global__ void transformKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, float theta);

// ���㵥�����ܶ�ͶӰ�ĺ˺���
__global__ void computeSinMatIndensityProjKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);

/* *************** Grid ********************* */ 
// ��դ�Ķ���ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);


// ��դ�ĵ����ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);
__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);


// ��դ�ĵ����ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardSinMatNoResponseProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

/* ************** No Grid *************/
// ��դ�Ķ���ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);

// ��դ�ĵ����ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num);
__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);

// ��դ�ĵ����ʶ�ɫ��̽Ԫ��Ӧ����Ͷ�˺���
__global__ void forwardSinMatNoResponseProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas);


// -----------------------------------------------------------------------------------

// ��ʼ������
void initTexture3D(cudaTextureObject_t& texObj, cudaArray_t &d_cuArray3D, float* h_data, cudaExtent volumeSize);

// ��������, ����cudaArray
void updateTex(cudaArray_t &d_cuArray3D, float* h_data, cudaExtent volumeSize);




