

#include "CBCTPolyProjGPU.cuh"

#include <iostream>


cudaTextureObject_t texObj = 0;
cudaArray_t d_cuArray3D;
cudaExtent volumeSize;

void forwardProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// Select Gpu
	/*cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	}*/

	//initDeviceVar(d_mPolyForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mPolyForwardProj);

	//cudaTextureObject_t texObj = 0;
	// ģ���������
	updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// դ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);
	// ̽������Ӧ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);


	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		forwardProjGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);

	//// ���Դ���
	//FILE* fp;
	//fp = fopen("test.raw", "wb");
	//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
	//fclose(fp);
}


//// ��դ��������ͶӰ
//void forwardSinMatProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
//{
//	//// ģ���������
//	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);  // �����ܶ�
//
//	// դ���� ���� -> �Դ�
//	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);
//	
//	// ̽������Ӧ���� ���� -> �Դ�
//	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);
//
//
//	// GPU��ʱ
//	cudaEvent_t g_start, g_stop;
//	cudaEventCreate(&g_start);
//	cudaEventCreate(&g_stop);
//	cudaEventRecord(g_start, 0);
//
//	// Kernel parameters
//	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
//	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);
//
//	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
//	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);
//
//	float angle = 0.0f;
//	for (size_t i = 0; i < mCTScanParas.projNum; i++)
//	{
//
//		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������
//
//		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);
//
//		cudaError_t cudaStatus = cudaGetLastError();
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//		}cudaStatus = cudaDeviceSynchronize();
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
//		}
//
//		//// ���Դ���
//		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
//		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
//		//FILE* fp;
//		//fp = fopen("test.raw", "wb");
//		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
//		//fclose(fp);
//
//		forwardSinMatProjGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);
//
//		cudaStatus = cudaGetLastError();
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//		}cudaStatus = cudaDeviceSynchronize();
//		if (cudaStatus != cudaSuccess) {
//			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
//		}
//
//		//// ���Դ���
//		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//		//// ���
//		////FILE* fp;
//		//fp = fopen("test.raw", "wb");
//		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
//		//fclose(fp);
//	}
//
//	// I  �Դ� -> ����
//	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//
//
//	// ��ʱ
//	cudaEventRecord(g_stop, 0);
//	cudaEventSynchronize(g_stop);
//	float elapsedTime = 0;
//	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
//	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
//	cudaEventDestroy(g_start);
//	cudaEventDestroy(g_stop);
//
//}

// ��դ��������ͶӰ
void forwardSinMatProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	//// ģ���������
	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);  // �����ܶ�

	// դ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// ̽������Ӧ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeProj(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1, (mCTScanParas.projNum - 1) / blockSizeProj.y + 1);

	forwardSinMatProjGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas);

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "forwardSinMatProjGridKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardSinMatProjGridKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);

}

void forwardSinMatNoResponseProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// դ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// ̽������Ӧ���� ���� -> �Դ�
	// cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeProj(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1, (mCTScanParas.projNum - 1) / blockSizeProj.y + 1);

	forwardSinMatNoResponseProjGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas);

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "forwardSinMatNoResponseProjGridKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardSinMatNoResponseProjGridKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}


// ��դ��ͶӰ
void forwardProjNoGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	//// Select Gpu
	//cudaError_t cudaStatus;
	//cudaStatus = cudaSetDevice(0);
	//if (cudaStatus != cudaSuccess) {
	//	fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	//}

	//initDeviceVar(d_mPolyForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mPolyForwardProj);

	//cudaTextureObject_t texObj = 0;
	// ģ���������
	updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// ̽������Ӧ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);


	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		forwardProjNoGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}


// ��դ��������ͶӰ
void forwardSinMatProjNoGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{


	//// ģ���������
	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// ̽������Ӧ���� ���� -> �Դ�
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeProj(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1, (mCTScanParas.projNum - 1) / blockSizeProj.y + 1);

	forwardSinMatProjNoGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas);

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "forwardSinMatProjNoGridKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardSinMatProjNoGridKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}

void forwardSinMatNoResponseProjNoGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// դ���� ���� -> �Դ�
	//cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// ̽������Ӧ���� ���� -> �Դ�
	// cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeProj(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1, (mCTScanParas.projNum - 1) / blockSizeProj.y + 1);

	forwardSinMatNoResponseProjNoGridKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas);

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "forwardSinMatNoResponseProjNoGridKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardSinMatNoResponseProjNoGridKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// I  �Դ� -> ����
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}



void initDeviceGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}
	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
}

void initDeviceSinMatGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}


	}

	// ���Դ���
	//h_mPolyForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	//cudaMemcpy(h_mPolyForwardProj.proj, d_mPolyForwardProj.proj, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	//// ���
	//FILE* fp;
	//fp = fopen("test.raw", "wb");
	//fwrite(h_mPolyForwardProj.proj, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
	//fclose(fp);

	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
}

// ��ʼ���������ʡ���դ������Ӧ���н���
void initDeviceSinMatFoSpSiGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	// ������ֵ�����
	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ����ƫ�����ݷ���
	size_t sizeFoSpOffset = mCTScanParas.projNum * sizeof(float);

	cudaStatus = cudaMalloc(&d_mPolyForwardProj.foSpOffsetU, sizeFoSpOffset);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetU cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mPolyForwardProj.foSpOffsetV, sizeFoSpOffset);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetV cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ��������
	cudaStatus = cudaMemcpy(d_mPolyForwardProj.foSpOffsetU, h_mPolyForwardProj.foSpOffsetU, sizeFoSpOffset, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetU cudaMemcpy failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMemcpy(d_mPolyForwardProj.foSpOffsetV, h_mPolyForwardProj.foSpOffsetV, sizeFoSpOffset, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetV cudaMemcpy failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		//angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformFocalSpotKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, i);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformFocalSpotKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformFocalSpotKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// ����ͶӰ����ƫ�ƣ���������ƫ�ƣ��൱�������̽����ͬʱ��ͬһ����ƫ��
	createTexture3D(texObj, d_mPolyForwardProj.proj, mCTScanParas.dNumU, mCTScanParas.dNumV, mCTScanParas.projNum);

	dim3 blockSizePOffsetMatch(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizePOffsetMatch((mCTScanParas.dNumU - 1) / blockSizePOffsetMatch.x + 1, (mCTScanParas.dNumV - 1) / blockSizePOffsetMatch.y + 1, (mCTScanParas.projNum - 1) / blockSizePOffsetMatch.z + 1);

	projOffsetMatchKernel << <gridSizePOffsetMatch, blockSizePOffsetMatch >> > (d_mPolyForwardProj, texObj, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "projOffsetMatchKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching projOffsetMatchKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));

}

void initDeviceSinMatNoResponseGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

}

void initDeviceNoGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// ��ʼ��device
	cudaDeviceSynchronize();
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();

	// Select Gpu
	//cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	}

	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed!");
	}
	cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Variable d_mPolyForwardProj->phantom cudaMalloc failed!");
	}
	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// ̽������Ӧ����
	//mallocDetResponse(d_mPolyForwardProj, mCTScanParas);
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
}


void initDeviceSinMatNoGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));


}

void initDeviceSinMatFoSpSiNoGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ����ƫ�����ݷ���
	size_t sizeFoSpOffset = mCTScanParas.projNum * sizeof(float);

	cudaStatus = cudaMalloc(&d_mPolyForwardProj.foSpOffsetU, sizeFoSpOffset);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetU cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mPolyForwardProj.foSpOffsetV, sizeFoSpOffset);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetV cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ��������
	cudaStatus = cudaMemcpy(d_mPolyForwardProj.foSpOffsetU, h_mPolyForwardProj.foSpOffsetU, sizeFoSpOffset, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetU cudaMemcpy failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMemcpy(d_mPolyForwardProj.foSpOffsetV, h_mPolyForwardProj.foSpOffsetV, sizeFoSpOffset, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.foSpOffsetV cudaMemcpy failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		//angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformFocalSpotKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, i);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformFocalSpotKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformFocalSpotKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// ����ͶӰ����ƫ�ƣ���������ƫ�ƣ��൱�������̽����ͬʱ��ͬһ����ƫ��
	createTexture3D(texObj, d_mPolyForwardProj.proj, mCTScanParas.dNumU, mCTScanParas.dNumV, mCTScanParas.projNum);

	dim3 blockSizePOffsetMatch(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizePOffsetMatch((mCTScanParas.dNumU - 1) / blockSizePOffsetMatch.x + 1, (mCTScanParas.dNumV - 1) / blockSizePOffsetMatch.y + 1, (mCTScanParas.projNum - 1) / blockSizePOffsetMatch.z + 1);

	projOffsetMatchKernel << <gridSizePOffsetMatch, blockSizePOffsetMatch >> > (d_mPolyForwardProj, texObj, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "projOffsetMatchKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching projOffsetMatchKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));

}

void initDeviceSinMatNoResponseNoGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// ��ʼ��device
	cudaDeviceSynchronize();    // ȷ����CUDA��������
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // ������


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// ���ֵ�����
	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntX, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->x cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntY, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->y cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	cudaStatus = cudaMalloc(&d_mCoordinate.imgIntZ, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable imageCoordinate->z cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	/*cudaMalloc(&d_mCoordinate->detU, mCTScanParas.dNumU * sizeof(float));
	cudaMalloc(&d_mCoordinate->detV, mCTScanParas.dNumV * sizeof(float));*/

	//cudaMallocManaged();

	dim3 blockSizeIPC(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeIPC((mCTScanSystemInfo.intNum - 1) / blockSizeIPC.x + 10, (mCTScanParas.dNumU - 1) / blockSizeIPC.y + 10, (mCTScanParas.dNumV - 1) / blockSizeIPC.z + 10);

	computeIntPointCoordinatesKernel << <gridSizeIPC, blockSizeIPC >> > (d_mCoordinate, mCTScanSystemInfo, mCTScanParas);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "computeIntPointCoordinatesKernel launch failed: %s! \n\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching computeIntPointCoordinatesKernel!\n Error: %s!\n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// ��ʼ������
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // �洢���ֵ㴦����ֵ��ת�������
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// ���� �ܶ�ͼ���Proj
	// �ܶ�ͼ���Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU��ʱ
	cudaEvent_t g_start, g_stop;
	cudaEventCreate(&g_start);
	cudaEventCreate(&g_stop);
	cudaEventRecord(g_start, 0);

	// Kernel parameters
	dim3 blockSizeTrans(BLOCKSIZEX, BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeTrans((mCTScanSystemInfo.intNum - 1) / blockSizeTrans.x + 1, (mCTScanParas.dNumU - 1) / blockSizeTrans.y + 1, (mCTScanParas.dNumV - 1) / blockSizeTrans.z + 1);

	dim3 blockSizeProj(BLOCKSIZEY, BLOCKSIZEZ);
	dim3 gridSizeProj((mCTScanParas.dNumU - 1) / blockSizeProj.x + 1, (mCTScanParas.dNumV - 1) / blockSizeProj.y + 1);

	float angle = 0.0f;
	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // ���㵱ǰ��ת�ǣ���ת��Ϊ������

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
		//cudaMemcpy(tempIntX, d_mPolyForwardProj->phantom, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);
		//FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(tempIntX, 1, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), fp);
		//fclose(fp);

		computeSinMatIndensityProjKernel << <gridSizeProj, blockSizeProj >> > (d_mPolyForwardProj, mCTScanSystemInfo, mCTScanParas, i);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "forwardProjKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching forwardProjKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// ���Դ���
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// ���
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// ��ʱ
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>�ܶ�ͼ��ͶӰ��ʱ(GPU)��" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// դ����
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

}

void mallocDetResponse(PolyForwardProj& d_mPolyForwardProj, CTScanParas mCTScanParas)
{
	// ̽������Ӧ����
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
}



void freeDeviceMemory(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate)
{
	CUDAFREE(d_mCoordinate.imgIntX);
	CUDAFREE(d_mCoordinate.imgIntY);
	CUDAFREE(d_mCoordinate.imgIntZ);

	CUDAFREE(d_mPolyForwardProj.I);
	CUDAFREE(d_mPolyForwardProj.I0);
	CUDAFREE(d_mPolyForwardProj.IAbsorb);
	CUDAFREE(d_mPolyForwardProj.proj);
	CUDAFREE(d_mPolyForwardProj.phantom);
	CUDAFREE(d_mPolyForwardProj.phantomMassAtten);
	CUDAFREE(d_mPolyForwardProj.spectrumNormal);
	CUDAFREE(d_mPolyForwardProj.detResponse);
	CUDAFREE(d_mPolyForwardProj.grid);
	CUDAFREE(d_mPolyForwardProj.gridLinearAtten);
	CUDAFREE(d_mPolyForwardProj.scintillatorLineAtten);
	CUDAFREE(d_mPolyForwardProj.scintillatorPerThickness);
	CUDAFREE(d_mPolyForwardProj.foSpOffsetU);
	CUDAFREE(d_mPolyForwardProj.foSpOffsetV);

	cudaDestroyTextureObject(texObj);
	cudaFreeArray(d_cuArray3D);
}

// --------------------------------------------------------------------------------------------------

__global__ void computeIntPointCoordinatesKernel(Coordinate d_Coordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// ��������
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{

		// ̽������Ԫ����
		// ����ģ����ϵ��������任
		float detU = mCTScanSystemInfo.dHalfLU - mCTScanParas.dSize / 2 - y * mCTScanParas.dSize;
		float detV = mCTScanSystemInfo.dHalfLV - mCTScanParas.dSize / 2 - z * mCTScanParas.dSize;



		// ֻ�Ǵ˴����õ���ʱ��������˲���Ҫ�洢
		/*d_Coordinate.detU[y] = mCTScanSystemInfo.dHalfLU - mCTScanParas.dSize / 2 - y * mCTScanParas.dSize;
		d_Coordinate.detV[z] = mCTScanSystemInfo.dHalfLV - mCTScanParas.dSize / 2 - z * mCTScanParas.dSize;*/

		// Compute trigonometric value of Gamma angle and Beta angle
		// Gamma represents the angle between the ray and the xoy plane
		// Beta represents the angle between the ray and the xoz plane
		float sinGamma = detV / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));
		float cosGamma = sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));

		float sinBeta = detU / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));
		float cosBeta = sqrtf(powf(mCTScanParas.sdd, 2) + powf(detV, 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));


		/*float sinGamma = d_Coordinate.detV[z] / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detU[y], 2) + powf(d_Coordinate.detV[z], 2));
		float cosGamma = sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detU[y], 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detU[y], 2) + powf(d_Coordinate.detV[z], 2));

		float sinBeta = d_Coordinate.detU[y] / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detU[y], 2) + powf(d_Coordinate.detV[z], 2));
		float cosBeta = sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detV[z], 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate.detU[y], 2) + powf(d_Coordinate.detV[z], 2));*/


		/*float sinBeta = detY[y] / sqrtf(powf(sdd, 2) + powf(detY[y], 2) + powf(detZ[z], 2));
		float cosBeta = sqrtf(powf(sdd, 2) + powf(detZ[z], 2)) / sqrtf(powf(sdd, 2) + powf(detY[y], 2) + powf(detZ[z], 2));*/


		// Compute integration point coordinates
		size_t index = z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + x;

		d_Coordinate.imgIntX[index] = cosGamma * cosBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx) - mCTScanParas.sod;  // ��
		d_Coordinate.imgIntY[index] = cosGamma * sinBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // ��
		d_Coordinate.imgIntZ[index] = sinGamma * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // ҳ
	}
}


//__global__ void computeIntPointCoordinatesKernel(Coordinate* d_Coordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
//{
//	// ��������
//	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
//	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
//	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;
//
//	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
//	{
//
//		// ̽������Ԫ����
//		// ����ģ����ϵ��������任
//		/*float detU = mCTScanSystemInfo.dHalfLU / 2 - y * mCTScanParas.dSize;
//		float detV = mCTScanSystemInfo.dHalfLV / 2 - z * mCTScanParas.dSize;*/
//
//		d_Coordinate->detU[y] = mCTScanSystemInfo.dHalfLU / 2 - y * mCTScanParas.dSize;
//		d_Coordinate->detV[z] = mCTScanSystemInfo.dHalfLV / 2 - z * mCTScanParas.dSize;
//
//		// Compute trigonometric value of Gamma angle and Beta angle
//		// Gamma represents the angle between the ray and the xoy plane
//		// Beta represents the angle between the ray and the xoz plane
//		/*float sinGamma = detV / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));
//		float cosGamma = sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));
//
//		float sinBeta = detU / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));
//		float cosBeta = sqrtf(powf(mCTScanParas.sdd, 2) + powf(detV, 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(detU, 2) + powf(detV, 2));*/
//
//
//		float sinGamma = d_Coordinate->detU[y] / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detU[y], 2) + powf(d_Coordinate->detV[z], 2));
//		float cosGamma = sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detU[y], 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detU[y], 2) + powf(d_Coordinate->detV[z], 2));
//
//		float sinBeta = d_Coordinate->detU[y] / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detU[y], 2) + powf(d_Coordinate->detV[z], 2));
//		float cosBeta = sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detV[z], 2)) / sqrtf(powf(mCTScanParas.sdd, 2) + powf(d_Coordinate->detU[y], 2) + powf(d_Coordinate->detV[z], 2));
//
//
//		/*float sinBeta = detY[y] / sqrtf(powf(sdd, 2) + powf(detY[y], 2) + powf(detZ[z], 2));
//		float cosBeta = sqrtf(powf(sdd, 2) + powf(detZ[z], 2)) / sqrtf(powf(sdd, 2) + powf(detY[y], 2) + powf(detZ[z], 2));*/
//
//
//		// Compute integration point coordinates
//		size_t index = z * mCTScanParas.dNumU * mCTScanParas.dNumV + y * mCTScanParas.dNumU + x;
//
//		d_Coordinate->imgIntX[index] = cosGamma * cosBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx) - mCTScanParas.sod;  // ��
//		d_Coordinate->imgIntY[index] = cosGamma * sinBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // ��
//		d_Coordinate->imgIntZ[index] = sinGamma * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // ҳ
//	}
//}


__global__ void transformKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, float theta)
{
	// ������������
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		// �ǹ�һ��
		size_t index = z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + x;


		float tu = d_Coordinate.imgIntX[index] * cosf(theta) - d_Coordinate.imgIntY[index] * sinf(theta) + mCTScanSystemInfo.pHalfX;// mCTScanSystemInfo.FOVR;
		float tv = -(d_Coordinate.imgIntX[index] * sinf(theta) + d_Coordinate.imgIntY[index] * cosf(theta)) + mCTScanSystemInfo.pHalfY;//mCTScanSystemInfo.FOVR;     // ��ͶӰģ���У�ʹ������ϵ���˴��ı�y�᷽�򣬱������ʹ�õ�����ϵ
		float tw = -d_Coordinate.imgIntZ[index] + mCTScanSystemInfo.pHalfZ; //mCTScanSystemInfo.FOVH;   // �ı�Z�᷽�򣬷�������Ϊ������Ԥ�ڽ�ģһ��

		// �������ж�ȡ��д��ȫ�ִ洢
		d_mPolyForwardProj.phantom[index] = tex3D<float>(texObj, tu / mCTScanSystemInfo.pSizeX + 0.5, tv / mCTScanSystemInfo.pSizeY + 0.5, tw / mCTScanSystemInfo.pSizeZ + 0.5); //  �Ƿ���Ҫ��0.5����������   �������0.5Ҫ��Ϊ���ش�С
		// ����Y���ּ���������������ת�任ʱӦ��ȥ����ֵ����ͳһ�����ڲ�ֵ����ʱ��Ҫ���������ֵ�����Ե��������������ʡ�ԣ�X����Ϊ�˷���Ҳ��Y����ͬ��ʡ��
		// ���
	}
}

// ������ƫ��
__global__ void transformFocalSpotKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, int angleIndex)
{
	// ������������
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	float theta = mCTScanSystemInfo.rotatedDirection * angleIndex * mCTScanSystemInfo.thetaStep / 180.0f * PI;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		// �ǹ�һ��
		size_t index = z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + x;


		float tx = d_Coordinate.imgIntX[index] * cosf(theta) - d_Coordinate.imgIntY[index] * sinf(theta) + mCTScanSystemInfo.pHalfX;// mCTScanSystemInfo.FOVR;
		float ty = -(d_Coordinate.imgIntX[index] * sinf(theta) + d_Coordinate.imgIntY[index] * cosf(theta)) + mCTScanSystemInfo.pHalfY + d_mPolyForwardProj.foSpOffsetU[angleIndex];//mCTScanSystemInfo.FOVR;     // ��ͶӰģ���У�ʹ������ϵ���˴��ı�y�᷽�򣬱������ʹ�õ�����ϵ
		float tz = -d_Coordinate.imgIntZ[index] + mCTScanSystemInfo.pHalfZ + d_mPolyForwardProj.foSpOffsetV[angleIndex]; //mCTScanSystemInfo.FOVH;   // �ı�Z�᷽�򣬷�������Ϊ������Ԥ�ڽ�ģһ��
		// foSpOffset ����ƫ����


		// �������ж�ȡ��д��ȫ�ִ洢
		d_mPolyForwardProj.phantom[index] = tex3D<float>(texObj, tx / mCTScanSystemInfo.pSizeX + 0.5, ty / mCTScanSystemInfo.pSizeY + 0.5, tz / mCTScanSystemInfo.pSizeZ + 0.5);
		// ����Y���ּ���������������ת�任ʱӦ��ȥ����ֵ����ͳһ�����ڲ�ֵ����ʱ��Ҫ���������ֵ�����Ե��������������ʡ�ԣ�X����Ϊ�˷���Ҳ��Y����ͬ��ʡ��
		// ���
	}

}

__global__ void projOffsetMatchKernel(PolyForwardProj d_mPolyForwardProj, cudaTextureObject_t texObj, CTScanParas mCTScanParas)
{
	// ������������
	unsigned int u = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int v = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		// �ǹ�һ��

		size_t index = pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u;

		// foSpOffset ����ƫ����

		// �������ж�ȡ��д��ȫ�ִ洢
		d_mPolyForwardProj.phantom[index]
			= tex3D<float>(texObj
				, u + d_mPolyForwardProj.foSpOffsetU[pn] / mCTScanParas.dSize + 0.5
				, v + d_mPolyForwardProj.foSpOffsetV[pn] / mCTScanParas.dSize + 0.5
				, pn + 0.5);
	}
}

__global__ void forwardProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// ��X���������
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-temp * mCTScanSystemInfo.dx) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-temp * mCTScanSystemInfo.dx)) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void computeSinMatIndensityProjKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// ��X���������
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[v * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + u * mCTScanSystemInfo.intNum + i];
		}
		d_mPolyForwardProj.proj[num * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = temp * mCTScanSystemInfo.dx;    // �ܶȺ� �� ��������
	}
}

__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// ��X���������
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx)) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];
	}
}

__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// ��X���������
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.grid[u] * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.grid[u] * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];
	}
}

__global__ void forwardSinMatNoResponseProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// ��X���������
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.grid[u];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.grid[u];
	}
}

__global__ void forwardProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// ��X���������
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* expf(-temp * mCTScanSystemInfo.dx)
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* (1 - expf(-temp * mCTScanSystemInfo.dx))
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// ��X���������
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx)
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];         // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* (1 - expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx))
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// ��X���������
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];

	}

}

__global__ void forwardSinMatNoResponseProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// ��X���������
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]);  // ��X�������, ÿһ���Ƕ���Ϊһ��, Y����Ϊ�У� Z����Ϊ��
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]));
	}
}



// ��ʼ����ά����, �������Դ����, ��������ʱֻ�����cuArray, ���ݴ����� Host to Device.
// texObj -- �������, d_cuArray3D -- Device�д洢���ݵ�ָ��, data -- Դ����, volumeSize -- �����С.
void initTexture3D(cudaTextureObject_t& texObj, cudaArray_t& d_cuArray3D, float* h_data, cudaExtent volumeSize)
{
	//cudaExtent volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	//cudaArray* d_cuArray3D;
	cudaMalloc3DArray(&d_cuArray3D, &channelDesc, volumeSize);

	// �������ݵ�CUDA array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)h_data, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_cuArray3D;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);

	// ������Դ������
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(cudaResourceDesc));

	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_cuArray3D;

	// ��������������
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 0;

	// �����������
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);
}

// ��������
void updateTex(cudaArray_t& d_cuArray3D, float* h_data, cudaExtent volumeSize)
{
	// �������ݵ�CUDA array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)h_data, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_cuArray3D;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);
}

void createTexture3D(cudaTextureObject_t& texObj, float* h_data, size_t width, size_t height, size_t depth)
{
	cudaExtent volumeSize = make_cudaExtent(width, height, depth);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* d_cuArray3D;
	cudaMalloc3DArray(&d_cuArray3D, &channelDesc, volumeSize);

	// �������ݵ�CUDA array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)h_data, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_cuArray3D;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);

	// ������Դ������
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(cudaResourceDesc));

	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_cuArray3D;

	// ��������������
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 0;

	// �����������
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);
}




