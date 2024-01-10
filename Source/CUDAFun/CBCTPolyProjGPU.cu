

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
	// 模体纹理更新
	updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// 栅数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);
	// 探测器响应数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);

	//// 调试代码
	//FILE* fp;
	//fp = fopen("test.raw", "wb");
	//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
	//fclose(fp);
}


//// 带栅单材质正投影
//void forwardSinMatProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
//{
//	//// 模体纹理更新
//	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);  // 更新密度
//
//	// 栅数据 主存 -> 显存
//	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);
//	
//	// 探测器响应数据 主存 -> 显存
//	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);
//
//
//	// GPU计时
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
//		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制
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
//		//// 调试代码
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
//		//// 调试代码
//		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//		//// 检测
//		////FILE* fp;
//		//fp = fopen("test.raw", "wb");
//		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
//		//fclose(fp);
//	}
//
//	// I  显存 -> 主存
//	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
//
//
//	// 计时
//	cudaEventRecord(g_stop, 0);
//	cudaEventSynchronize(g_stop);
//	float elapsedTime = 0;
//	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
//	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
//	cudaEventDestroy(g_start);
//	cudaEventDestroy(g_stop);
//
//}

// 带栅单材质正投影
void forwardSinMatProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	//// 模体纹理更新
	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);  // 更新密度

	// 栅数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// 探测器响应数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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

	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);

}

void forwardSinMatNoResponseProjGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// 栅数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// 探测器响应数据 主存 -> 显存
	// cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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

	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}


// 无栅正投影
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
	// 模体纹理更新
	updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// 探测器响应数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}


// 无栅单材质正投影
void forwardSinMatProjNoGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{


	//// 模体纹理更新
	//updateTex(d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);

	// 探测器响应数据 主存 -> 显存
	cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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


	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);
}

void forwardSinMatNoResponseProjNoGridGPU(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// 栅数据 主存 -> 显存
	//cudaMemcpy(d_mPolyForwardProj.grid, h_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float), cudaMemcpyHostToDevice);

	// 探测器响应数据 主存 -> 显存
	// cudaMemcpy(d_mPolyForwardProj.detResponse, h_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyHostToDevice);


	// GPU计时
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

	// I  显存 -> 主存
	cudaMemcpy(h_mPolyForwardProj.I, d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_mPolyForwardProj.IAbsorb, d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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


	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}
	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// 探测器响应数据
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

	// 调试代码
	//h_mPolyForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	//cudaMemcpy(h_mPolyForwardProj.proj, d_mPolyForwardProj.proj, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
	//// 检测
	//FILE* fp;
	//fp = fopen("test.raw", "wb");
	//fwrite(h_mPolyForwardProj.proj, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
	//fclose(fp);

	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// 探测器响应数据
	cudaMalloc(&d_mPolyForwardProj.detResponse, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
}

// 初始化，单材质、有栅、有响应、有焦斑
void initDeviceSinMatFoSpSiGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	cudaError_t cudaStatus;

	// Select Gpu
	cudaStatus = cudaSetDevice(GPUINDEX);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed? \n Error: %s\n", cudaGetErrorString(cudaStatus));
	}

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 计算积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 焦点偏移数据发送
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

	// 传输数据
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


	// GPU计时
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

		//angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformFocalSpotKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, i);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformFocalSpotKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformFocalSpotKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// 最终投影进行偏移，做到焦点偏移，相当于物体和探测器同时朝同一方向偏移
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

	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// 探测器响应数据
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

}

void initDeviceNoGrid(PolyForwardProj& d_mPolyForwardProj, Coordinate& d_mCoordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, PolyForwardProj& h_mPolyForwardProj)
{
	// 初始化device
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


	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Variable d_mPolyForwardProj->phantom cudaMalloc failed!");
	}
	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 探测器响应数据
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// 探测器响应数据
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 焦点偏移数据发送
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

	// 传输数据
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

	// GPU计时
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

		//angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformFocalSpotKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, i);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformFocalSpotKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformFocalSpotKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// 最终投影进行偏移，做到焦点偏移，相当于物体和探测器同时朝同一方向偏移
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


	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

	// 探测器响应数据
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

	// 初始化device
	cudaDeviceSynchronize();    // 确保无CUDA程序运行
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		printf("Error: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaDeviceReset();  // 清理缓存


	// Allocate GPU buffers
	size_t sizeIntPoint = mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float);

	// 积分点坐标
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

	// 检查
	/*float* tempIntX = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntX, d_mCoordinate->imgIntX, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntY = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntY, d_mCoordinate->imgIntY, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);

	float* tempIntZ = new float[mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV];
	cudaMemcpy(tempIntZ, d_mCoordinate->imgIntZ, mCTScanSystemInfo.intNum * mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float), cudaMemcpyDeviceToHost);*/


	// 初始化纹理
	volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	initTexture3D(texObj, d_cuArray3D, h_mPolyForwardProj.phantom, volumeSize);


	// Allocate GPU buffers for temporary output variables
	cudaStatus = cudaMalloc(&d_mPolyForwardProj.phantom, sizeIntPoint);   // 存储积分点处像素值旋转后的数据
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.phantom cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// 计算 密度图像的Proj
	// 密度图像的Proj
	cudaMalloc(&d_mPolyForwardProj.proj, sizeIntPoint);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "Variable d_mPolyForwardProj.proj cudaMalloc failed! \n Error Code: %d --- %s! \n\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}

	// GPU计时
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

		angle = mCTScanSystemInfo.rotatedDirection * i * mCTScanSystemInfo.thetaStep / 180 * PI;   // 计算当前旋转角，并转换为弧度制

		transformKernel << <gridSizeTrans, blockSizeTrans >> > (d_mPolyForwardProj, d_mCoordinate, texObj, mCTScanSystemInfo, mCTScanParas, angle);

		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "transformKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		}cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching transformKernel!\nError: %s!\n", cudaStatus, cudaGetErrorString(cudaStatus));
		}

		//// 调试代码
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

		//// 调试代码
		//cudaMemcpy(h_mPolyForwardProj->I, d_mPolyForwardProj->I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), cudaMemcpyDeviceToHost);
		//// 检测
		////FILE* fp;
		//fp = fopen("test.raw", "wb");
		//fwrite(h_mPolyForwardProj->I, 1, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float), fp);
		//fclose(fp);
	}

	// 计时
	cudaEventRecord(g_stop, 0);
	cudaEventSynchronize(g_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime, g_start, g_stop);
	std::cout << "==>>密度图像投影耗时(GPU)：" << elapsedTime / 1000.0f << " s" << std::endl;
	cudaEventDestroy(g_start);
	cudaEventDestroy(g_stop);



	// Allocate GPU buffers and host
	cudaMalloc(&d_mPolyForwardProj.I, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	cudaMalloc(&d_mPolyForwardProj.IAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));

	// 栅数据
	//cudaMalloc(&d_mPolyForwardProj.grid, mCTScanParas.dNumU * sizeof(float));

}

void mallocDetResponse(PolyForwardProj& d_mPolyForwardProj, CTScanParas mCTScanParas)
{
	// 探测器响应数据
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
	// 计算坐标
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{

		// 探测器单元坐标
		// 依建模坐标系进行坐标变换
		float detU = mCTScanSystemInfo.dHalfLU - mCTScanParas.dSize / 2 - y * mCTScanParas.dSize;
		float detV = mCTScanSystemInfo.dHalfLV - mCTScanParas.dSize / 2 - z * mCTScanParas.dSize;



		// 只是此处所用的临时变量，因此不需要存储
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

		d_Coordinate.imgIntX[index] = cosGamma * cosBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx) - mCTScanParas.sod;  // 列
		d_Coordinate.imgIntY[index] = cosGamma * sinBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // 行
		d_Coordinate.imgIntZ[index] = sinGamma * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // 页
	}
}


//__global__ void computeIntPointCoordinatesKernel(Coordinate* d_Coordinate, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
//{
//	// 计算坐标
//	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
//	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
//	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;
//
//	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
//	{
//
//		// 探测器单元坐标
//		// 依建模坐标系进行坐标变换
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
//		d_Coordinate->imgIntX[index] = cosGamma * cosBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx) - mCTScanParas.sod;  // 列
//		d_Coordinate->imgIntY[index] = cosGamma * sinBeta * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // 行
//		d_Coordinate->imgIntZ[index] = sinGamma * (mCTScanParas.sod - mCTScanSystemInfo.FOVR + mCTScanSystemInfo.dx / 2 + x * mCTScanSystemInfo.dx);  // 页
//	}
//}


__global__ void transformKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, float theta)
{
	// 计算纹理坐标
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		// 非归一化
		size_t index = z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + x;


		float tu = d_Coordinate.imgIntX[index] * cosf(theta) - d_Coordinate.imgIntY[index] * sinf(theta) + mCTScanSystemInfo.pHalfX;// mCTScanSystemInfo.FOVR;
		float tv = -(d_Coordinate.imgIntX[index] * sinf(theta) + d_Coordinate.imgIntY[index] * cosf(theta)) + mCTScanSystemInfo.pHalfY;//mCTScanSystemInfo.FOVR;     // 正投影模型中，使用右手系，此处改变y轴方向，变回纹理使用的坐标系
		float tw = -d_Coordinate.imgIntZ[index] + mCTScanSystemInfo.pHalfZ; //mCTScanSystemInfo.FOVH;   // 改变Z轴方向，符合向上为正，与预期建模一致

		// 从纹理中读取并写入全局存储
		d_mPolyForwardProj.phantom[index] = tex3D<float>(texObj, tu / mCTScanSystemInfo.pSizeX + 0.5, tv / mCTScanSystemInfo.pSizeY + 0.5, tw / mCTScanSystemInfo.pSizeZ + 0.5); //  是否需要加0.5？！！？？   测试完后0.5要改为像素大小
		// 由于Y发现间隔不定，因此在旋转变换时应减去的数值不能统一，而在插值索引时需要加上这个数值，可以抵消，因此两处都省略，X方向为了方便也和Y方向同样省略
		// 检测
	}
}

// 含焦点偏移
__global__ void transformFocalSpotKernel(PolyForwardProj d_mPolyForwardProj, Coordinate d_Coordinate, cudaTextureObject_t texObj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, int angleIndex)
{
	// 计算纹理坐标
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	float theta = mCTScanSystemInfo.rotatedDirection * angleIndex * mCTScanSystemInfo.thetaStep / 180.0f * PI;

	if (x < mCTScanSystemInfo.intNum && y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		// 非归一化
		size_t index = z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + x;


		float tx = d_Coordinate.imgIntX[index] * cosf(theta) - d_Coordinate.imgIntY[index] * sinf(theta) + mCTScanSystemInfo.pHalfX;// mCTScanSystemInfo.FOVR;
		float ty = -(d_Coordinate.imgIntX[index] * sinf(theta) + d_Coordinate.imgIntY[index] * cosf(theta)) + mCTScanSystemInfo.pHalfY + d_mPolyForwardProj.foSpOffsetU[angleIndex];//mCTScanSystemInfo.FOVR;     // 正投影模型中，使用右手系，此处改变y轴方向，变回纹理使用的坐标系
		float tz = -d_Coordinate.imgIntZ[index] + mCTScanSystemInfo.pHalfZ + d_mPolyForwardProj.foSpOffsetV[angleIndex]; //mCTScanSystemInfo.FOVH;   // 改变Z轴方向，符合向上为正，与预期建模一致
		// foSpOffset 焦点偏移量


		// 从纹理中读取并写入全局存储
		d_mPolyForwardProj.phantom[index] = tex3D<float>(texObj, tx / mCTScanSystemInfo.pSizeX + 0.5, ty / mCTScanSystemInfo.pSizeY + 0.5, tz / mCTScanSystemInfo.pSizeZ + 0.5);
		// 由于Y发现间隔不定，因此在旋转变换时应减去的数值不能统一，而在插值索引时需要加上这个数值，可以抵消，因此两处都省略，X方向为了方便也和Y方向同样省略
		// 检测
	}

}

__global__ void projOffsetMatchKernel(PolyForwardProj d_mPolyForwardProj, cudaTextureObject_t texObj, CTScanParas mCTScanParas)
{
	// 计算纹理坐标
	unsigned int u = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int v = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		// 非归一化

		size_t index = pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u;

		// foSpOffset 焦点偏移量

		// 从纹理中读取并写入全局存储
		d_mPolyForwardProj.phantom[index]
			= tex3D<float>(texObj
				, u + d_mPolyForwardProj.foSpOffsetU[pn] / mCTScanParas.dSize + 0.5
				, v + d_mPolyForwardProj.foSpOffsetV[pn] / mCTScanParas.dSize + 0.5
				, pn + 0.5);
	}
}

__global__ void forwardProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// 沿X正方向积分
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-temp * mCTScanSystemInfo.dx) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-temp * mCTScanSystemInfo.dx)) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void computeSinMatIndensityProjKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// 沿X正方向积分
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[v * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + u * mCTScanSystemInfo.intNum + i];
		}
		d_mPolyForwardProj.proj[num * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = temp * mCTScanSystemInfo.dx;    // 密度和 乘 积分区间
	}
}

__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// 沿X正方向积分
	size_t y = blockIdx.x * blockDim.x + threadIdx.x;
	size_t z = blockIdx.y * blockDim.y + threadIdx.y;

	if (y < mCTScanParas.dNumU && z < mCTScanParas.dNumV)
	{
		float temp = 0.0f;
		for (size_t i = 0; i < mCTScanSystemInfo.intNum; i++)
		{
			temp += d_mPolyForwardProj.phantom[z * mCTScanSystemInfo.intNum * mCTScanParas.dNumU + y * mCTScanSystemInfo.intNum + i];
		}

		d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx)) * d_mPolyForwardProj.grid[y] * d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];
	}
}

__global__ void forwardSinMatProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// 沿X正方向积分
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.grid[u] * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.grid[u] * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];
	}
}

__global__ void forwardSinMatNoResponseProjGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// 沿X正方向积分
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.grid[u];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.grid[u];
	}
}

__global__ void forwardProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// 沿X正方向积分
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
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* (1 - expf(-temp * mCTScanSystemInfo.dx))
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas, size_t num)
{
	// 沿X正方向积分
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
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];         // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y]
			= mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal
			* (1 - expf(-mCTScanSystemInfo.phantomMAtten * temp * mCTScanSystemInfo.dx))
			* d_mPolyForwardProj.detResponse[z * mCTScanParas.dNumU + y];

		// 
		//temp = mCTScanParas.I0Val * expf(-temp * mCTScanSystemInfo.dx);
		//d_mPolyForwardProj.I[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * temp * d_mPolyForwardProj.grid[y];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		//d_mPolyForwardProj.IAbsorb[num * mCTScanParas.dNumU * mCTScanParas.dNumV + z * mCTScanParas.dNumU + y] = mCTScanSystemInfo.spectrumVal * (mCTScanParas.I0Val - temp) * d_mPolyForwardProj.grid[y];

	}
}

__global__ void forwardSinMatProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// 沿X正方向积分
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]) * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u])) * d_mPolyForwardProj.detResponse[v * mCTScanParas.dNumU + u];

	}

}

__global__ void forwardSinMatNoResponseProjNoGridKernel(PolyForwardProj d_mPolyForwardProj, CTScanSystemInfo mCTScanSystemInfo, CTScanParas mCTScanParas)
{
	// 沿X正方向积分
	size_t u = blockIdx.x * blockDim.x + threadIdx.x;
	size_t v = blockIdx.y * blockDim.y + threadIdx.y;
	size_t pn = blockIdx.z * blockDim.z + threadIdx.z;

	if (u < mCTScanParas.dNumU && v < mCTScanParas.dNumV && pn < mCTScanParas.projNum)
	{
		d_mPolyForwardProj.I[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]);  // 沿X方向积分, 每一个角度作为一层, Y方向为列， Z方向为行
		d_mPolyForwardProj.IAbsorb[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u] = mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * (1 - expf(-mCTScanSystemInfo.phantomMAtten * d_mPolyForwardProj.proj[pn * mCTScanParas.dNumU * mCTScanParas.dNumV + v * mCTScanParas.dNumU + u]));
	}
}



// 初始化三维纹理, 绑定纹理和源数据, 更新纹理时只需更新cuArray, 数据传输是 Host to Device.
// texObj -- 纹理对象, d_cuArray3D -- Device中存储数据的指针, data -- 源数据, volumeSize -- 纹理大小.
void initTexture3D(cudaTextureObject_t& texObj, cudaArray_t& d_cuArray3D, float* h_data, cudaExtent volumeSize)
{
	//cudaExtent volumeSize = make_cudaExtent(mCTScanParas.pNumX, mCTScanParas.pNumY, mCTScanParas.pNumZ);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	//cudaArray* d_cuArray3D;
	cudaMalloc3DArray(&d_cuArray3D, &channelDesc, volumeSize);

	// 拷贝数据到CUDA array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)h_data, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_cuArray3D;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);

	// 定义资源描述符
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(cudaResourceDesc));

	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_cuArray3D;

	// 定义纹理对象参数
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 0;

	// 生产纹理对象
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);
}

// 更新纹理
void updateTex(cudaArray_t& d_cuArray3D, float* h_data, cudaExtent volumeSize)
{
	// 拷贝数据到CUDA array
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

	// 拷贝数据到CUDA array
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr((void*)h_data, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
	copyParams.dstArray = d_cuArray3D;
	copyParams.extent = volumeSize;
	copyParams.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams);

	// 定义资源描述符
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(cudaResourceDesc));

	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_cuArray3D;

	// 定义纹理对象参数
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModeLinear;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 0;

	// 创建纹理对象
	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);
}




