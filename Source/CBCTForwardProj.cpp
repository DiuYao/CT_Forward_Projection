
#include <iostream>
#include <fstream>
#include <cmath>
#include <direct.h>
#include <corecrt_io.h>

#include "CBCTForwardProj.h"



using namespace std;

CBCTForwardProj::CBCTForwardProj()
{
}

CBCTForwardProj::~CBCTForwardProj()
{

}

void CBCTForwardProj::computeParas()
{
	mCTScanSystemInfo.rotatedDirection = mCTScanParas.rotatedDirection;

	// Compute parameters
	mCTScanSystemInfo.thetaStep = (float)360 / mCTScanParas.projNum;
	mCTScanSystemInfo.dHalfLU = mCTScanParas.dNumU * mCTScanParas.dSize / 2;
	mCTScanSystemInfo.dHalfLV = mCTScanParas.dNumV * mCTScanParas.dSize / 2;

	// FOV
	mCTScanSystemInfo.FOVR = mCTScanParas.sod * mCTScanSystemInfo.dHalfLU
		/ sqrtf(mCTScanSystemInfo.dHalfLU * mCTScanSystemInfo.dHalfLU + mCTScanParas.sdd * mCTScanParas.sdd);
	mCTScanSystemInfo.FOVH = (mCTScanParas.sod - mCTScanSystemInfo.FOVR)
		* mCTScanSystemInfo.dHalfLV / mCTScanParas.sdd;   // 圆柱高
	// mCTScanSystemInfo.FOVH = mCTScanParas.sod * mCTScanSystemInfo.dHalfLV / sqrt(mCTScanSystemInfo.dHalfLV * mCTScanSystemInfo.dHalfLV + mCTScanParas.sdd * mCTScanParas.sdd);                  // 垂直视野圆半径

	// Pixel size
	
	/*mCTScanSystemInfo.pSizeX = 2 * mCTScanSystemInfo.FOVR / mCTScanParas.pNumX;
	mCTScanSystemInfo.pSizeY = 2 * mCTScanSystemInfo.FOVR / mCTScanParas.pNumY;
	mCTScanSystemInfo.pSizeZ = 2 * mCTScanSystemInfo.FOVH / mCTScanParas.pNumZ;*/

	mCTScanSystemInfo.pSizeX = mCTScanParas.pLengthX / mCTScanParas.pNumX;
	mCTScanSystemInfo.pSizeY = mCTScanParas.pLengthY / mCTScanParas.pNumY;
	mCTScanSystemInfo.pSizeZ = mCTScanParas.pLengthZ / mCTScanParas.pNumZ;

	// Phantom Half Length
	mCTScanSystemInfo.pHalfX = mCTScanParas.pLengthX / 2;
	mCTScanSystemInfo.pHalfY = mCTScanParas.pLengthY / 2;
	mCTScanSystemInfo.pHalfZ = mCTScanParas.pLengthZ / 2;

	// Integral calculus
	mCTScanSystemInfo.dx = 0.5 * mCTScanSystemInfo.pSizeX;
	mCTScanSystemInfo.dy = 0.5 * mCTScanSystemInfo.pSizeY;
	mCTScanSystemInfo.dz = 0.5 * mCTScanSystemInfo.pSizeZ;

	// Compute detector coordinates
	mCTScanSystemInfo.intNum = roundf(2 * mCTScanSystemInfo.FOVR / mCTScanSystemInfo.dx);        // 某条射线上的积分点个数

}

// 不同能量下的剩余光子和物体吸收光子分别累加
void CBCTForwardProj::addIandIAbsorb()
{
	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum; i++)
	{
		I[i] += h_mForwardProj.I[i];
	}

	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum; i++)
	{
		IAbsorb[i] += h_mForwardProj.IAbsorb[i];
	}
}

// 剩余光子和物体吸收光子数据类型转换 float —> uint32
void CBCTForwardProj::convertIandIAbsorbDt()
{
	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum; i++)
	{
		IPolyenergetic[i] = roundf(I[i]);
	}

	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum; i++)
	{
		IPolyenergeticAbsorb[i] = roundf(IAbsorb[i]);
	}
}

// 生成闪烁体厚度，模拟工艺导致厚度不一致
void CBCTForwardProj::computePerScinThinckness()
{
	srand(100);  // 随机数种子
	float tempErr = 0.0f;
	int beg = -mCTScanParas.mScintilltorInfo.scintillatorThicknessErr * 10000;  // 放大为整数
	int endd = mCTScanParas.mScintilltorInfo.scintillatorThicknessErr * 10000;

	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV; i++)
	{
		tempErr = -mCTScanParas.mScintilltorInfo.scintillatorThicknessErr + (rand() % (endd - beg + 1)) / 10000.0f;  // 产生随机误差 
		h_mForwardProj.scintillatorPerThickness[i] = mCTScanParas.mScintilltorInfo.scintillatorThickness + tempErr;
	}
}

// 更新探测器响应，随能量变化
void CBCTForwardProj::updateDetResponse()
{
	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV; i++)
	{
		h_mForwardProj.detResponse[i] = mCTScanParas.mScintilltorInfo.detResponseFactor * (1 - expf(-h_mForwardProj.scintillatorLineAtten[specIndex] * h_mForwardProj.scintillatorPerThickness[i]));
	}
}

void CBCTForwardProj::readSpecrtumNorm()
{
	ifstream ifs;
	ifs.open(mFilePath.spectrumPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "文件打开失败" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		h_mForwardProj.spectrumNormal[ii] = buf;
		ii++;
	}

	ifs.close();
}

// 读取滤线栅的质量衰减系数(cm^2/g)，转换为线性衰减系数(1/mm)
void CBCTForwardProj::readGridMassAttu()
{
	ifstream ifs;
	ifs.open(mFilePath.gridPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "栅衰减系数文件打开失败" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		// 读取的是质量衰减系数，需要转换为线性衰减系数，并将单位转化为 1/mm
		h_mForwardProj.gridLineAtten[ii] = buf;
		// 转换
		h_mForwardProj.gridLineAtten[ii] = h_mForwardProj.gridLineAtten[ii] * mGridInfo.rhoGS / 10;    // 1/mm
		ii++;
	}

	ifs.close();
}

// 读取闪烁体的质量衰减系数(cm^2/g)，转换为线性衰减系数(1/mm)
void CBCTForwardProj::readScintillatorMassAttu()
{
	ifstream ifs;
	ifs.open(mFilePath.scintillatorPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "闪烁体质量衰减系数文件打开失败！" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		// 读取的是质量衰减系数，需要转换为线性衰减系数，并将单位转化为 1/mm
		h_mForwardProj.scintillatorLineAtten[ii] = buf;
		// 转换
		h_mForwardProj.scintillatorLineAtten[ii] = h_mForwardProj.scintillatorLineAtten[ii] * mCTScanParas.mScintilltorInfo.scintillatorDensity / 10;  // 1/mm
		ii++;
	}

	ifs.close();
}

// 获得中间位置探元的响应
void CBCTForwardProj::getADetResponse()
{
	aDetResponse[specIndex] = h_mForwardProj.detResponse[lroundl((mCTScanParas.dNumU + mCTScanParas.dNumV) / 2)];
}


void CBCTForwardProj::saveADetResponse()
{
	string directory = "OutputResult"; // 指定目录名
	string filename = "ADetResponse.txt"; // 文件名
	
	if (_access(directory.c_str(), 0) == -1)
	{
		// if this folder not exist, create a new one.
		// 创建目录
		if (_mkdir(directory.c_str()) == -1)
		{
			cerr << "Error: Creating ADetResponse file is fair！" << endl;
		}
		// 返回 0 表示创建成功，-1 表示失败
	}
	


	// 打开文件
	ofstream ofs;
	ofs.open(directory + "/" + filename, ios::out);

	if (!ofs.is_open()) {
		cerr << "Error: Cannot open ADetResponse file" << endl;
		return ;
	}

	// 写入数据
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		ofs << aDetResponse[i] << endl;
	}

	// 关闭文件
	ofs.close();
}

void CBCTForwardProj::showScanSystemInfo()
{
	cout.setf(ios::left);
	cout << "*******************成像系统信息******************" << endl;
	cout.width(15);
	cout << "探测器尺寸:" << mCTScanParas.dNumU << "×" << mCTScanParas.dNumV << endl;
	cout.width(15);
	cout << "探测器分辨率:"; cout.width(8); cout << mCTScanParas.dSize << " mm" << endl;
	cout.width(15);
	cout << "SDD:"; cout.width(8); cout << mCTScanParas.sdd << " mm" << endl;
	cout.width(15);
	cout << "SOD:"; cout.width(8); cout << mCTScanParas.sod << " mm" << endl;
	cout.width(15);
	cout << "FOV(直径×长):"; cout.width(8); cout << 2 * mCTScanSystemInfo.FOVR << "×" << 2 * mCTScanSystemInfo.FOVH << " mm" << endl;
	cout.width(15);
	cout << "模体像素大小:"; cout << mCTScanSystemInfo.pSizeX << "×" << mCTScanSystemInfo.pSizeY << "×" << mCTScanSystemInfo.pSizeZ << " mm" << endl;
	cout << endl;

	/*cout.setf(ios::left);
	cout << "*******************成像系统信息******************" << endl;
	cout.width(15);
	cout << "探测器尺寸:" << dNums << "x" << dNums << endl;
	cout.width(15);
	cout << "探测器分辨率:"; cout.width(8); cout << dSize << " um" << endl;

	cout << endl;

	cout << "*****************聚焦型滤线栅信息****************" << endl;
	cout.width(15);
	cout << "栅条材质:" << materialGridStrip << endl;
	cout.width(15);
	cout << "线性衰减系数:"; cout.width(10); cout << uGridStrip * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "间隙物材质:" << materialGridInterspacer << endl;
	cout.width(15);
	cout << "线性衰减系数:"; cout.width(10); cout << uInterspacer * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "周期:"; cout.width(10); cout << gridPeriod << " um" << endl;
	cout.width(15);
	cout << "栅条宽度:"; cout.width(10); cout << leadStripsWidth << " um" << endl;
	cout.width(15);
	cout << "间隙距离:"; cout.width(10); cout << leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "厚度:"; cout.width(10); cout << h / 1000 << " mm" << endl;
	cout.width(15);
	cout << "焦距:"; cout.width(10); cout << FD / 1000 << " mm" << endl;
	cout << "*************************************************" << endl;*/
}

void CBCTForwardProj::showGridmInfo()
{
	cout << endl;
	cout << "*****************聚焦型滤线栅信息****************" << endl;
	cout.width(15);
	cout << "栅条材质:" << mGridInfo.materialGridStrip << endl;
	//cout.width(15);
	//cout << "线性衰减系数:"; cout.width(10); cout << uGridStrip * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "间隙物材质:" << mGridInfo.materialGridInterspacer << endl;
	//cout.width(15);
	//cout << "线性衰减系数:"; cout.width(10); cout << uInterspacer * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "周期:"; cout.width(10); cout << mGridInfo.leadStripsWidth + mGridInfo.leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "栅条宽度:"; cout.width(10); cout << mGridInfo.leadStripsWidth << " um" << endl;
	cout.width(15);
	cout << "间隙距离:"; cout.width(10); cout << mGridInfo.leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "厚度:"; cout.width(10); cout << mGridInfo.h << " mm" << endl;
	cout.width(15);
	cout << "焦距:"; cout.width(10); cout << mGridInfo.FD << " mm" << endl;
	cout << "*************************************************" << endl;

}

void CBCTForwardProj::computeProj()
{
	for (size_t j = 0; j < mCTScanParas.projNum; j++)
	{
		for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV; i++)
		{
			h_mForwardProj.proj[j * mCTScanParas.dNumU * mCTScanParas.dNumV + i]
				= logf(h_mForwardProj.I0[i] / IPolyenergetic[j * mCTScanParas.dNumU * mCTScanParas.dNumV + i]);
		}
	}
	
}
