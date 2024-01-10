
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "CBCTForwardProj.h"


#define PI acosf(-1)

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
		IPolyenergetic[i] = floorf(I[i]);
	}

	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum; i++)
	{
		IPolyenergeticAbsorb[i] = floorf(IAbsorb[i]);
	}
}

void CBCTForwardProj::addStandardI0(float& I0Standard)
{
	I0Standard += mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * mCTScanParas.mScintilltorInfo.detResponseFactor * (1 - expf(-h_mForwardProj.scintillatorLineAtten[specIndex] * mCTScanParas.mScintilltorInfo.scintillatorThickness));
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

void CBCTForwardProj::computeFoSpOffset()
{
	//srand(time(nullptr));

	srand(10);  // 便于实验结果对比，固定种子，实际中应该不固定

	float foSpRadius = mCTScanParas.focalSpotSize / 2.0f;
	int tempFoSpSize = mCTScanParas.focalSpotSize * 10000;

	float offsetU = 0.0f;
	float offsetV = 0.0f;

	for (size_t i = 0; i < mCTScanParas.projNum; i++)
	{
		offsetU = -foSpRadius + (rand() % tempFoSpSize) / 10000.0f;
		offsetV = -foSpRadius + (rand() % tempFoSpSize) / 10000.0f;
		while (offsetU * offsetU + offsetV * offsetV > PI * foSpRadius * foSpRadius)
		{
			offsetU = -foSpRadius + (rand() % tempFoSpSize) / 10000.0f;
			offsetV = -foSpRadius + (rand() % tempFoSpSize) / 10000.0f;
		}

		h_mForwardProj.foSpOffsetU[i] = offsetU;
		h_mForwardProj.foSpOffsetV[i] = offsetV;
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
	while (ifs >> buf && ii < mCTScanParas.specEnergyNum)
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

	// 读取栅条质量衰减系数
	ifs.open(mGridInfo.gridStripMAPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "栅条质量衰减系数文件打开失败！" << endl;
		system("pause");
		exit(0);
	}
	else
	{
		int ii = 0;
		float buf;
		while (ifs >> buf && ii < mCTScanParas.specEnergyNum)
		{
			// 读取的是质量衰减系数，需要转换为线性衰减系数，并将单位转化为 1/mm
			//h_mForwardProj.gridLinearAtten[ii] = buf;
			// 转换
			h_mForwardProj.gridLinearAtten[ii] = buf * mGridInfo.gridStripDensity / 10;    // 1/mm
			ii++;
		}
		ifs.close();
	}

	// 读取间隔填充物质量衰减系数
	ifs.open(mGridInfo.interspaceMAPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "间隔填充物质量衰减系数文件打开失败！" << endl;
		system("pause");
		exit(0);
	}
	else
	{
		int ii = 0;
		float buf;
		while (ifs >> buf && ii < mCTScanParas.specEnergyNum)
		{
			// 读取的是质量衰减系数，需要转换为线性衰减系数，并将单位转化为 1/mm
			//h_mForwardProj.gridLinearAtten[ii] = buf;
			// 转换
			h_mForwardProj.interspaceLinearAtten[ii] = buf * mGridInfo.interspaceDensity / 10;    // 1/mm
			ii++;
		}
		ifs.close();
	}
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
	while (ifs >> buf && ii < mCTScanParas.specEnergyNum)
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

	string tempPath = mFilePath.outputFolder + "/ADetResponse.txt";

	if (_access(mFilePath.outputFolder.c_str(), 0) == -1)
	{
		// if this folder not exist, create a new one.
		// 创建目录
		if (_mkdir(mFilePath.outputFolder.c_str()) == -1)
		{
			cerr << "Error: Creating ADetResponse file is fair！" << endl;
		}
		// 返回 0 表示创建成功，-1 表示失败
	}


	// 打开文件
	ofstream ofs;
	ofs.open(tempPath, ios::out);

	if (!ofs.is_open())
	{
		cerr << "Error: Cannot open ADetResponse file" << endl;
		return;
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
	cout.width(15);
	cout << "程序中初始光子数(非I0):"; cout << mCTScanParas.I0Val << endl;
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
	cout << "栅条材质:" << gridStripMaterial << endl;
	cout.width(15);
	cout << "线性衰减系数:"; cout.width(10); cout << gridStripLinearAtten * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "间隙物材质:" << interspaceMaterial << endl;
	cout.width(15);
	cout << "线性衰减系数:"; cout.width(10); cout << interspaceLinearAtten * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "周期:"; cout.width(10); cout << gridPeriod << " um" << endl;
	cout.width(15);
	cout << "栅条宽度:"; cout.width(10); cout << gridStripsWidth << " um" << endl;
	cout.width(15);
	cout << "间隙距离:"; cout.width(10); cout << interspcaceWidth << " um" << endl;
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
	cout << "栅条材质:" << mGridInfo.gridStripMaterial << endl;
	//cout.width(15);
	//cout << "线性衰减系数:"; cout.width(10); cout << gridStripLinearAtten * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "间隙物材质:" << mGridInfo.interspaceMaterial << endl;
	//cout.width(15);
	//cout << "线性衰减系数:"; cout.width(10); cout << interspaceLinearAtten * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "周期:"; cout.width(10); cout << mGridInfo.gridStripsWidth + mGridInfo.interspcaceWidth << " um" << endl;
	cout.width(15);
	cout << "栅条宽度:"; cout.width(10); cout << mGridInfo.gridStripsWidth << " um" << endl;
	cout.width(15);
	cout << "间隙距离:"; cout.width(10); cout << mGridInfo.interspcaceWidth << " um" << endl;
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

void CBCTForwardProj::creatOutputFolder()
{
	// 指定要创建的多层文件夹路径
	filesystem::path folder_path = mFilePath.outputFolder;

	// 检查是否是目录
	//if (filesystem::is_directory(folder_path)) // 删除整个文件夹，好像不需要判断是否是目录，不存在创建即可。若是删除目录下的文件，应该需要先判断此目录存在否
	//{
		// 检查文件夹是否已经存在
		if (filesystem::exists(folder_path))
		{
			// 删除已存在文件夹
			_removeOutputFolder();
			//cout << "Folder already exists: " << folder_path << endl;
		}

		// 创建
		try {
			// 创建多层文件夹
			filesystem::create_directories(folder_path);
			//cout << "Folder created: " << folder_path << endl;
		}
		catch (const std::exception& e)
		{
			cerr << "Error creating directories(" << mFilePath.outputFolder << "): " << e.what() << endl;
		}

	//}
	/*else
	{
		cerr << "Error: (" << mFilePath.outputFolder << ") Invalid directory path." << endl;
		system("pause");
		exit(1);
		return;
	}*/

}

void CBCTForwardProj::_removeOutputFolder()
{
	filesystem::path folder_path = mFilePath.outputFolder;

	try 
	{
		filesystem::remove_all(folder_path);
	}
	catch (const std::exception& e) {
		cerr << "Error removing directory(" << mFilePath.outputFolder << "): " << e.what() << endl;
		system("pause");
		exit(1);
		return;
	}
}
