
#include <iostream>
#include <fstream>


#include "CBCTPolyForwardProjNoGrid.h"


using namespace std;

CBCTPolyForwardProjNoGrid::CBCTPolyForwardProjNoGrid(CTScanParas inCTScanParas, FilePath inFileInfo)
{
	specIndex = 0;
	mCTScanParas = { 0 };
	mCTScanSystemInfo = { 0 };
	mFilePath = {  };
	mGridInfo = { 0 };

	h_mForwardProj = { nullptr };
	d_mForwardProj = { nullptr };
	d_mCoordinate = { nullptr };

	// 扫描参数初始化
	mCTScanParas.pNumX = inCTScanParas.pNumX;
	mCTScanParas.pNumY = inCTScanParas.pNumY;
	mCTScanParas.pNumZ = inCTScanParas.pNumZ;
	mCTScanParas.pLengthX = inCTScanParas.pLengthX;
	mCTScanParas.pLengthY = inCTScanParas.pLengthY;
	mCTScanParas.pLengthZ = inCTScanParas.pLengthZ;

	mCTScanParas.dNumU = inCTScanParas.dNumU;
	mCTScanParas.dNumV = inCTScanParas.dNumV;
	mCTScanParas.dSize = inCTScanParas.dSize;

	mCTScanParas.sdd = inCTScanParas.sdd;
	mCTScanParas.sod = inCTScanParas.sod;
	mCTScanParas.projNum = inCTScanParas.projNum;
	mCTScanParas.rotatedDirection = inCTScanParas.rotatedDirection;

	mCTScanParas.spectrumStep = inCTScanParas.spectrumStep;
	mCTScanParas.specEnergyNum = inCTScanParas.specEnergyNum;

	mCTScanParas.I0Val = inCTScanParas.I0Val;
	mCTScanParas.mScintilltorInfo.scintillatorName = inCTScanParas.mScintilltorInfo.scintillatorName;
	mCTScanParas.mScintilltorInfo.scintillatorThickness = inCTScanParas.mScintilltorInfo.scintillatorThickness;
	mCTScanParas.mScintilltorInfo.scintillatorThicknessErr = inCTScanParas.mScintilltorInfo.scintillatorThicknessErr;
	mCTScanParas.mScintilltorInfo.scintillatorDensity = inCTScanParas.mScintilltorInfo.scintillatorDensity;
	mCTScanParas.mScintilltorInfo.detResponseFactor = inCTScanParas.mScintilltorInfo.detResponseFactor;

	// 路径初始化
	mFilePath.spectrumPath = inFileInfo.spectrumPath;
	mFilePath.phantomPath = inFileInfo.phantomPath;
	mFilePath.phantomMAttenPath = inFileInfo.phantomMAttenPath;
	mFilePath.gridPath = inFileInfo.gridPath;
	mFilePath.IEPath = inFileInfo.IEPath;
	mFilePath.scintillatorPath = inFileInfo.scintillatorPath;
	mFilePath.convKernelGridPath = inFileInfo.convKernelGridPath;
}

CBCTPolyForwardProjNoGrid::~CBCTPolyForwardProjNoGrid()
{
	DELETEARR(h_mForwardProj.grid);
	DELETEARR(h_mForwardProj.I);
	DELETEARR(h_mForwardProj.I0);
	DELETEARR(h_mForwardProj.phantom);
	DELETEARR(h_mForwardProj.phantomMassAtten);
	DELETEARR(h_mForwardProj.spectrumNormal);
	DELETEARR(h_mForwardProj.IAbsorb);
	DELETEARR(h_mForwardProj.detResponse);
	DELETEARR(h_mForwardProj.gridLineAtten);
	DELETEARR(h_mForwardProj.scintillatorLineAtten);
	DELETEARR(h_mForwardProj.scintillatorPerThickness);
	DELETEARR(h_mForwardProj.proj);

	DELETEARR(I);
	DELETEARR(IPolyenergetic);
	DELETEARR(IAbsorb);
	DELETEARR(IPolyenergeticAbsorb);
	DELETEARR(aDetResponse);
}

void CBCTPolyForwardProjNoGrid::computePolyForwProj()
{
	cout << "无滤线栅的正过程成像开始 ====>>>>" << endl;

	// 初始化主存
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // 初始化I为0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果


	computeParas();

	// 读取能谱信息
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// 计算闪烁体厚度
	h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	computePerScinThinckness();

	// 读取闪烁体质量衰减系数
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();


	// 初始化 GPU
	initDeviceNoGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// 获得当前能量对应的概率，即谱中不同能量的对应值
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];

		// 读取模体数据
		readPhantom();

		// 更新探测器响应
		updateDetResponse();

		showProcessInfo();

		// GPU计算
		forwardProjNoGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// 计算结果叠加
		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	// 释放GPU内存
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	// 数据类型转换，并保存
	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	convertIandIAbsorbDt();
	saveIandIAbsorb();
	saveI0();


	//scatterSimulGrid();

	//
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "加滤线栅的正过程成像完成!"
		<< "-------------------------------------" << endl;
}

void CBCTPolyForwardProjNoGrid::initGridInfo()
{

}

void CBCTPolyForwardProjNoGrid::saveIandIAbsorb()
{
	// 保存探测器接收的光强
	string tempPath = mFilePath.IEPath + "/PhoImgPro_NoGrid_Poly_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "x" + to_string(mCTScanParas.projNum) + "_uint32.raw";

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);

	ofs.write((const char*)IPolyenergetic, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(uint32));
	ofs.close();


	// 保存物体吸收的光强
	tempPath = mFilePath.IEPath + "/PhoImgAbsorb_NoGrid_Poly_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "x" + to_string(mCTScanParas.projNum) + "_uint32.raw";

	ofs.open(tempPath, ios::out | ios::binary);

	ofs.write((const char*)IPolyenergeticAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(uint32));
	ofs.close();
}

void CBCTPolyForwardProjNoGrid::addI0()
{
	for (size_t i = 0; i < mCTScanParas.dNumV; i++)
	{
		for (size_t j = 0; j < mCTScanParas.dNumU; j++)
		{
			h_mForwardProj.I0[i * mCTScanParas.dNumU + j] += mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * h_mForwardProj.detResponse[i * mCTScanParas.dNumU + j];
		}
	}
}

void CBCTPolyForwardProjNoGrid::addI0NoResponse()
{
	for (size_t i = 0; i < mCTScanParas.dNumV; i++)
	{
		for (size_t j = 0; j < mCTScanParas.dNumU; j++)
		{
			h_mForwardProj.I0[i * mCTScanParas.dNumU + j] += mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal;
		}
	}
}

void CBCTPolyForwardProjNoGrid::saveI0()
{
	// 保存I0
	string tempPath = mFilePath.IEPath + "/BgImgPro_NoGrid_Poly_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "_float.raw";

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);

	ofs.write((const char*)h_mForwardProj.I0, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
	ofs.close();
}

void CBCTPolyForwardProjNoGrid::readPhantom()
{
	ifstream ifs;
	string tempPath;

	tempPath = mFilePath.phantomPath + "_" + to_string((specIndex + 1) * mCTScanParas.spectrumStep) + "keV.raw";

	ifs.open(tempPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << (specIndex + 1) * mCTScanParas.spectrumStep << "keV的模体数据打开失败" << endl;
		system("pause");
		exit(0);
	}

	ifs.read((char*)h_mForwardProj.phantom, mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ * sizeof(float));
	ifs.close();
}

void CBCTPolyForwardProjNoGrid::showProcessInfo()
{
	cout << endl;
	cout << (specIndex + 1) * mCTScanParas.spectrumStep << "keV正过程计算中 =====>>>>>" << endl;
	cout << "归一化能谱信息: " << h_mForwardProj.spectrumNormal[specIndex] << endl;
	cout << "模体质量衰减系数: " << h_mForwardProj.phantomMassAtten[specIndex] << " cm^2/g" << endl;
	cout << "闪烁体线性衰减系数: " << h_mForwardProj.scintillatorLineAtten[specIndex] << " 1/mm" << endl;
}

void CBCTPolyForwardProjNoGrid::showProcessInfoNoResponse()
{
	cout << endl;
	cout << (specIndex + 1) * mCTScanParas.spectrumStep << "keV正过程计算中 =====>>>>>" << endl;
	cout << "归一化能谱信息: " << h_mForwardProj.spectrumNormal[specIndex] << endl;
	cout << "模体质量衰减系数: " << h_mForwardProj.phantomMassAtten[specIndex] << " cm^2/g" << endl;
}

void CBCTPolyForwardProjNoGrid::scatterSimulGrid()
{
}

void CBCTPolyForwardProjNoGrid::saveProj()
{
	// 保存Proj
	string tempPath = mFilePath.IEPath + "/Projection_NoGrid_Poly_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "_float.raw";

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);

	ofs.write((const char*)h_mForwardProj.proj, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	ofs.close();
}

