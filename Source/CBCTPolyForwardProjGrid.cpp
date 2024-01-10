#include <fstream>
#include <iostream>


#include "CBCTPolyForwardProjGrid.h"

using namespace std;

CBCTPolyForwardProjGrid::CBCTPolyForwardProjGrid(CTScanParas inCTScanParas, GridInfo inGridInfo, FilePath inFileInfo, const PhantomMaterial inPhantomMaterial)
{
	specIndex = 0;
	mCTScanParas = { 0 };
	mCTScanSystemInfo = { 0 };
	mFilePath = {  };
	mGridInfo = { 0 };
	mGridAndDetectorSystem = { 0 };
	_mPhantomMaterial = { 0 };

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

	mCTScanParas.spectrumStep = 1;// inCTScanParas.spectrumStep;
	mCTScanParas.specEnergyNum = inCTScanParas.specEnergyNum;

	mCTScanParas.I0Val = inCTScanParas.I0Val;
	mCTScanParas.focalSpotSize = inCTScanParas.focalSpotSize;

	mCTScanParas.mScintilltorInfo.scintillatorName = inCTScanParas.mScintilltorInfo.scintillatorName;
	mCTScanParas.mScintilltorInfo.scintillatorThickness = inCTScanParas.mScintilltorInfo.scintillatorThickness;
	mCTScanParas.mScintilltorInfo.scintillatorThicknessErr = inCTScanParas.mScintilltorInfo.scintillatorThicknessErr;
	mCTScanParas.mScintilltorInfo.scintillatorDensity = inCTScanParas.mScintilltorInfo.scintillatorDensity;
	mCTScanParas.mScintilltorInfo.detResponseFactor = inCTScanParas.mScintilltorInfo.detResponseFactor;

	// 栅信息初始化
	mGridInfo.FD = inGridInfo.FD;
	mGridInfo.h = inGridInfo.h;
	mGridInfo.interspcaceWidth = inGridInfo.interspcaceWidth;
	mGridInfo.gridStripsWidth = inGridInfo.gridStripsWidth;
	mGridInfo.interspaceMaterial = inGridInfo.interspaceMaterial;
	mGridInfo.gridStripMaterial = inGridInfo.gridStripMaterial;
	mGridInfo.gridStripDensity = inGridInfo.gridStripDensity;
	mGridInfo.gridStripLinearAtten = inGridInfo.gridStripLinearAtten;
	mGridInfo.interspaceLinearAtten = inGridInfo.interspaceLinearAtten;
	mGridInfo.gridStripMAPath = inGridInfo.gridStripMAPath;
	mGridInfo.interspaceMAPath = inGridInfo.interspaceMAPath;

	// 路径初始化
	mFilePath.spectrumPath = inFileInfo.spectrumPath;
	mFilePath.phantomPath = inFileInfo.phantomPath;
	mFilePath.phantomMAttenPath = inFileInfo.phantomMAttenPath;
	mFilePath.gridPath = inFileInfo.gridPath;
	mFilePath.outputFolder = inFileInfo.outputFolder;
	mFilePath.scintillatorPath = inFileInfo.scintillatorPath;
	mFilePath.convKernelGridPath = inFileInfo.convKernelGridPath;

	// 模体材质信息初始化
	_mPhantomMaterial.materialNum = inPhantomMaterial.materialNum;
	_mPhantomMaterial.structureMarkPath = inPhantomMaterial.structureMarkPath;
	_mPhantomMaterial.density = inPhantomMaterial.density;
	_mPhantomMaterial.massAttenPath = inPhantomMaterial.massAttenPath;
}

CBCTPolyForwardProjGrid::~CBCTPolyForwardProjGrid()
{
	DELETEARR(h_mForwardProj.grid);
	DELETEARR(h_mForwardProj.I);
	DELETEARR(h_mForwardProj.I0);
	DELETEARR(h_mForwardProj.phantom);
	DELETEARR(h_mForwardProj.phantomMassAtten);
	DELETEARR(h_mForwardProj.spectrumNormal);
	DELETEARR(h_mForwardProj.IAbsorb);
	DELETEARR(h_mForwardProj.detResponse);
	DELETEARR(h_mForwardProj.gridLinearAtten);
	DELETEARR(h_mForwardProj.interspaceLinearAtten);
	DELETEARR(h_mForwardProj.scintillatorLineAtten);
	DELETEARR(h_mForwardProj.scintillatorPerThickness);
	DELETEARR(h_mForwardProj.proj);
	DELETEARR(h_mForwardProj.foSpOffsetU);
	DELETEARR(h_mForwardProj.foSpOffsetV);

	DELETEARR(I);
	DELETEARR(IPolyenergetic);
	DELETEARR(IAbsorb);
	DELETEARR(IPolyenergeticAbsorb);
	DELETEARR(aDetResponse);
}



void CBCTPolyForwardProjGrid::computePolyForwProj()
{
	cout << "加滤线栅的正过程成像开始 ====>>>>" << endl;
	
	// 初始化主存
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // 初始化I为0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果


	computeParas();

	// 显示成像信息
	showScanSystemInfo();

	// 读取能谱信息
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// 模体构造
	mPhantomObject = new PhantomObject(mCTScanParas, _mPhantomMaterial);
	mPhantomObject->readData();

	// 计算闪烁体厚度
	h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	computePerScinThinckness();

	// 读取闪烁体质量衰减系数
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();

	// 读取栅信息
	h_mForwardProj.gridLinearAtten = new float[mCTScanParas.specEnergyNum];
	h_mForwardProj.interspaceLinearAtten = new float[mCTScanParas.specEnergyNum];
	readGridMassAttu();

	initGridInfo();
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	mGrid = new Grid(h_mForwardProj.grid, mGridAndDetectorSystem);
	mGrid->defGrid();


	initDeviceGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		//readPhantom();
		mPhantomObject->updataLineAtten(h_mForwardProj.phantom, specIndex);

		updateDetResponse();

		mGrid->updataGrid(h_mForwardProj.gridLinearAtten[i], h_mForwardProj.interspaceLinearAtten[i]);

		showProcessInfo();

		forwardProjGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	convertIandIAbsorbDt();

	creatOutputFolder();
	saveIandIAbsorb();
	saveI0();
	

	mGrid->getPeriodIncex(gridPeriodIndex);


	//scatterSimulGrid();

	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	DELETE(mGrid);

	cout << "加滤线栅的正过程成像完成!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTPolyForwardProjGrid::testI0()
{
	cout << "I0测试计算开始 ====>>>>" << endl;

	// 初始化主存
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();


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

	// 读取栅信息
	h_mForwardProj.gridLinearAtten = new float[mCTScanParas.specEnergyNum];
	h_mForwardProj.interspaceLinearAtten = new float[mCTScanParas.specEnergyNum];
	readGridMassAttu();

	initGridInfo();
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	mGrid = new Grid(h_mForwardProj.grid, mGridAndDetectorSystem);
	mGrid->defGrid();


	//initDeviceGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		//readPhantom();

		updateDetResponse();

		mGrid->updataGrid(h_mForwardProj.gridLinearAtten[i], h_mForwardProj.interspaceLinearAtten[i]);

		//showProcessInfo();

		//forwardProjGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		//addIandIAbsorb();
		addI0();

		specIndex++;
	}
	//freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	//IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	//IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	//convertIandIAbsorbDt();
	//saveIandIAbsorb();
	saveI0();


	mGrid->getPeriodIncex(gridPeriodIndex);


	//scatterSimulGrid();

	//h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	//computeProj();
	//saveProj();

	DELETE(mGrid);

	cout << "I0测试完成!" << endl
		<< "-------------------------------------" << endl;
}

// 初始化栅信息
void CBCTPolyForwardProjGrid::initGridInfo()
{
	mGridAndDetectorSystem.angleNums = mCTScanParas.projNum;
	mGridAndDetectorSystem.dNums = mCTScanParas.dNumU;
	mGridAndDetectorSystem.dSize = mCTScanParas.dSize * 1000; // 单位换算 mm -> um
	
	mGridAndDetectorSystem.FD = mGridInfo.FD;
	mGridAndDetectorSystem.h = mGridInfo.h;
	mGridAndDetectorSystem.gridStripsWidth = mGridInfo.gridStripsWidth;
	mGridAndDetectorSystem.interspcaceWidth = mGridInfo.interspcaceWidth;
	mGridAndDetectorSystem.gridStripMaterial = mGridInfo.gridStripMaterial;
	mGridAndDetectorSystem.interspaceMaterial = mGridInfo.interspaceMaterial;
	mGridAndDetectorSystem.gridStripDensity = mGridInfo.gridStripDensity;
	mGridAndDetectorSystem.interspaceDensity = mGridInfo.interspaceDensity;
	mGridAndDetectorSystem.gridStripMaterial = mGridInfo.gridStripMAPath;
	mGridAndDetectorSystem.interspaceMAPath = mGridInfo.interspaceMAPath;

	mGridAndDetectorSystem.mI0Photons = mCTScanParas.I0Val;
	mGridAndDetectorSystem.gridStripLinearAtten = mGridInfo.gridStripLinearAtten;
	mGridAndDetectorSystem.interspaceLinearAtten = mGridInfo.interspaceLinearAtten;

}

void CBCTPolyForwardProjGrid::saveIandIAbsorb()
{
	// 保存探测器接收的光强数据
	string tempPath = mFilePath.outputFolder + "/PhoImgPro_Grid_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "x" + to_string(mCTScanParas.projNum) + "_uint32.raw";

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);
	ofs.write((const char*)IPolyenergetic, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(uint32));
	ofs.close();


	// 保存物体吸收的光强
	tempPath = mFilePath.outputFolder + "/PhoImgAbsorb_Grid_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "x" + to_string(mCTScanParas.projNum) + "_uint32.raw";

	ofs.open(tempPath, ios::out | ios::binary);
	ofs.write((const char*)IPolyenergeticAbsorb, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(uint32));
	ofs.close();
}

// 不同能量下的I0累加
void CBCTPolyForwardProjGrid::addI0()
{
	for (size_t i = 0; i < mCTScanParas.dNumV; i++)
	{
		for (size_t j = 0; j < mCTScanParas.dNumU; j++)
		{
			h_mForwardProj.I0[i * mCTScanParas.dNumU + j] += mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * h_mForwardProj.grid[j] * h_mForwardProj.detResponse[i * mCTScanParas.dNumU + j];
		}
	}
}

void CBCTPolyForwardProjGrid::addI0NoResponse()
{
	for (size_t i = 0; i < mCTScanParas.dNumV; i++)
	{
		for (size_t j = 0; j < mCTScanParas.dNumU; j++)
		{
			h_mForwardProj.I0[i * mCTScanParas.dNumU + j] += mCTScanParas.I0Val * mCTScanSystemInfo.spectrumVal * h_mForwardProj.grid[j];
		}
	}
}

void CBCTPolyForwardProjGrid::saveI0()
{
	// 保存I0
	string tempPath = mFilePath.outputFolder + "/BgImgPro_Grid_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "_float.raw";

	//deleteFile(tempPath);

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);
	ofs.write((const char*)h_mForwardProj.I0, mCTScanParas.dNumU * mCTScanParas.dNumV * sizeof(float));
	ofs.close();
}

void CBCTPolyForwardProjGrid::readPhantom()
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

void CBCTPolyForwardProjGrid::showProcessInfo()
{
	cout << endl;
	cout << (specIndex + 1) * mCTScanParas.spectrumStep << "keV正过程计算中 =====>>>>>" << endl;
	cout << "归一化能谱信息: " << h_mForwardProj.spectrumNormal[specIndex] << endl;
	//cout << "模体质量衰减系数: " << h_mForwardProj.phantomMassAtten[specIndex] << " cm^2/g" << endl;
	cout << "闪烁体线性衰减系数: " << h_mForwardProj.scintillatorLineAtten[specIndex] << " 1/mm" << "\t";
	cout << "滤线栅线性衰减系数: " << h_mForwardProj.gridLinearAtten[specIndex] << " 1/mm" << endl;
}

void CBCTPolyForwardProjGrid::showProcessInfoNoResponse()
{
	cout << endl;
	cout << (specIndex + 1) * mCTScanParas.spectrumStep << "keV正过程计算中 =====>>>>>" << endl;
	cout << "归一化能谱信息: " << h_mForwardProj.spectrumNormal[specIndex] << endl;
	//cout << "模体质量衰减系数: " << h_mForwardProj.phantomMassAtten[specIndex] << " cm^2/g" << endl;
	//cout << "闪烁体线性衰减系数: " << h_mForwardProj.scintillatorLineAtten[specIndex] << " 1/mm" << "\t";
	cout << "滤线栅线性衰减系数: " << h_mForwardProj.gridLinearAtten[specIndex] << " 1/mm" << endl;
}


void CBCTPolyForwardProjGrid::scatterSimulGrid()
{
}

void CBCTPolyForwardProjGrid::saveProj()
{
	// 保存Proj
	string tempPath = mFilePath.outputFolder + "/Projection_Grid_" + to_string(mCTScanParas.dNumU) + "x" + to_string(mCTScanParas.dNumV) + "x" + to_string(mCTScanParas.projNum) + "_float.raw";

	ofstream ofs;
	ofs.open(tempPath, ios::out | ios::binary);
	ofs.write((const char*)h_mForwardProj.proj, mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum * sizeof(float));
	ofs.close();
}
