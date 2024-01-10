#include <iostream>
#include <fstream>

#include "CBCTSinMatPolyForwardProjNoGrid.h"

CBCTSinMatPolyForwardProjNoGrid::~CBCTSinMatPolyForwardProjNoGrid()
{
}

void CBCTSinMatPolyForwardProjNoGrid::computePolyForwProj()
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

	aDetResponse = new float[mCTScanParas.specEnergyNum];   // 存放某一探元响应

	computeParas();

	// 显示成像信息
	showScanSystemInfo();

	// 读取能谱信息
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// 读取单材质模体质量衰减系数
	h_mForwardProj.phantomMassAtten = new float[mCTScanParas.specEnergyNum];
	//readPhantomMassAtten();
	// 获得读取单材质模体密度
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ];
	//readPhantom();
	mPhantomObject = new PhantomObject(mCTScanParas, _mPhantomMaterial);
	mPhantomObject->readData();
	mPhantomObject->computeDesity(h_mForwardProj.phantom);

	// 计算闪烁体厚度
	h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	computePerScinThinckness();

	// 读取闪烁体质量衰减系数
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();


	// 初始化 GPU
	cout << endl << "密度图像投影计算中 ===>>>" << endl;
	initDeviceSinMatNoGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// 获得当前能量对应的概率，即谱中不同能量的对应值
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// 获得单材质模体的质量衰减系数
		//mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;   // 除以10保证单位后面计算单位统一
		mPhantomObject->updataMassAtten(mCTScanSystemInfo.phantomMAtten, specIndex);

		// 更新探测器响应
		updateDetResponse();

		// 获得某一探元的响应
		getADetResponse();

		// 显示进程
		showProcessInfo();

		// 
		forwardSinMatProjNoGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// 计算结果叠加
		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	// 数据类型转换，并保存
	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// 数据类型转换和数据保存
	convertIandIAbsorbDt();

	creatOutputFolder();
	saveIandIAbsorb();
	saveI0();


	//scatterSimulGrid();

	// 保存中间探元响应
	saveADetResponse();

	// 计算proj并保存
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "无滤线栅的正过程成像完成!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTSinMatPolyForwardProjNoGrid::computePolyForwProjNoResponse()
{
	cout << "无滤线栅无探元响应的正过程成像开始 ====>>>>" << endl;

	// 初始化主存
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // 初始化I为0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	//h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果

	//aDetResponse = new float[mCTScanParas.specEnergyNum];   // 存放某一探元响应

	computeParas();

	// 显示成像信息
	showScanSystemInfo();

	// 读取能谱信息
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// 读取单材质模体质量衰减系数
	h_mForwardProj.phantomMassAtten = new float[mCTScanParas.specEnergyNum];
	//readPhantomMassAtten();
	// 读取单材质模体密度
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ];
	//readPhantom();
	mPhantomObject = new PhantomObject(mCTScanParas, _mPhantomMaterial);
	mPhantomObject->readData();
	mPhantomObject->computeDesity(h_mForwardProj.phantom);

	// 计算闪烁体厚度
	//h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	//computePerScinThinckness();

	// 读取闪烁体质量衰减系数
	//h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	//readScintillatorMassAttu();


	// 初始化 GPU
	cout << endl << "密度图像投影计算中 ===>>>" << endl;
	initDeviceSinMatNoResponseNoGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// 获得当前能量对应的概率，即谱中不同能量的对应值
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// 或者单材质模体的质量衰减系数
		//mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;   // 除以10保证单位后面计算单位统一
		mPhantomObject->updataMassAtten(mCTScanSystemInfo.phantomMAtten, specIndex);

		// 更新探测器响应
		//updateDetResponse();

		// 获得某一探元的响应
		//getADetResponse();

		// 显示进程
		showProcessInfoNoResponse();

		// 
		forwardSinMatNoResponseProjNoGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// 计算结果叠加
		addIandIAbsorb();
		addI0NoResponse();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	// 数据类型转换，并保存
	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// 数据类型转换和数据保存
	convertIandIAbsorbDt();

	creatOutputFolder();
	saveIandIAbsorb();
	saveI0();


	//scatterSimulGrid();

	// 保存中间探元响应
	//saveADetResponse();

	// 计算proj并保存
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "无滤线栅无探元响应的正过程成像完成!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTSinMatPolyForwardProjNoGrid::computePolyForwProjFoSp()
{
	cout << "无滤线栅有探元响应有焦斑大小的正过程成像开始 ====>>>>" << endl;

	// 初始化主存
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // 初始化I为0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // 存放累加结果

	aDetResponse = new float[mCTScanParas.specEnergyNum];   // 存放某一探元响应

	computeParas();

	// 显示成像信息
	showScanSystemInfo();

	// 计算焦点偏移
	h_mForwardProj.foSpOffsetU = new float[mCTScanParas.projNum];
	h_mForwardProj.foSpOffsetV = new float[mCTScanParas.projNum];
	computeFoSpOffset();

	// 读取能谱信息
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// // 读取单材质模体质量衰减系数
	h_mForwardProj.phantomMassAtten = new float[mCTScanParas.specEnergyNum];
	//readPhantomMassAtten();
	// 读取单材质模体密度
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ];
	readPhantom();

	// 计算闪烁体厚度
	h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	computePerScinThinckness();

	// 读取闪烁体质量衰减系数
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();


	// 初始化 GPU
	cout << endl << "密度图像投影计算中 ===>>>" << endl;
	initDeviceSinMatFoSpSiNoGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// 获得当前能量对应的概率，即谱中不同能量的对应值
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// 或者单材质模体的质量衰减系数
		//mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;   // 除以10保证单位后面计算单位统一
		mPhantomObject->updataMassAtten(mCTScanSystemInfo.phantomMAtten, specIndex);

		// 更新探测器响应
		updateDetResponse();

		// 获得某一探元的响应
		getADetResponse();

		// 显示进程
		showProcessInfo();

		// 
		forwardSinMatProjNoGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// 计算结果叠加
		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	// 数据类型转换，并保存
	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// 数据类型转换和数据保存
	convertIandIAbsorbDt();
	saveIandIAbsorb();
	saveI0();

	//scatterSimulGrid();

	// 保存中间探元响应
	saveADetResponse();

	// 计算proj并保存
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "无滤线栅有探元响应有焦斑大小的正过程成像完成!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTSinMatPolyForwardProjNoGrid::readPhantom()
{
	ifstream ifs;

	ifs.open(mFilePath.phantomPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << "单材质模体密度数据打开失败!" << endl;
		system("pause");
		exit(0);
	}

	ifs.read((char*)h_mForwardProj.phantom, mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ * sizeof(float));
	ifs.close();
}

void CBCTSinMatPolyForwardProjNoGrid::readPhantomMassAtten()
{
	ifstream ifs;

	ifs.open(mFilePath.phantomMAttenPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << "单材质模体质量衰减系数数据打开失败!" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf && ii < mCTScanParas.specEnergyNum)
	{
		// 读取的是质量衰减系数
		h_mForwardProj.phantomMassAtten[ii] = buf;
		ii++;
	}

	ifs.close();
}
