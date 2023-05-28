#include "CBCTSinMatPolyForwardProjGrid.h"

#include <iostream>
#include <fstream>

using namespace std;






CBCTSinMatPolyForwardProjGrid::~CBCTSinMatPolyForwardProjGrid()
{
}

void CBCTSinMatPolyForwardProjGrid::computePolyForwProj()
{
	cout << "������դ�������̳���ʼ ====>>>>" << endl;

	// ��ʼ������
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // ��ʼ��IΪ0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��
	
	aDetResponse = new float[mCTScanParas.specEnergyNum];       // ���ĳһ̽Ԫ��Ӧ

	computeParas();

	// ��ʾ������Ϣ
	showScanSystemInfo();
	showGridmInfo();

	// ��ȡ������Ϣ
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// ��ȡ������ģ������˥��ϵ��
	h_mForwardProj.phantomMassAtten = new float[mCTScanParas.specEnergyNum];
	readPhantomMassAtten();
	// ��ȡ������ģ���ܶ�
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ];
	readPhantom();

	// ������˸����
	h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	computePerScinThinckness();

	// ��ȡ��˸������˥��ϵ�� ������������˥��ϵ��
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();

	// ��ȡդ��Ϣ
	h_mForwardProj.gridLineAtten = new float[mCTScanParas.specEnergyNum];
	readGridMassAttu();
	// ��ʼ��դ��Ϣ
	initGridInfo();
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	mGrid = new Grid(h_mForwardProj.grid, mGridAndDetectorSystem);
	mGrid->defGrid();

	// ��ʼ��GPU
	cout << endl << "�ܶ�ͼ��ͶӰ������ ===>>>" << endl;
	initDeviceSinMatGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// ��õ�ǰ������Ӧ�ĸ��ʣ������в�ͬ�����Ķ�Ӧֵ
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// ���ߵ�����ģ�������˥��ϵ��
		mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;

		// ����̽������Ӧ
		updateDetResponse();

		// ���ĳһ̽Ԫ����Ӧ
		getADetResponse();

		// ����դ
		mGrid->updataGrid(h_mForwardProj.gridLineAtten[i]);

		// ��ʾ����
		showProcessInfo();

		forwardSinMatProjGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// ����
		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// ��������ת�������ݱ���
	convertIandIAbsorbDt();
	saveIandIAbsorb();
	saveI0();

	//showI0Grid();

	// �������������Ϣ
	mGrid->getPeriodIncex(gridPeriodIndex);


	//scatterSimulGrid();

	// �����м�̽Ԫ��Ӧ
	saveADetResponse();

	DELETE(mGrid);

	// ����proj������
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "������դ�������̳������!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTSinMatPolyForwardProjGrid::computePolyForwProjNoResponse()
{
	cout << "������դ��̽������Ӧ�������̳���ʼ ====>>>>" << endl;

	// ��ʼ������
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // ��ʼ��IΪ0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	//h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��

	//aDetResponse = new float[mCTScanParas.specEnergyNum];       // ���ĳһ̽Ԫ��Ӧ

	computeParas();

	// ��ʾ������Ϣ
	showScanSystemInfo();
	showGridmInfo();

	// ��ȡ������Ϣ
	h_mForwardProj.spectrumNormal = new float[mCTScanParas.specEnergyNum];
	readSpecrtumNorm();

	// ��ȡ������ģ������˥��ϵ��
	h_mForwardProj.phantomMassAtten = new float[mCTScanParas.specEnergyNum];
	readPhantomMassAtten();
	// ��ȡ������ģ���ܶ�
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ];
	readPhantom();

	// ������˸����
	//h_mForwardProj.scintillatorPerThickness = new float[mCTScanParas.dNumU * mCTScanParas.dNumV];
	//computePerScinThinckness();

	// ��ȡ��˸������˥��ϵ��
	//h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	//readScintillatorMassAttu();

	// ��ȡդ��Ϣ
	h_mForwardProj.gridLineAtten = new float[mCTScanParas.specEnergyNum];
	readGridMassAttu();
	// ��ʼ��դ��Ϣ
	initGridInfo();
	h_mForwardProj.grid = new float[mCTScanParas.dNumU];
	mGrid = new Grid(h_mForwardProj.grid, mGridAndDetectorSystem);
	mGrid->defGrid();

	// ��ʼ��GPU
	initDeviceSinMatNoResponseGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// ��õ�ǰ������Ӧ�ĸ��ʣ������в�ͬ�����Ķ�Ӧֵ
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// ���ߵ�����ģ�������˥��ϵ��
		mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;

		// ����̽������Ӧ
		//updateDetResponse();

		// ���ĳһ̽Ԫ����Ӧ
		//getADetResponse();

		// ����դ
		mGrid->updataGrid(h_mForwardProj.gridLineAtten[i]);

		// ��ʾ����
		showProcessInfoNoResponse();

		forwardSinMatNoResponseProjGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// ����
		addIandIAbsorb();
		addI0NoResponse();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// ��������ת�������ݱ���
	convertIandIAbsorbDt();
	saveIandIAbsorb();
	saveI0();

	//showI0Grid();

	// �������������Ϣ
	mGrid->getPeriodIncex(gridPeriodIndex);


	//scatterSimulGrid();

	// �����м�̽Ԫ��Ӧ
	//saveADetResponse();

	DELETE(mGrid);

	// ����proj������
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "������դ��̽������Ӧ�������̳������!" << endl
		<< "-------------------------------------" << endl;
}

// �����ʣ���ȡ�ܶ� g/cn^3
void CBCTSinMatPolyForwardProjGrid::readPhantom()
{
	ifstream ifs;

	ifs.open(mFilePath.phantomPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << "������ģ���ܶ����ݴ�ʧ��!" << endl;
		system("pause");
		exit(0);
	}

	ifs.read((char*)h_mForwardProj.phantom, mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ * sizeof(float));
	ifs.close();
}

// ��ȡģ������˥��ϵ�����ı��ļ�
void CBCTSinMatPolyForwardProjGrid::readPhantomMassAtten()
{
	ifstream ifs;

	ifs.open(mFilePath.phantomMAttenPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << "������ģ���ܶ����ݴ�ʧ��!" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		// ��ȡ��������˥��ϵ��
		h_mForwardProj.phantomMassAtten[ii] = buf;
		ii++;
	}

	ifs.close();
}
