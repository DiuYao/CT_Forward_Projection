#include <iostream>
#include <fstream>

#include "CBCTSinMatPolyForwardProjNoGrid.h"

CBCTSinMatPolyForwardProjNoGrid::~CBCTSinMatPolyForwardProjNoGrid()
{
}

void CBCTSinMatPolyForwardProjNoGrid::computePolyForwProj()
{
	cout << "������դ�������̳���ʼ ====>>>>" << endl;

	// ��ʼ������
	h_mForwardProj.phantom = new float[mCTScanParas.pNumX * mCTScanParas.pNumY * mCTScanParas.pNumZ]();
	h_mForwardProj.I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum](); // ��ʼ��IΪ0
	h_mForwardProj.IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	h_mForwardProj.I0 = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();
	h_mForwardProj.detResponse = new float[mCTScanParas.dNumU * mCTScanParas.dNumV]();

	I = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��
	IAbsorb = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();  // ����ۼӽ��

	aDetResponse = new float[mCTScanParas.specEnergyNum];   // ���ĳһ̽Ԫ��Ӧ

	computeParas();

	// ��ʾ������Ϣ
	showScanSystemInfo();

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

	// ��ȡ��˸������˥��ϵ��
	h_mForwardProj.scintillatorLineAtten = new float[mCTScanParas.specEnergyNum];
	readScintillatorMassAttu();


	// ��ʼ�� GPU
	initDeviceNoGrid(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

	specIndex = 0;
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		// ��õ�ǰ������Ӧ�ĸ��ʣ������в�ͬ�����Ķ�Ӧֵ
		mCTScanSystemInfo.spectrumVal = h_mForwardProj.spectrumNormal[i];
		// ���ߵ�����ģ�������˥��ϵ��
		mCTScanSystemInfo.phantomMAtten = h_mForwardProj.phantomMassAtten[i] / 10.0f;   // ����10��֤��λ������㵥λͳһ

		// ����̽������Ӧ
		updateDetResponse();

		// ���ĳһ̽Ԫ����Ӧ
		getADetResponse();

		// ��ʾ����
		showProcessInfo();

		// 
		forwardSinMatProjNoGridGPU(d_mForwardProj, d_mCoordinate, mCTScanSystemInfo, mCTScanParas, h_mForwardProj);

		// ����������
		addIandIAbsorb();
		addI0();

		specIndex++;
	}
	freeDeviceMemory(d_mForwardProj, d_mCoordinate);

	// ��������ת����������
	IPolyenergetic = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];
	IPolyenergeticAbsorb = new uint32[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum];

	// ��������ת�������ݱ���
	convertIandIAbsorbDt();
	saveIandIAbsorb();
	saveI0();


	//scatterSimulGrid();

	// �����м�̽Ԫ��Ӧ
	saveADetResponse();

	// ����proj������
	h_mForwardProj.proj = new float[mCTScanParas.dNumU * mCTScanParas.dNumV * mCTScanParas.projNum]();
	computeProj();
	saveProj();

	cout << "������դ�������̳������!" << endl
		<< "-------------------------------------" << endl;
}

void CBCTSinMatPolyForwardProjNoGrid::computePolyForwProjNoResponse()
{
	cout << "�պ�����" << endl;
}

void CBCTSinMatPolyForwardProjNoGrid::readPhantom()
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

void CBCTSinMatPolyForwardProjNoGrid::readPhantomMassAtten()
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