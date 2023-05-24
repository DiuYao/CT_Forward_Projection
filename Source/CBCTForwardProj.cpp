
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
		* mCTScanSystemInfo.dHalfLV / mCTScanParas.sdd;   // Բ����
	// mCTScanSystemInfo.FOVH = mCTScanParas.sod * mCTScanSystemInfo.dHalfLV / sqrt(mCTScanSystemInfo.dHalfLV * mCTScanSystemInfo.dHalfLV + mCTScanParas.sdd * mCTScanParas.sdd);                  // ��ֱ��ҰԲ�뾶

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
	mCTScanSystemInfo.intNum = roundf(2 * mCTScanSystemInfo.FOVR / mCTScanSystemInfo.dx);        // ĳ�������ϵĻ��ֵ����

}

// ��ͬ�����µ�ʣ����Ӻ��������չ��ӷֱ��ۼ�
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

// ʣ����Ӻ��������չ�����������ת�� float ��> uint32
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

// ������˸���ȣ�ģ�⹤�յ��º�Ȳ�һ��
void CBCTForwardProj::computePerScinThinckness()
{
	srand(100);  // ���������
	float tempErr = 0.0f;
	int beg = -mCTScanParas.mScintilltorInfo.scintillatorThicknessErr * 10000;  // �Ŵ�Ϊ����
	int endd = mCTScanParas.mScintilltorInfo.scintillatorThicknessErr * 10000;

	for (size_t i = 0; i < mCTScanParas.dNumU * mCTScanParas.dNumV; i++)
	{
		tempErr = -mCTScanParas.mScintilltorInfo.scintillatorThicknessErr + (rand() % (endd - beg + 1)) / 10000.0f;  // ���������� 
		h_mForwardProj.scintillatorPerThickness[i] = mCTScanParas.mScintilltorInfo.scintillatorThickness + tempErr;
	}
}

// ����̽������Ӧ���������仯
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
		cout << "�ļ���ʧ��" << endl;
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

// ��ȡ����դ������˥��ϵ��(cm^2/g)��ת��Ϊ����˥��ϵ��(1/mm)
void CBCTForwardProj::readGridMassAttu()
{
	ifstream ifs;
	ifs.open(mFilePath.gridPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "դ˥��ϵ���ļ���ʧ��" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		// ��ȡ��������˥��ϵ������Ҫת��Ϊ����˥��ϵ����������λת��Ϊ 1/mm
		h_mForwardProj.gridLineAtten[ii] = buf;
		// ת��
		h_mForwardProj.gridLineAtten[ii] = h_mForwardProj.gridLineAtten[ii] * mGridInfo.rhoGS / 10;    // 1/mm
		ii++;
	}

	ifs.close();
}

// ��ȡ��˸�������˥��ϵ��(cm^2/g)��ת��Ϊ����˥��ϵ��(1/mm)
void CBCTForwardProj::readScintillatorMassAttu()
{
	ifstream ifs;
	ifs.open(mFilePath.scintillatorPath, ios::in);

	if (!ifs.is_open())
	{
		cout << "��˸������˥��ϵ���ļ���ʧ�ܣ�" << endl;
		system("pause");
		exit(0);
	}

	int ii = 0;
	float buf;
	while (ifs >> buf)
	{
		// ��ȡ��������˥��ϵ������Ҫת��Ϊ����˥��ϵ����������λת��Ϊ 1/mm
		h_mForwardProj.scintillatorLineAtten[ii] = buf;
		// ת��
		h_mForwardProj.scintillatorLineAtten[ii] = h_mForwardProj.scintillatorLineAtten[ii] * mCTScanParas.mScintilltorInfo.scintillatorDensity / 10;  // 1/mm
		ii++;
	}

	ifs.close();
}

// ����м�λ��̽Ԫ����Ӧ
void CBCTForwardProj::getADetResponse()
{
	aDetResponse[specIndex] = h_mForwardProj.detResponse[lroundl((mCTScanParas.dNumU + mCTScanParas.dNumV) / 2)];
}


void CBCTForwardProj::saveADetResponse()
{
	string directory = "OutputResult"; // ָ��Ŀ¼��
	string filename = "ADetResponse.txt"; // �ļ���
	
	if (_access(directory.c_str(), 0) == -1)
	{
		// if this folder not exist, create a new one.
		// ����Ŀ¼
		if (_mkdir(directory.c_str()) == -1)
		{
			cerr << "Error: Creating ADetResponse file is fair��" << endl;
		}
		// ���� 0 ��ʾ�����ɹ���-1 ��ʾʧ��
	}
	


	// ���ļ�
	ofstream ofs;
	ofs.open(directory + "/" + filename, ios::out);

	if (!ofs.is_open()) {
		cerr << "Error: Cannot open ADetResponse file" << endl;
		return ;
	}

	// д������
	for (size_t i = 0; i < mCTScanParas.specEnergyNum; i++)
	{
		ofs << aDetResponse[i] << endl;
	}

	// �ر��ļ�
	ofs.close();
}

void CBCTForwardProj::showScanSystemInfo()
{
	cout.setf(ios::left);
	cout << "*******************����ϵͳ��Ϣ******************" << endl;
	cout.width(15);
	cout << "̽�����ߴ�:" << mCTScanParas.dNumU << "��" << mCTScanParas.dNumV << endl;
	cout.width(15);
	cout << "̽�����ֱ���:"; cout.width(8); cout << mCTScanParas.dSize << " mm" << endl;
	cout.width(15);
	cout << "SDD:"; cout.width(8); cout << mCTScanParas.sdd << " mm" << endl;
	cout.width(15);
	cout << "SOD:"; cout.width(8); cout << mCTScanParas.sod << " mm" << endl;
	cout.width(15);
	cout << "FOV(ֱ������):"; cout.width(8); cout << 2 * mCTScanSystemInfo.FOVR << "��" << 2 * mCTScanSystemInfo.FOVH << " mm" << endl;
	cout.width(15);
	cout << "ģ�����ش�С:"; cout << mCTScanSystemInfo.pSizeX << "��" << mCTScanSystemInfo.pSizeY << "��" << mCTScanSystemInfo.pSizeZ << " mm" << endl;
	cout << endl;

	/*cout.setf(ios::left);
	cout << "*******************����ϵͳ��Ϣ******************" << endl;
	cout.width(15);
	cout << "̽�����ߴ�:" << dNums << "x" << dNums << endl;
	cout.width(15);
	cout << "̽�����ֱ���:"; cout.width(8); cout << dSize << " um" << endl;

	cout << endl;

	cout << "*****************�۽�������դ��Ϣ****************" << endl;
	cout.width(15);
	cout << "դ������:" << materialGridStrip << endl;
	cout.width(15);
	cout << "����˥��ϵ��:"; cout.width(10); cout << uGridStrip * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "��϶�����:" << materialGridInterspacer << endl;
	cout.width(15);
	cout << "����˥��ϵ��:"; cout.width(10); cout << uInterspacer * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "����:"; cout.width(10); cout << gridPeriod << " um" << endl;
	cout.width(15);
	cout << "դ�����:"; cout.width(10); cout << leadStripsWidth << " um" << endl;
	cout.width(15);
	cout << "��϶����:"; cout.width(10); cout << leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "���:"; cout.width(10); cout << h / 1000 << " mm" << endl;
	cout.width(15);
	cout << "����:"; cout.width(10); cout << FD / 1000 << " mm" << endl;
	cout << "*************************************************" << endl;*/
}

void CBCTForwardProj::showGridmInfo()
{
	cout << endl;
	cout << "*****************�۽�������դ��Ϣ****************" << endl;
	cout.width(15);
	cout << "դ������:" << mGridInfo.materialGridStrip << endl;
	//cout.width(15);
	//cout << "����˥��ϵ��:"; cout.width(10); cout << uGridStrip * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "��϶�����:" << mGridInfo.materialGridInterspacer << endl;
	//cout.width(15);
	//cout << "����˥��ϵ��:"; cout.width(10); cout << uInterspacer * 1000 << " 1/mm" << endl;
	cout.width(15);
	cout << "����:"; cout.width(10); cout << mGridInfo.leadStripsWidth + mGridInfo.leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "դ�����:"; cout.width(10); cout << mGridInfo.leadStripsWidth << " um" << endl;
	cout.width(15);
	cout << "��϶����:"; cout.width(10); cout << mGridInfo.leadStripsDistance << " um" << endl;
	cout.width(15);
	cout << "���:"; cout.width(10); cout << mGridInfo.h << " mm" << endl;
	cout.width(15);
	cout << "����:"; cout.width(10); cout << mGridInfo.FD << " mm" << endl;
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
