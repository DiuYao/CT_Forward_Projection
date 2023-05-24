#pragma once

#include "DataType.h"

class CBCTForwardProj
{
public:
	CBCTForwardProj();
	virtual ~CBCTForwardProj() = 0;

	virtual void computePolyForwProj() = 0;

	virtual void computeParas();
	virtual void initGridInfo() = 0;

	virtual void addIandIAbsorb();
	virtual void convertIandIAbsorbDt();
	virtual void saveIandIAbsorb() = 0;

	
	virtual void addI0() = 0;
	virtual void saveI0() = 0;


	virtual void computePerScinThinckness();
	virtual void updateDetResponse();

	// ��ȡ����
	virtual void readSpecrtumNorm();
	virtual void readGridMassAttu();

	// ��ȡ��˸������˥��ϵ������ת��Ϊ����˥��ϵ��
	virtual void readScintillatorMassAttu();
	//
	virtual void getADetResponse();
	virtual void saveADetResponse();


	virtual void readPhantom() = 0;
	

	// ɨ��ϵͳ��Ϣ
	virtual void showScanSystemInfo();
	// ����դ��Ϣ
	virtual void showGridmInfo();


	// ���������ʾ
	virtual void showProcessInfo() = 0;

	// ����Proj
	virtual void computeProj();
	virtual void saveProj() = 0;

	// ���
	virtual void scatterSimulGrid() = 0;

public:
	CTScanSystemInfo mCTScanSystemInfo;
	CTScanParas mCTScanParas;
	Coordinate d_mCoordinate;
	/*Coordinate d_mCoordinate;*/
	PolyForwardProj h_mForwardProj;
	PolyForwardProj d_mForwardProj;
	FilePath mFilePath;

	GridInfo mGridInfo;   // դ��Ϣ

	size_t specIndex;

	float* I;
	uint32* IPolyenergetic;

	float* IAbsorb;
	uint32* IPolyenergeticAbsorb;

	float* IEScatter;

	float* aDetResponse;
};

