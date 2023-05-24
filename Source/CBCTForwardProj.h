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

	// 读取数据
	virtual void readSpecrtumNorm();
	virtual void readGridMassAttu();

	// 读取闪烁体质量衰减系数，并转化为线性衰减系数
	virtual void readScintillatorMassAttu();
	//
	virtual void getADetResponse();
	virtual void saveADetResponse();


	virtual void readPhantom() = 0;
	

	// 扫描系统信息
	virtual void showScanSystemInfo();
	// 滤线栅信息
	virtual void showGridmInfo();


	// 程序进度显示
	virtual void showProcessInfo() = 0;

	// 计算Proj
	virtual void computeProj();
	virtual void saveProj() = 0;

	// 卷积
	virtual void scatterSimulGrid() = 0;

public:
	CTScanSystemInfo mCTScanSystemInfo;
	CTScanParas mCTScanParas;
	Coordinate d_mCoordinate;
	/*Coordinate d_mCoordinate;*/
	PolyForwardProj h_mForwardProj;
	PolyForwardProj d_mForwardProj;
	FilePath mFilePath;

	GridInfo mGridInfo;   // 栅信息

	size_t specIndex;

	float* I;
	uint32* IPolyenergetic;

	float* IAbsorb;
	uint32* IPolyenergeticAbsorb;

	float* IEScatter;

	float* aDetResponse;
};

