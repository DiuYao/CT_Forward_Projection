#pragma once

#include <direct.h>
#include <corecrt_io.h>
#include <filesystem>

#include "DataType.h"
#include "PhantomObject.h"

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

	
	virtual void testI0() = 0;
	virtual void addStandardI0(float& I0Standard);

	virtual void addI0() = 0;
	virtual void saveI0() = 0;


	virtual void computePerScinThinckness();
	virtual void updateDetResponse();

	virtual void computeFoSpOffset();

	// 读取数据
	virtual void readSpecrtumNorm();
	// 读取滤线栅的质量衰减系数(cm^2/g)，转换为线性衰减系数(1/mm)
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

	// 输出文件夹检测与创建
	virtual void creatOutputFolder();
	// 保存数据前检测是否是空文件夹，非空，删除其中数据
	virtual void _removeOutputFolder();

	// 卷积
	virtual void scatterSimulGrid() = 0;

public:
	PhantomObject* mPhantomObject;

	CTScanSystemInfo mCTScanSystemInfo;
	CTScanParas mCTScanParas;
	Coordinate d_mCoordinate;
	/*Coordinate d_mCoordinate;*/
	PolyForwardProj h_mForwardProj;
	PolyForwardProj d_mForwardProj;
	FilePath mFilePath;

	PhantomMaterial _mPhantomMaterial;

	GridInfo mGridInfo;   // 栅信息

	size_t specIndex;

	float* I;
	uint32* IPolyenergetic;

	float* IAbsorb;
	uint32* IPolyenergeticAbsorb;

	float* IEScatter;

	float* aDetResponse;
};

