#pragma once

#include "DataType.h"

// 读取的结构模体数据，根据标签计算对应的线性衰减系数



class PhantomObject
{
public:
	PhantomObject(const CTScanParas mCTScanParas, const PhantomMaterial mPhantomMaterial);
	~PhantomObject();

	void readData();

	// 多材质更新线性衰减系数，单位1/mm
	void updataLineAtten(float* lineAttenCoef, size_t specIndex);
	// 单材质更新质量衰减系数，数值除以10使得线性衰减系数单位为1/mm
	void updataMassAtten(float& massAttenCoef, size_t specIndex);

	// 单材质计算密度数据，单位g/cm^3
	void computeDesity(float* density);

	
private:
	void _readPhantomStructure();
	void _readMassAtten();

private:

	PhantomMaterial _mPhantomMaterial;
	
	uint8* _structureMark;	// 模体结构标记
	
	
	float** _massAttenCoef;

	// 能量个数
	int _specEnergyNum;
	// 模体大小
	size_t _pNumX, _pNumY, _pNumZ;


	

};

