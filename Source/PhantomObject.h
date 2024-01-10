#pragma once

#include "DataType.h"

// ��ȡ�Ľṹģ�����ݣ����ݱ�ǩ�����Ӧ������˥��ϵ��



class PhantomObject
{
public:
	PhantomObject(const CTScanParas mCTScanParas, const PhantomMaterial mPhantomMaterial);
	~PhantomObject();

	void readData();

	// ����ʸ�������˥��ϵ������λ1/mm
	void updataLineAtten(float* lineAttenCoef, size_t specIndex);
	// �����ʸ�������˥��ϵ������ֵ����10ʹ������˥��ϵ����λΪ1/mm
	void updataMassAtten(float& massAttenCoef, size_t specIndex);

	// �����ʼ����ܶ����ݣ���λg/cm^3
	void computeDesity(float* density);

	
private:
	void _readPhantomStructure();
	void _readMassAtten();

private:

	PhantomMaterial _mPhantomMaterial;
	
	uint8* _structureMark;	// ģ��ṹ���
	
	
	float** _massAttenCoef;

	// ��������
	int _specEnergyNum;
	// ģ���С
	size_t _pNumX, _pNumY, _pNumZ;


	

};

