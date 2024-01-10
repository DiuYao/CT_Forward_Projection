#include <fstream>
#include <iostream>

#include "PhantomObject.h"

using namespace std;

PhantomObject::PhantomObject(CTScanParas mCTScanParas, PhantomMaterial mPhantomMaterial)
{
	_structureMark = nullptr;
	_specEnergyNum = mCTScanParas.specEnergyNum;
	_pNumX = mCTScanParas.pNumX;
	_pNumY = mCTScanParas.pNumY;
	_pNumZ = mCTScanParas.pNumZ;
	_mPhantomMaterial.density = mPhantomMaterial.density;
	_mPhantomMaterial.massAttenPath = mPhantomMaterial.massAttenPath;
	_mPhantomMaterial.materialNum = mPhantomMaterial.materialNum;
	_mPhantomMaterial.structureMarkPath = mPhantomMaterial.structureMarkPath;
}

PhantomObject::~PhantomObject()
{
	DELETEARR(_structureMark);
	// �ͷŶ�ά����
	for (int i = 0; i < _specEnergyNum; i++) {
		DELETEARR(_massAttenCoef[i]);
	}
	DELETEARR(_massAttenCoef);
}

void PhantomObject::readData()
{
	_readPhantomStructure();
	_readMassAtten();
}

// ����ʸ�������˥��ϵ������λ1/mm
void PhantomObject::updataLineAtten(float* lineAttenCoef, size_t specIndex)
{
	//readData();
	for (size_t i = 0; i < _pNumX * _pNumY * _pNumZ; i++)
	{
		for (uint8 j = 0; j < _mPhantomMaterial.materialNum; j++)
		{
			if (_structureMark[i] == j+1)
			{
				// ��������˥��ϵ����������λ��Ϊ1/mm
				lineAttenCoef[i] = _massAttenCoef[specIndex][j] * _mPhantomMaterial.density[j] / 10.f;
				break;
			}
			else if(_structureMark[i] == 0)
			{
				lineAttenCoef[i] = 0;
				break;
			}
		}
	}
}

void PhantomObject::updataMassAtten(float& massAttenCoef, size_t specIndex)
{
	massAttenCoef = _massAttenCoef[specIndex][0] / 10.f;
}

void PhantomObject::computeDesity(float* density)
{
	for (size_t i = 0; i < _pNumX * _pNumY * _pNumZ; i++)
	{
		density[i] = _structureMark[i] * _mPhantomMaterial.density[0];  // ������ʱ���ṹ���ֻ��0��1
	}
}



void PhantomObject::_readPhantomStructure()
{
	_structureMark = new uint8[_pNumX * _pNumY * _pNumZ];

	ifstream ifs;
	
	ifs.open(_mPhantomMaterial.structureMarkPath, ios::in | ios::binary);
	if (!ifs.is_open())
	{
		cout << "ģ��ṹ������ݴ�ʧ��" << endl;
		system("pause");
		exit(0);
	}

	ifs.read((char*)_structureMark, _pNumX * _pNumY * _pNumZ * sizeof(uint8));
	ifs.close();
}

void PhantomObject::_readMassAtten()
{
	cout << endl << "[����˥��ϵ��]" << endl;

	// _massAttenCoef[��������][���ʸ���], ��������������
	_massAttenCoef = new float* [_specEnergyNum];
	for (int i = 0; i < _specEnergyNum; i++)
	{
		_massAttenCoef[i] = new float[_mPhantomMaterial.materialNum];
	}

	ifstream ifs;
	int j;
	float buf;
	for (int i = 0; i < _mPhantomMaterial.materialNum; i++)
	{
		cout << "Material" << i + 1 << ": "<< _mPhantomMaterial.massAttenPath[i] << endl;
		ifs.open(_mPhantomMaterial.massAttenPath[i], ios::in);
		if (!ifs.is_open())
		{
			cout << "��" << i + 1 << "�����ʵ�����˥��ϵ�����ݴ�ʧ��!" << endl;
			system("pause");
			exit(0);
		}

		j = 0;
		while (ifs >> buf && j < _specEnergyNum)  // ���Ӻ����������Ϊ�˶�ȡ����Ҫ�����ݣ��뿪�ٿռ��С��ͬ����ֹ�ڴ�й¶
		{
			// ��ȡ��������˥��ϵ��
			_massAttenCoef[j][i] = buf;
			j++;
		}

		ifs.close();

		
	}
	cout << "======== ��ȡ�ɹ���========" << endl;
}