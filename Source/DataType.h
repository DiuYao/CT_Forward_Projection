#pragma once

#include <string>

using namespace std;

#define DELETEARR(arrP)\
if (arrP != nullptr){\
delete[] arrP;\
arrP = nullptr;\
}

#define DELETE(var)\
if (var != nullptr){\
delete var;\
var = nullptr;\
}

typedef unsigned int uint32;

// ��˸����Ϣ
struct ScintilltorInfo
{
	string scintillatorName;				// ��˸�����
	float scintillatorDensity;				// ��˸���ܶ�  ��λ g/cm^3
	float scintillatorThickness;			// ��˸����  ��λ mm
	float scintillatorThicknessErr;			// ��˸��ӹ����  ��λ mm
	float detResponseFactor;				// ��Ӧϵ��
};

// ɨ�����
struct CTScanParas
{
	size_t pNumX, pNumY, pNumZ;					// ģ�������� ��λ pix
	float pLengthX, pLengthY, pLengthZ;			// ģ���С ��λ mm

	size_t projNum;							// ɨ��Ƕ���
	size_t dNumU, dNumV;					// ̽������С, dNumU-����, dNumV--����
	float sod, sdd;							// ��λ mm
	int rotatedDirection;					// ɨ��ʱ����ת����

	float dSize;							 // ̽�����ֱ���   mm

	uint32 I0Val;							// ����
	int spectrumStep;						// ���׼��
	int specEnergyNum;						// ��ɢ�����׸���

	ScintilltorInfo mScintilltorInfo;		// ��˸����Ϣ

	//GridAndDetectorSystem mGridAndDetectorSystem;
};

struct GridInfo
{
	// um
	int leadStripsDistance;         // դ�����	��λ: um
	int leadStripsWidth;            // դ�����	��λ: um
	float h;                        // ���		��λ: mm
	float FD;                       // ����		��λ: mm

	std::string materialGridStrip;            // դ������
	float rhoGS;                              // դ���ܶ�   g/cm^3 
	std::string materialGridInterspacer;      // ��������

	float uGridStrip;                        // դ�� Line attenuation coefficient  ��λ: 1/mm
	float uInterspacer;                      // ����� Line attenuation coefficient  ��λ: 1/mm
};

// ϵͳ����
struct CTScanSystemInfo
{
	float totalTime = 0.0f;

	int rotatedDirection;

	// Compute parameters
	float thetaStep;
	float dHalfLU;
	float dHalfLV;
	float FOVR;		// ˮƽ��ҰԲ�뾶	��λ mm
	float FOVH;		// ��ֱ��Ұ���		��λ mm

	// ģ��볤
	float pHalfX;	// ��λ mm
	float pHalfY;	// ��λ mm
	float pHalfZ;	// ��λ mm

	// Pixel size
	float pSizeX;	// X����	��λ mm
	float pSizeY;	// Y����	��λ mm
	float pSizeZ;	// Z����	��λ mm

	// Integral calculus
	float dx;
	float dy;
	float dz;

	// Compute detector coordinates
	size_t intNum;        // ĳ�������ϵĻ��ֵ����

	float spectrumVal;		// ����ֵ
	float phantomMAtten;	// ������ʱ������˥��ϵ�� ��λ cm^2/g
};

// �����
struct Coordinate
{
	float* imgIntX, * imgIntY, * imgIntZ;
	//float* detU, * detV;
};

// 
struct PolyForwardProj
{
	float* phantom;							// ������ʱ���ܶ�(g/cm^3)�������ʱ������˥��ϵ��
	float* phantomMassAtten;				// ģ�������˥��ϵ����������ʱʹ�ô˲��� ��λ cm^2/g
	float* spectrumNormal;					// ��һ������
	float* grid;							// դ����
	float* gridLineAtten;					// դ������˥��ϵ��
	float* I, * I0, * IAbsorb;				// �����еĴ�͸I����������IAbsorb������I0
	float* proj;							// Proj

	float* scintillatorLineAtten;			// ��˸������˥��ϵ��  ��λ cm^2/g
	float* scintillatorPerThickness;		// ÿ����˸����			��λ mm
	float* detResponse;						// ̽������Ԫ��Ӧ
};
typedef PolyForwardProj ForwardProj;

// ·��
struct FilePath
{
	std::string spectrumPath;
	std::string gridPath;
	std::string phantomPath;				// ģ������˥��ϵ��ͼ��·��,������·��, �����Գ���Ҫ��Ϊ��; ������ʱ���ܶ�ͼ��·��,��������·��
	std::string phantomMAttenPath;			// ������ʱ��ģ������˥��ϵ��ͼ��·��, ����·��
	std::string IEPath;
	std::string scintillatorPath;
	std::string convKernelGridPath;
};










