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

typedef unsigned char uint8;
typedef unsigned int uint32;

// 闪烁体信息
struct ScintilltorInfo
{
	string scintillatorName;				// 闪烁体材质
	float scintillatorDensity;				// 闪烁体密度  单位 g/cm^3
	float scintillatorThickness;			// 闪烁体厚度  单位 mm
	float scintillatorThicknessErr;			// 闪烁体加工误差  单位 mm
	float detResponseFactor;				// 响应系数
};

// 扫描参数
struct CTScanParas
{
	size_t pNumX, pNumY, pNumZ;					// 模体像素数 单位 pix
	float pLengthX, pLengthY, pLengthZ;			// 模体大小 单位 mm

	size_t projNum;							// 扫描角度数
	size_t dNumU, dNumV;					// 探测器大小, dNumU-横向, dNumV--纵向
	float sod, sdd;							// 单位 mm
	int rotatedDirection;					// 扫描时的旋转方向
	
	float dSize;							 // 探测器分辨率   mm

	size_t I0Val;							// 程序中的光子初始数量，注意，在有探元响应时和I0不等
	
	float focalSpotSize;					// 焦斑大小, 单位 mm

	int spectrumStep;						// 能谱间隔
	int specEnergyNum;						// 离散化能谱个数

	ScintilltorInfo mScintilltorInfo;		// 闪烁体信息

	//GridAndDetectorSystem mGridAndDetectorSystem;
};

struct GridInfo
{
	// um
	int interspcaceWidth;         // 栅条间隔	单位: um
	int gridStripsWidth;            // 栅条宽度	单位: um
	float h;                        // 厚度		单位: mm
	float FD;                       // 焦距		单位: mm

	std::string gridStripMaterial;            // 栅条材料
	float gridStripDensity;                   // 栅条密度   g/cm^3 
	std::string gridStripMAPath;				// 栅条质量衰减系数路径

	std::string interspaceMaterial;			// 间隔物材料
	float interspaceDensity;
	std::string interspaceMAPath;				// 间隔物衰减系数路径

	float gridStripLinearAtten;                        // 栅条 Line attenuation coefficient  单位: 1/mm
	float interspaceLinearAtten;                      // 间隔物 Line attenuation coefficient  单位: 1/mm
};

// 系统参数
struct CTScanSystemInfo
{
	float totalTime = 0.0f;

	int rotatedDirection;

	// Compute parameters
	float thetaStep;
	float dHalfLU;
	float dHalfLV;
	float FOVR;		// 水平视野圆半径	单位 mm
	float FOVH;		// 垂直视野半高		单位 mm

	// 模体半长
	float pHalfX;	// 单位 mm
	float pHalfY;	// 单位 mm
	float pHalfZ;	// 单位 mm

	// Pixel size
	float pSizeX;	// X方向	单位 mm
	float pSizeY;	// Y方向	单位 mm
	float pSizeZ;	// Z方向	单位 mm

	// Integral calculus
	float dx;
	float dy;
	float dz;

	// Compute detector coordinates
	size_t intNum;        // 某条射线上的积分点个数

	float spectrumVal;		// 能谱值
	float phantomMAtten;	// 单材质时的质量衰减系数 单位 cm^2/g
};

// 焦斑，焦点移动
//struct FocalSpotOffset
//{
//	float* offsetU, * offsetV;
//};

// 坐标点
struct Coordinate
{
	float* imgIntX, * imgIntY, * imgIntZ;
	//FocalSpotOffset mFoSpotOffset;
	//float* detU, * detV;
};

// 
struct PolyForwardProj
{
	uint8* phantomStructure;
	float* phantom;							// 单材质时是密度(g/cm^3)，多材质时是线性衰减系数
	float* phantomMassAtten;				// 模体的质量衰减系数，单材质时使用此参数 单位 cm^2/g
	float* spectrumNormal;					// 归一化能谱
	float* grid;							// 栅数据
	float* gridLinearAtten;					// 栅条线性衰减系数
	float* interspaceLinearAtten;			// 间隔填充物线性衰减系数
	float* I, * I0, * IAbsorb;				// 过程中的穿透I，物体吸收IAbsorb，亮场I0
	float* proj;							// Proj

	float* scintillatorLineAtten;			// 闪烁体线性衰减系数  单位 cm^2/g
	float* scintillatorPerThickness;		// 每个闪烁体厚度			单位 mm
	float* detResponse;						// 探测器单元响应

	float* foSpOffsetU, * foSpOffsetV;		// 焦点偏移量
};
typedef PolyForwardProj ForwardProj;

// 路径
struct FilePath
{
	std::string spectrumPath;
	std::string gridPath;
	std::string phantomPath;				// 模体线性衰减系数图像路径,非完整路径, 具体以程序要求为主; 单材质时是密度图像路径,输入完整路径
	std::string phantomMAttenPath;			// 单材质时的模体质量衰减系数图像路径, 完整路径
	std::string outputFolder;
	std::string scintillatorPath;
	std::string convKernelGridPath;
};


struct PhantomMaterial
{
	int materialNum;
	std::string structureMarkPath;

	float* density;
	std::string* massAttenPath;
};



enum ScanMode
{
	TestI0,
	GridNoResponse,
	GridResponse,
	NoGridNoResponse,
	NoGridResponse
	
};



