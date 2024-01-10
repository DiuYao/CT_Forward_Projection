#pragma once
#include <cstddef>
#include <iostream>
#include <string>

/*
* 滤线栅和探测器参数须为整数，单位 um
*/

/* ****************** V1 *********************
* 
******************************************* */

/* ****************** V2 *********************
* 考虑栅条的厚度，即考虑了栅条对光子的吸收
******************************************* */

/* ****************** V3 *********************
* 考虑聚焦型滤线栅栅条的厚度变换，即离中心越远，厚度越大，进而考虑不同位置对光子吸收的不同。
* 程序默认离散标准是1um。
* 具体的说，该程序加入了影响因子：
* 滤线栅的焦距、周期、材料、栅条宽度、栅条距离、栅条材料，光子能量，线性衰减系数(材料，能量)
******************************************* */


#include <vector>

#define DELETENEW(val) \
if(val != nullptr){\
delete[] val;\
val = nullptr;}


#define DELETE(var)\
if (var != nullptr){\
delete var;\
var = nullptr;\
}

using namespace std;

struct GridAndDetectorSystem
{
    // um
    size_t dNums;                   // 单位: pix
    float dSize;                   // 单位: um
    int interspcaceWidth;         // 单位: um
    int gridStripsWidth;            // 单位: um
    float h;                        // 单位: mm
    float FD;                       // 焦距, 单位: mm

    std::string gridStripMaterial;              // 栅条材料
    float gridStripDensity;                     // 栅条密度   g/cm^3 
    std::string gridStripMAPath;				    // 栅条质量衰减系数路径

    std::string interspaceMaterial;         // 间隔物材料
    float interspaceDensity;
    std::string interspaceMAPath;				// 间隔物衰减系数路径

    float gridStripLinearAtten;                        // 栅条 Line attenuation coefficient  单位: 1/mm
    float interspaceLinearAtten;                      // 间隔物 Line attenuation coefficient  单位: 1/mm

    size_t mI0Photons;              // I0的光子数

    size_t angleNums;
};

class Grid
{
public:
    Grid(GridAndDetectorSystem mGridAndDetectorSystem);
    Grid(float* gridOutput, GridAndDetectorSystem mGridAndDetectorSystem);
    ~Grid();

    
    void creatGrid();
    void defGrid();
    // 更新栅信息，没有保存到本地的功能
    void updataGrid(float gridStripLinearAtten, float interspaceLinearAtten);
    
    void showGridH();
    void showPeriodIncex();
    void getPeriodIncex(vector<size_t>& gridPeriodIndex);

    void creatBackground();
    void addGrid(char* fileNameI);

private:
    void computeunitH();

    // 离散单元赋值, 单元总数Even，单位间隙单元总数Even, 即leadStripsDistance = Even um
    void createDiscreteGrid_uEdE();  

    // 离散单元赋值, 单元总数Even，单位间隙单元总数Odd, 即leadStripsDistance = Odd um
    void createDiscreteGrid_uEdO();

    // 离散单元赋值, 单元总数Odd，单位间隙单元总数Even, 即leadStripsDistance = Even um
    void createDiscreteGrid_uOdE();

    // 离散单元赋值, 单元总数Odd，单位间隙单元总数Odd, 即leadStripsDistance = Odd um
    void createDiscreteGrid_uOdO();

    void discrateUnitBinning();
    void binning1D();


    // 打印模拟条件和滤线栅信息
    void printSimuInfo();


    int greatestCommonDivisor();
    void loadIntensity(char* fileName);

    void savemI0();
    void saveIntensity();


    void RetrievePeriodIndex();
    void RetrievePeriodIndex(vector<size_t>& gridPeriodIndex);
    void saveAsPeriodIndex();
    void saveAsPeriodIndex(size_t* gridPeriodIndex);

private:
    size_t dNums;
    size_t angleNums;
    size_t dSize;   // um

    int interspcaceWidth;     // um
    int gridStripsWidth;        // um
    int gridPeriod;             // um
    float h;                    // 栅厚 um
    float FD;                   // 焦距 um

   
    std::string gridStripMaterial;          // 栅条的材料
    float gridStripLinearAtten;             // 线性衰减系数 1/um
    float gridStripDensity;                 // 栅条密度   g/cm^3 
    
    std::string interspaceMaterial;         // 间隔物的材料
    float interspaceLinearAtten;           // 线性衰减系数 1/um
    float interspaceDensity;               // 栅条密度   g/cm^3 

    size_t discreteUnitNums;            // 离散单元总数
    
    char* fileName;

    float* unitH;               // 离散化后的单元厚度  um
    float* distanceCenter;      // 到滤线栅中心的距离   

    float* grid1D;              // 1D栅
    float* grid;                // 绑定输入参数的地址
    size_t* mI;                 // 接收光子数是短整型
    float* mI0;
    size_t mI0Photons;          // I0的光子数

    // 检索周期索引所需变量
    vector<size_t> periodIndex;
};

