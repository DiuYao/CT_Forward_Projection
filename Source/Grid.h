#pragma once
#include <cstddef>
#include <iostream>
#include <string>

/*
* ����դ��̽����������Ϊ��������λ um
*/

/* ****************** V1 *********************
* 
******************************************* */

/* ****************** V2 *********************
* ����դ���ĺ�ȣ���������դ���Թ��ӵ�����
******************************************* */

/* ****************** V3 *********************
* ���Ǿ۽�������դդ���ĺ�ȱ任����������ԽԶ�����Խ�󣬽������ǲ�ͬλ�öԹ������յĲ�ͬ��
* ����Ĭ����ɢ��׼��1um��
* �����˵���ó��������Ӱ�����ӣ�
* ����դ�Ľ��ࡢ���ڡ����ϡ�դ����ȡ�դ�����롢դ�����ϣ���������������˥��ϵ��(���ϣ�����)
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
    size_t dNums;                   // ��λ: pix
    float dSize;                   // ��λ: um
    int leadStripsDistance;         // ��λ: um
    int leadStripsWidth;            // ��λ: um
    float h;                        // ��λ: mm
    float FD;                       // ����, ��λ: mm

    std::string materialGridStrip;            // դ������
    float rhoGS;                              // դ���ܶ�   g/cm^3 
    std::string materialGridInterspacer;      // ��������

    float uGridStrip;                        // դ�� Line attenuation coefficient  ��λ: 1/mm
    float uInterspacer;                      // ����� Line attenuation coefficient  ��λ: 1/mm

    size_t mI0Photons;              // I0�Ĺ�����

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
    // ����դ��Ϣ��û�б��浽���صĹ���
    void updataGrid(float uGridStrip);
    void showGridH();
    void showPeriodIncex();
    void getPeriodIncex(vector<size_t>& gridPeriodIndex);

    void creatBackground();
    void addGrid(char* fileNameI);

private:
    void computeunitH();

    // ��ɢ��Ԫ��ֵ, ��Ԫ����Even����λ��϶��Ԫ����Even, ��leadStripsDistance = Even um
    void createDiscreteGrid_uEdE();  

    // ��ɢ��Ԫ��ֵ, ��Ԫ����Even����λ��϶��Ԫ����Odd, ��leadStripsDistance = Odd um
    void createDiscreteGrid_uEdO();

    // ��ɢ��Ԫ��ֵ, ��Ԫ����Odd����λ��϶��Ԫ����Even, ��leadStripsDistance = Even um
    void createDiscreteGrid_uOdE();

    // ��ɢ��Ԫ��ֵ, ��Ԫ����Odd����λ��϶��Ԫ����Odd, ��leadStripsDistance = Odd um
    void createDiscreteGrid_uOdO();

    void discrateUnitBinning();
    void binning1D();


    // ��ӡģ������������դ��Ϣ
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

    int leadStripsDistance;     // um
    int leadStripsWidth;        // um
    int gridPeriod;             // um
    float h;                    // դ�� um
    float FD;                   // ���� um

   
    std::string materialGridStrip;         // դ���Ĳ���
    float uGridStrip;                      // ����˥��ϵ�� 1/um
    float rhoGS;                              // դ���ܶ�   g/cm^3 
    std::string materialGridInterspacer;   // �����Ĳ���
    float uInterspacer;                    // ����˥��ϵ�� 1/um

    size_t discreteUnitNums;    // ��ɢ��Ԫ����
    
    char* fileName;

    float* unitH;               // ��ɢ����ĵ�Ԫ���  um
    float* distanceCenter;      // ������դ���ĵľ���   

    float* grid1D;              // 1Dդ
    float* grid;                // ����������ĵ�ַ
    size_t* mI;                 // ���չ������Ƕ�����
    float* mI0;
    size_t mI0Photons;          // I0�Ĺ�����

    // �������������������
    vector<size_t> periodIndex;
};

