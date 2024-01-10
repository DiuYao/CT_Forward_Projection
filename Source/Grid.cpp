
#include "Grid.h"

#include <fstream>
#include <iostream>



Grid::Grid(GridAndDetectorSystem mGridAndDetectorSystem)
    : dNums(mGridAndDetectorSystem.dNums)
    , dSize(mGridAndDetectorSystem.dSize)
    , interspcaceWidth(mGridAndDetectorSystem.interspcaceWidth)
    , gridStripsWidth(mGridAndDetectorSystem.gridStripsWidth)
    , angleNums(mGridAndDetectorSystem.angleNums)
    , mI0Photons(mGridAndDetectorSystem.mI0Photons)
    , h(mGridAndDetectorSystem.h * 1000)    // mm --> um
    , FD(mGridAndDetectorSystem.FD *1000)   // mm --> um
    , gridStripMaterial(mGridAndDetectorSystem.gridStripMaterial)
    , gridStripLinearAtten(mGridAndDetectorSystem.gridStripLinearAtten / 1000)    // 1/mm --> 1/um
    , interspaceMaterial(mGridAndDetectorSystem.interspaceMaterial)
    , interspaceLinearAtten(mGridAndDetectorSystem.interspaceLinearAtten / 1000) // 1/mm --> 1/um  
    
{
    gridPeriod = 0;
    discreteUnitNums = 0;
    fileName = NULL;
    unitH = nullptr;
    distanceCenter = nullptr;
    grid1D = nullptr;
    grid = nullptr;
    mI0 = nullptr;
    mI = nullptr;

}

Grid::Grid(float* gridOutput, GridAndDetectorSystem mGridAndDetectorSystem)
    : dNums(mGridAndDetectorSystem.dNums)
    , dSize(mGridAndDetectorSystem.dSize)
    , interspcaceWidth(mGridAndDetectorSystem.interspcaceWidth)
    , gridStripsWidth(mGridAndDetectorSystem.gridStripsWidth)
    , angleNums(mGridAndDetectorSystem.angleNums)
    , mI0Photons(mGridAndDetectorSystem.mI0Photons)
    , h(mGridAndDetectorSystem.h * 1000)    // mm --> um
    , FD(mGridAndDetectorSystem.FD * 1000)   // mm --> um
    , gridStripMaterial(mGridAndDetectorSystem.gridStripMaterial)
    , interspaceMaterial(mGridAndDetectorSystem.interspaceMaterial)
    , gridStripLinearAtten(mGridAndDetectorSystem.gridStripLinearAtten / 1000)    // 1/mm --> 1/um
    , interspaceLinearAtten(mGridAndDetectorSystem.interspaceLinearAtten / 1000) // 1/mm --> 1/um 
    , gridStripDensity(mGridAndDetectorSystem.gridStripDensity)
    , interspaceDensity(mGridAndDetectorSystem.interspaceDensity)
    , grid(gridOutput)
{

    gridPeriod = 0;
    discreteUnitNums = 0;
    
    fileName = NULL;
    unitH = nullptr;
    distanceCenter = nullptr;
    grid1D = nullptr;
    mI0 = nullptr;
    mI = nullptr;
}

Grid::~Grid()
{
    DELETENEW(grid1D);
    DELETENEW(distanceCenter);
    DELETENEW(unitH);
    
    // DELETE(grid);
}



void Grid::computeunitH()
{
    for (size_t i = 0; i < discreteUnitNums; i++)
    {
        distanceCenter[i] = i + 0.5 - (float)discreteUnitNums / 2;  // ��0.5����Ϊ��ɢ����λ��1um
    }

    for (size_t i = 0; i < discreteUnitNums; i++)
    {
        unitH[i] = h / FD * sqrtf(FD * FD + distanceCenter[i] * distanceCenter[i]);   // ��������λ: um
    }
}

// ��ɢ��Ԫ��ֵ, ��Ԫ����Even����λ��϶��Ԫ����Even, ��leadStripsDistance = Even um
void Grid::createDiscreteGrid_uEdE()
{
    // �м�λ������ discreteUnitNums / 2 - 1, ����һ��դ����϶���ֵ�һ������ [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 - 1]
    for (size_t i = discreteUnitNums / 2 - 1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]);
    }

    // ���µİ�դ��ֵ, ��ʼ�� discreteUnitNums / 2 - interspcaceWidth / 2 - 1
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // ͸��ϵ��
        }
        for (int interspacerIndex = 0; interspacerIndex < interspcaceWidth; interspacerIndex++)
        {
            if (i - gridStripsWidth - interspacerIndex < 0)
            {
                break;
            }
            grid1D[i - gridStripsWidth - interspacerIndex] = expf(-interspaceLinearAtten * unitH[i]);
        }

        if (i - gridPeriod < 0)
        {
            //break;
        }

    }

    // �ԳƲ���
    for (size_t i = discreteUnitNums / 2; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i - 1];
    }
}

// ��ɢ��Ԫ��ֵ, ��Ԫ����Even����λ��϶��Ԫ����Odd, ��leadStripsDistance = Odd um
void Grid::createDiscreteGrid_uEdO()
{
    // �м�λ������ discreteUnitNums / 2 - 1, �м��϶������ɢ����Ԫ��������ȡ�е� discreteUnitNums / 2 ��������ֵ
    // ���е��⣬ʣ�ಿ�ֵ�һ��դ����϶���ֵ�һ������ [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 - 1]�� leadStripsDistance��ɢ��Ԫ�� >= 3
    grid1D[discreteUnitNums / 2] = expf(-interspaceLinearAtten * unitH[discreteUnitNums / 2]);
    for (size_t i = discreteUnitNums / 2 - 1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // ���µİ�դ��ֵ, ��ʼ�� discreteUnitNums / 2 - interspcaceWidth / 2 - 1
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // ͸��ϵ��
        }
        for (int interspacerIndex = 0; interspacerIndex < interspcaceWidth; interspacerIndex++)
        {
            if (i - gridStripsWidth - interspacerIndex < 0)
            {
                break;
            }
            grid1D[i - gridStripsWidth - interspacerIndex] = expf(-interspaceLinearAtten * unitH[i]);
        }

        if (i - gridPeriod < 0)
        {
            break;
        }

    }
    // �ԳƲ�����ǰһ��Ⱥ�һ���һ����ɢ��Ԫ�����м��϶��ɢ��Ԫ�������е�ռ1 -- ������ discreteUnitNums / 2��
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i];
    }
}

// ��ɢ��Ԫ��ֵ, ��Ԫ����Odd����λ��϶��Ԫ����Even, ��leadStripsDistance = Even um
void Grid::createDiscreteGrid_uOdE()
{
    // �м�λ������ȡ discreteUnitNums / 2, (��Ϊ int / 2���������ȡ��)
    // ����һ��դ����϶���ֵ�һ������ [discreteUnitNums / 2 - interspcaceWidth / 2 + 1, discreteUnitNums / 2]
    for (size_t i = discreteUnitNums / 2; i > discreteUnitNums / 2 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // ���µİ�դ��ֵ, ��ʼ�� discreteUnitNums / 2 - interspcaceWidth / 2
    for (int i = discreteUnitNums / 2 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // ͸��ϵ��
        }
        for (int interspacerIndex = 0; interspacerIndex < interspcaceWidth; interspacerIndex++)
        {
            if (i - gridStripsWidth - interspacerIndex < 0)
            {
                break;
            }
            grid1D[i - gridStripsWidth - interspacerIndex] = expf(-interspaceLinearAtten * unitH[i]);
        }

        if (i - gridPeriod < 0)
        {
            break;
        }

    }

    // �ԳƲ�����ǰһ��Ⱥ�һ���һ����ɢ��Ԫ����Ϊ�ܵ�Ԫ������������϶��Ԫ����ż����ȡ�ܵ���������Ϊ discreteUnitNums / 2Ϊ��϶��Ԫ��ǰһ��Ľ�β��
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i];
    }

}

// ��ɢ��Ԫ��ֵ, ��Ԫ����Odd����λ��϶��Ԫ����Odd, ��leadStripsDistance = Odd um
void Grid::createDiscreteGrid_uOdO()
{
    // �м�λ������ȡ discreteUnitNums / 2(�е�), (��Ϊ int / 2���������ȡ��)
    // ����һ��դ����϶���ֵ�һ������ [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 -1]�� interspcaceWidth >= 3
    grid1D[discreteUnitNums / 2] = expf(-interspaceLinearAtten * unitH[discreteUnitNums / 2]);
    for (size_t i = discreteUnitNums / 2 -1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // ���µİ�դ��ֵ, ��ʼ�� discreteUnitNums / 2 - interspcaceWidth / 2
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // ͸��ϵ��
        }
        for (int interspacerIndex = 0; interspacerIndex < interspcaceWidth; interspacerIndex++)
        {
            if (i - gridStripsWidth - interspacerIndex < 0)
            {
                break;
            }
            grid1D[i - gridStripsWidth - interspacerIndex] = expf(-interspaceLinearAtten * unitH[i]);
        }

        if (i - gridPeriod < 0)
        {
            break;
        }

    }

    // �ԳƲ�����ȡ�ܵ���������Ϊ discreteUnitNums / 2 ���Գ�
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i - 1];
    }

    // ��һ�ַ�������
    /*for (size_t i = 0; i < discreteUnitNums / 2 + 1; i++)
    {
        grid1D[discreteUnitNums / 2 + i] = grid1D[discreteUnitNums / 2 - i];
    }*/
}


void Grid::discrateUnitBinning()
{
    float sum = 0.0f;
    for (size_t i = 0; i < discreteUnitNums; i += dSize)
    {
        sum = 0.0f;
        for (size_t j = 0; j < dSize; j++)
        {
            sum += grid1D[i + j];
        }
        grid[i / dSize] = sum / dSize;
    }


    // ��һ����������̽������
    for (size_t i = 1; i < dNums; i++)
    {
        for (size_t j = 0; j < dNums; j++)
        {
            grid[i * dNums + j] = grid[j];
        }
    }
}

void Grid::binning1D()
{
    double sum = 0.0f;
    for (size_t i = 0; i < discreteUnitNums; i += dSize)
    {
        sum = 0.0f;
        for (size_t j = 0; j < dSize; j++)
        {
            sum += grid1D[i + j];
        }
        grid[i / dSize] = sum / dSize;
    }
}

void Grid::printSimuInfo()
{
    cout.setf(ios::left);
    cout << "*******************����ϵͳ��Ϣ******************" << endl;
    cout.width(15);
    cout << "̽�����ߴ�:" << dNums << "x" << dNums << endl;
    cout.width(15);
    cout << "̽�����ֱ���:"; cout.width(8); cout << dSize << " um" << endl;

    cout << endl;

    cout << "*****************�۽�������դ��Ϣ****************" << endl;
    cout.width(15);
    cout << "դ������:" << gridStripMaterial << endl;
    cout.width(15);
    cout << "����˥��ϵ��:"; cout.width(10); cout << gridStripLinearAtten * 1000 << " 1/mm" << endl;
    cout.width(15);
    cout << "��϶�����:" << interspaceMaterial << endl;
    cout.width(15);
    cout << "����˥��ϵ��:"; cout.width(10); cout << interspaceLinearAtten * 1000 << " 1/mm" << endl;
    cout.width(15);
    cout << "����:"; cout.width(10); cout << gridPeriod << " um" << endl;
    cout.width(15);
    cout << "դ�����:"; cout.width(10); cout << gridStripsWidth << " um" << endl;
    cout.width(15);
    cout << "��϶����:"; cout.width(10); cout << interspcaceWidth << " um" << endl;
    cout.width(15);
    cout << "���:"; cout.width(10); cout << h / 1000 << " mm" << endl;
    cout.width(15);
    cout << "����:"; cout.width(10); cout << FD / 1000 << " mm" << endl;
    cout << "*************************************************" << endl;

}

int Grid::greatestCommonDivisor()
{
	int x = interspcaceWidth;
	int y = gridStripsWidth;
	int z = dSize;

	int gcdSub;
	int gcd;

	gcdSub = y;
	while (x % y != 0)
	{
		gcdSub = x % y;
		x = y;
		y = gcdSub;
	}

	gcd = gcdSub;
	while (z % gcdSub != 0)
	{
		gcd = z % gcdSub;
		z = gcdSub;
        gcdSub = gcd;
	}

	return gcd;
}

void Grid::defGrid()
{
    // ��ɢ����Ԫ����
    discreteUnitNums = dNums * dSize;   //dNums * dSize / 1
    
    unitH = new float[discreteUnitNums];
    distanceCenter = new float[discreteUnitNums];

    computeunitH();

    // ������ɢ�����һάդ
    grid1D = new float[discreteUnitNums];

    gridPeriod = interspcaceWidth + gridStripsWidth;

    

    // �ϲ����أ�һ�ŵ�����, ����copyΪ����դ
    /*grid = new float[dNums * dNums];
    discrateUnitBinning();*/
    
    

    /*char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "GridImageData/Grid_dNums%dx%d_dSize%dum_S%s_IS%s_gridD%dum_gridd%dum_h%.2fmm_float_FD%.1fmm.raw", dNums, dNums, dSize, gridStripMaterial, interspaceMaterial, interspcaceWidth, gridStripsWidth, h/1000, FD/1000);

    ofstream ofs;
    ofs.open(pathWrite, ios::out | ios::binary);
    ofs.write((const char*)grid, dNums * dNums * sizeof(float));
    ofs.close();*/

    //printSimuInfo();
}

// ����դ��˥��ϵ���ĵ�λΪ 1/mm.
void Grid::updataGrid(float gridStripLinearAtten, float interspaceLinearAtten)
{
    this->gridStripLinearAtten = gridStripLinearAtten / 1000;   // 1/mm -> 1/um 
    this->interspaceLinearAtten = interspaceLinearAtten / 1000;   // 1/mm -> 1/um 

    // ��ɢ����Ԫ������ż��, ��϶����ż��um
    int UnitNumsFlag = discreteUnitNums % 2;
    int interSpaceFlag = interspcaceWidth % 2;
    switch (UnitNumsFlag)
    {
    case 0:
        if (interSpaceFlag == 0)
        {
            createDiscreteGrid_uEdE();
        }
        else if (interSpaceFlag == 1)
        {
            createDiscreteGrid_uEdO();
        }
        break;
    case 1:
        if (interSpaceFlag == 0)
        {
            createDiscreteGrid_uOdE();
        }
        else if (interSpaceFlag == 1)
        {
            createDiscreteGrid_uOdO();
        }
        break;
    default:
        break;
    }

    binning1D();
}

void Grid::creatGrid()
{
    
    // ��ɢ����Ԫ����
    discreteUnitNums = dNums * dSize;   //dNums * dSize / 1

    unitH = new float[discreteUnitNums];
    distanceCenter = new float[discreteUnitNums];

    computeunitH();

    // ������ɢ�����һάդ
    grid1D = new float[discreteUnitNums];

    gridPeriod = interspcaceWidth + gridStripsWidth;

    // ��ɢ����Ԫ������ż��, ��϶����ż��um
    int UnitNumsFlag = discreteUnitNums % 2;
    int interSpaceFlag = interspcaceWidth % 2;
    switch (UnitNumsFlag)
    {
    case 0:
        if (interSpaceFlag == 0)
        {
            createDiscreteGrid_uEdE();
        }
        else if (interSpaceFlag == 1)
        {
            createDiscreteGrid_uEdO();
        }
        break;
    case 1:
        if (interSpaceFlag == 0)
        {
            createDiscreteGrid_uEdE();
        }
        else if (interSpaceFlag == 1)
        {
            createDiscreteGrid_uEdO();
        }
        break;
    default:
        break;
    }

    // �ϲ����أ�һ�ŵ�����, ����copyΪ����դ
    grid = new float[dNums * dNums];
    discrateUnitBinning();


    char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "GridImageData/Grid_dNums%dx%d_dSize%dum_S%s_IS%s_gridD%dum_gridd%dum_h%.2fmm_float_FD%.1fmm.raw", dNums, dNums, dSize, gridStripMaterial, interspaceMaterial, interspcaceWidth, gridStripsWidth, h/1000, FD/1000);

    ofstream ofs;
    ofs.open(pathWrite, ios::out | ios::binary);
    ofs.write((const char*)grid, dNums * dNums * sizeof(float));
    ofs.close();

    printSimuInfo();
}

void Grid::showGridH()
{
    /*float* recordH = new float[discreteUnitNums];
    for (int i = 0; i < discreteUnitNums; i++)
    {
        recordH[i] = unitH[i];
    }*/

    char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "GridImageData/Grid_dNums%dx%d_dSize%dum_S%s_IS%s_gridD%dum_gridd%dum_h%.2fmm_float_FD%.1fmm_Hight%d.raw", dNums, dNums, dSize, gridStripMaterial, interspaceMaterial, interspcaceWidth, gridStripsWidth, h / 1000, FD / 1000, discreteUnitNums);

    ofstream ofs;
    ofs.open(pathWrite, ios::out | ios::binary);
    ofs.write((const char*)unitH, discreteUnitNums * sizeof(float));
    ofs.close();


}

void Grid::showPeriodIncex()
{
    RetrievePeriodIndex();
    saveAsPeriodIndex();
}

void Grid::getPeriodIncex(vector<size_t>& gridPeriodIndex)
{
    RetrievePeriodIndex();
    saveAsPeriodIndex();
    gridPeriodIndex = periodIndex;
}

void Grid::addGrid(char* fileNameI)
{
    loadIntensity(fileNameI);

    for (size_t angle = 0; angle < angleNums; angle++)
    {
        for (size_t i = 0; i < dNums * dNums; i++)
        {
            mI[angle * dNums * dNums + i] = mI[angle * dNums * dNums + i] * grid[i];
        }
    }

    saveIntensity();
}

void Grid::creatBackground()
{
    mI0 = new float[dNums * dNums];
    for (size_t i = 0; i < dNums * dNums; i++)
    {
        mI0[i] = grid[i] * mI0Photons;
    }
    savemI0();
}

void Grid::loadIntensity(char* mFileName)
{
    fileName = mFileName;

    mI = new size_t[angleNums * dNums * dNums];

    char pathRead[128];
    sprintf_s(pathRead, sizeof(pathRead), "IntensityImageData/%s.raw", fileName);

    ifstream ifs;
    ifs.open(pathRead, ios::in | ios::binary);
    if (!ifs.is_open())
    {
        std::cout << "�ļ���ʧ�ܣ�" << std::endl;
        exit(0);
    }
    ifs.read((char*)mI, angleNums * dNums * dNums * sizeof(unsigned short));
    ifs.close();
}

void Grid::savemI0()
{
    char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "AddIntensityImageData/Background(Grid)_dNums%dx%d_dSize%dum_gridD%dum_gridd%dum_h%.2f_Photon%d.raw", dNums, dNums, dSize, interspcaceWidth, gridStripsWidth, h, mI0Photons);

    ofstream ofs;
    ofs.open(pathWrite, ios::out | ios::binary);
    ofs.write((const char*)mI0, dNums * dNums * sizeof(float));
    ofs.close();
}

void Grid::saveIntensity()
{
    char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "AddIntensityImageData/AddGrid_dSize%dum_gridD%dum_gridd%dum_%s.raw", dSize, interspcaceWidth, gridStripsWidth, fileName);

    ofstream ofs;
    ofs.open(pathWrite, ios::out | ios::binary);
    ofs.write((const char*)mI, angleNums * dNums * dNums * sizeof(unsigned short));
    ofs.close();
}

void Grid::RetrievePeriodIndex()
{
    
    int i = 0;
    periodIndex.push_back(i);
    i++;
    while (i < dNums - 1)
    {
        if (grid[i - 1] > grid[i] && grid[i] < grid[i + 1])
        {
            periodIndex.push_back(i);
        }
        i++;
    }
    periodIndex.push_back(dNums - 1);
}

void Grid::RetrievePeriodIndex(vector<size_t>& gridPeriodIndex)
{
    int i = 0;
    gridPeriodIndex.push_back(i);
    i++;
    while (i < dNums - 1)
    {
        if (grid[i - 1] > grid[i] && grid[i] < grid[i + 1])
        {
            gridPeriodIndex.push_back(i);
        }
        i++;
    }
    gridPeriodIndex.push_back(dNums - 1);
}


void Grid::saveAsPeriodIndex()
{
    char pathWrite[256];
    sprintf_s(pathWrite, sizeof(pathWrite), "OutputResult/ForwardProjection/Grid/Grid_dNums%dx%d_dSize%dum_S%s_IS%s_gridD%dum_gridd%dum_h%.2fmm_FD%.1fmm_PeriodIndex.txt", dNums, dNums, dSize, gridStripMaterial.c_str(), interspaceMaterial.c_str(), interspcaceWidth, gridStripsWidth, h / 1000, FD / 1000);

    ofstream ofs;
    ofs.open(pathWrite, ios::out);


    for (vector<size_t> ::iterator it = periodIndex.begin(); it != periodIndex.end(); it++)
    {
        ofs << *it << endl;
    }

    //for (size_t i = 0; i < periodIndex.size(); i++)
    //{
    //    ofs << periodIndex.at(i) << endl;
    //    // ofs << periodIndex[i] << endl;
    //}
}

void Grid::saveAsPeriodIndex(size_t* gridPeriodIndex)
{
    gridPeriodIndex = new size_t[periodIndex.size()];

    char pathWrite[128];
    sprintf_s(pathWrite, sizeof(pathWrite), "OutputResult/Grid/Grid_dNums%dx%d_dSize%dum_S%s_IS%s_gridD%dum_gridd%dum_h%.2fmm_FD%.1fmm_PeriodIndex.txt", dNums, dNums, dSize, gridStripMaterial.c_str(), interspaceMaterial.c_str(), interspcaceWidth, gridStripsWidth, h / 1000, FD / 1000);

    ofstream ofs;
    ofs.open(pathWrite, ios::out);


	for (size_t i = 0; i < periodIndex.size(); i++)
	{
		ofs << periodIndex.at(i) << endl;
		// ofs << periodIndex[i] << endl;
	}
}
