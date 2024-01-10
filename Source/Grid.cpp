
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
        distanceCenter[i] = i + 0.5 - (float)discreteUnitNums / 2;  // 加0.5是因为离散化单位是1um
    }

    for (size_t i = 0; i < discreteUnitNums; i++)
    {
        unitH[i] = h / FD * sqrtf(FD * FD + distanceCenter[i] * distanceCenter[i]);   // 计算结果单位: um
    }
}

// 离散单元赋值, 单元总数Even，单位间隙单元总数Even, 即leadStripsDistance = Even um
void Grid::createDiscreteGrid_uEdE()
{
    // 中间位置索引 discreteUnitNums / 2 - 1, 对于一半栅，空隙部分的一半索引 [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 - 1]
    for (size_t i = discreteUnitNums / 2 - 1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]);
    }

    // 余下的半栅赋值, 起始是 discreteUnitNums / 2 - interspcaceWidth / 2 - 1
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // 透过系数
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

    // 对称补满
    for (size_t i = discreteUnitNums / 2; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i - 1];
    }
}

// 离散单元赋值, 单元总数Even，单位间隙单元总数Odd, 即leadStripsDistance = Odd um
void Grid::createDiscreteGrid_uEdO()
{
    // 中间位置索引 discreteUnitNums / 2 - 1, 中间间隙部分离散化单元是奇数，取中点 discreteUnitNums / 2 ，单独赋值
    // 除中点外，剩余部分的一半栅，空隙部分的一半索引 [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 - 1]， leadStripsDistance离散后单元数 >= 3
    grid1D[discreteUnitNums / 2] = expf(-interspaceLinearAtten * unitH[discreteUnitNums / 2]);
    for (size_t i = discreteUnitNums / 2 - 1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // 余下的半栅赋值, 起始是 discreteUnitNums / 2 - interspcaceWidth / 2 - 1
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // 透过系数
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
    // 对称补满，前一半比后一半多一个离散单元，（中间间隙离散单元奇数，中点占1 -- 即索引 discreteUnitNums / 2）
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i];
    }
}

// 离散单元赋值, 单元总数Odd，单位间隙单元总数Even, 即leadStripsDistance = Even um
void Grid::createDiscreteGrid_uOdE()
{
    // 中间位置索引取 discreteUnitNums / 2, (因为 int / 2，结果向下取整)
    // 对于一半栅，空隙部分的一半索引 [discreteUnitNums / 2 - interspcaceWidth / 2 + 1, discreteUnitNums / 2]
    for (size_t i = discreteUnitNums / 2; i > discreteUnitNums / 2 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // 余下的半栅赋值, 起始是 discreteUnitNums / 2 - interspcaceWidth / 2
    for (int i = discreteUnitNums / 2 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // 透过系数
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

    // 对称补满，前一半比后一半多一个离散单元，因为总单元数是奇数，间隙单元数是偶数，取总的中心索引为 discreteUnitNums / 2为间隙单元的前一半的结尾。
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i];
    }

}

// 离散单元赋值, 单元总数Odd，单位间隙单元总数Odd, 即leadStripsDistance = Odd um
void Grid::createDiscreteGrid_uOdO()
{
    // 中间位置索引取 discreteUnitNums / 2(中点), (因为 int / 2，结果向下取整)
    // 对于一半栅，空隙部分的一半索引 [discreteUnitNums / 2 - interspcaceWidth / 2, discreteUnitNums / 2 -1]， interspcaceWidth >= 3
    grid1D[discreteUnitNums / 2] = expf(-interspaceLinearAtten * unitH[discreteUnitNums / 2]);
    for (size_t i = discreteUnitNums / 2 -1; i > discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i--)
    {
        grid1D[i] = expf(-interspaceLinearAtten * unitH[i]); // expf(-gridStripLinearAtten * unitH[i]);
    }

    // 余下的半栅赋值, 起始是 discreteUnitNums / 2 - interspcaceWidth / 2
    for (int i = discreteUnitNums / 2 - 1 - interspcaceWidth / 2; i >= 0; i -= gridPeriod)
    {
        for (int stripIndex = 0; stripIndex < gridStripsWidth; stripIndex++)
        {
            if (i - stripIndex < 0)
            {
                break;
            }
            grid1D[i - stripIndex] = expf(-gridStripLinearAtten * unitH[i - stripIndex]);    // 透过系数
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

    // 对称补满，取总的中心索引为 discreteUnitNums / 2 做对称
    for (size_t i = discreteUnitNums / 2 + 1; i < discreteUnitNums; i++)
    {
        grid1D[i] = grid1D[discreteUnitNums - i - 1];
    }

    // 另一种方法布满
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


    // 由一排生成整个探测器的
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
    cout << "*******************成像系统信息******************" << endl;
    cout.width(15);
    cout << "探测器尺寸:" << dNums << "x" << dNums << endl;
    cout.width(15);
    cout << "探测器分辨率:"; cout.width(8); cout << dSize << " um" << endl;

    cout << endl;

    cout << "*****************聚焦型滤线栅信息****************" << endl;
    cout.width(15);
    cout << "栅条材质:" << gridStripMaterial << endl;
    cout.width(15);
    cout << "线性衰减系数:"; cout.width(10); cout << gridStripLinearAtten * 1000 << " 1/mm" << endl;
    cout.width(15);
    cout << "间隙物材质:" << interspaceMaterial << endl;
    cout.width(15);
    cout << "线性衰减系数:"; cout.width(10); cout << interspaceLinearAtten * 1000 << " 1/mm" << endl;
    cout.width(15);
    cout << "周期:"; cout.width(10); cout << gridPeriod << " um" << endl;
    cout.width(15);
    cout << "栅条宽度:"; cout.width(10); cout << gridStripsWidth << " um" << endl;
    cout.width(15);
    cout << "间隙距离:"; cout.width(10); cout << interspcaceWidth << " um" << endl;
    cout.width(15);
    cout << "厚度:"; cout.width(10); cout << h / 1000 << " mm" << endl;
    cout.width(15);
    cout << "焦距:"; cout.width(10); cout << FD / 1000 << " mm" << endl;
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
    // 离散化单元总数
    discreteUnitNums = dNums * dSize;   //dNums * dSize / 1
    
    unitH = new float[discreteUnitNums];
    distanceCenter = new float[discreteUnitNums];

    computeunitH();

    // 创建离散化后的一维栅
    grid1D = new float[discreteUnitNums];

    gridPeriod = interspcaceWidth + gridStripsWidth;

    

    // 合并像素，一排的像素, 并且copy为整个栅
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

// 输入栅条衰减系数的单位为 1/mm.
void Grid::updataGrid(float gridStripLinearAtten, float interspaceLinearAtten)
{
    this->gridStripLinearAtten = gridStripLinearAtten / 1000;   // 1/mm -> 1/um 
    this->interspaceLinearAtten = interspaceLinearAtten / 1000;   // 1/mm -> 1/um 

    // 离散化单元总数是偶数, 间隙距离偶数um
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
    
    // 离散化单元总数
    discreteUnitNums = dNums * dSize;   //dNums * dSize / 1

    unitH = new float[discreteUnitNums];
    distanceCenter = new float[discreteUnitNums];

    computeunitH();

    // 创建离散化后的一维栅
    grid1D = new float[discreteUnitNums];

    gridPeriod = interspcaceWidth + gridStripsWidth;

    // 离散化单元总数是偶数, 间隙距离偶数um
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

    // 合并像素，一排的像素, 并且copy为整个栅
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
        std::cout << "文件打开失败！" << std::endl;
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
