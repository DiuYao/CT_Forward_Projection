#pragma once

#include <vector>

#include "CBCTForwardProj.h"
#include "Grid.h"
#include ".\CUDAFun\CBCTPolyProjGPU.cuh"


class CBCTPolyForwardProjGrid : public CBCTForwardProj
{
public:
	CBCTPolyForwardProjGrid(CTScanParas inCTScanParas, GridInfo inGridInfo, FilePath inFileInfo);
	~CBCTPolyForwardProjGrid();

	virtual void computePolyForwProj();

//private:
	virtual void initGridInfo();

	virtual void saveIandIAbsorb();

	virtual void addI0();
	virtual void addI0NoResponse();
	virtual void saveI0();

	virtual void readPhantom();
	// 程序进度显示
	virtual void showProcessInfo();
	virtual void showProcessInfoNoResponse();
	// 卷积
	virtual void scatterSimulGrid();

	virtual void saveProj();

public:
//private:
	Grid* mGrid;
	GridAndDetectorSystem mGridAndDetectorSystem;
	vector<size_t> gridPeriodIndex;
};

