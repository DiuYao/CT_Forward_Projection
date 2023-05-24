#pragma once

#include "CBCTForwardProj.h"

#include ".\CUDAFun\CBCTPolyProjGPU.cuh"
#include "CBCTForwardProj.h"


class CBCTPolyForwardProjNoGrid : public CBCTForwardProj
{
public:
	CBCTPolyForwardProjNoGrid(CTScanParas inCTScanParas, FilePath inFileInfo);
	~CBCTPolyForwardProjNoGrid();

	virtual void computePolyForwProj();

public:
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
};


