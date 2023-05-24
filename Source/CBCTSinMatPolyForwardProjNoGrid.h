#pragma once
#include "CBCTPolyForwardProjNoGrid.h"
class CBCTSinMatPolyForwardProjNoGrid : private CBCTPolyForwardProjNoGrid
{
public:
	CBCTSinMatPolyForwardProjNoGrid(CTScanParas inCTScanParas, GridInfo inGridInfo, FilePath inFileInfo) : CBCTPolyForwardProjNoGrid(inCTScanParas, inFileInfo) {};
	~CBCTSinMatPolyForwardProjNoGrid();

	virtual void computePolyForwProj();
	void computePolyForwProjNoResponse();

private:
	virtual void readPhantom();

	void readPhantomMassAtten();
};

