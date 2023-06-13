#pragma once

#include "CBCTPolyForwardProjGrid.h"

class CBCTSinMatPolyForwardProjGrid : private CBCTPolyForwardProjGrid
{
public:
	CBCTSinMatPolyForwardProjGrid(CTScanParas inCTScanParas, GridInfo inGridInfo, FilePath inFileInfo) : CBCTPolyForwardProjGrid(inCTScanParas, inGridInfo, inFileInfo) {};
	~CBCTSinMatPolyForwardProjGrid();

	virtual void computePolyForwProj();
	void computePolyForwProjNoResponse();
	void computePolyForwProjFoSp();

private:

	virtual void readPhantom();
	
	void readPhantomMassAtten();
	
};

