#pragma once

#include "CBCTPolyForwardProjGrid.h"

class CBCTSinMatPolyForwardProjGrid : private CBCTPolyForwardProjGrid
{
public:
	CBCTSinMatPolyForwardProjGrid(CTScanParas inCTScanParas, GridInfo inGridInfo, FilePath inFileInfo, const PhantomMaterial inphantomMaterial) : CBCTPolyForwardProjGrid(inCTScanParas, inGridInfo, inFileInfo, inphantomMaterial) {};
	~CBCTSinMatPolyForwardProjGrid();

	virtual void computePolyForwProj();
	void computePolyForwProjNoResponse();
	void computePolyForwProjFoSp();

private:

	virtual void readPhantom();
	
	void readPhantomMassAtten();
	
};

