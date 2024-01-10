#pragma once
#include "CBCTPolyForwardProjNoGrid.h"
class CBCTSinMatPolyForwardProjNoGrid : private CBCTPolyForwardProjNoGrid
{
public:
	CBCTSinMatPolyForwardProjNoGrid(CTScanParas inCTScanParas, FilePath inFileInfo, const PhantomMaterial inPhantomMaterial) : CBCTPolyForwardProjNoGrid(inCTScanParas, inFileInfo, inPhantomMaterial) {};
	~CBCTSinMatPolyForwardProjNoGrid();

	virtual void computePolyForwProj();
	void computePolyForwProjNoResponse();
	void computePolyForwProjFoSp();

private:
	virtual void readPhantom();

	void readPhantomMassAtten();
};

