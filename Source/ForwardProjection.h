#pragma once

#include "DataType.h"
#include "CBCTPolyForwardProjGrid.h"
#include "CBCTPolyForwardProjNoGrid.h"
#include "CBCTSinMatPolyForwardProjGrid.h"
#include "CBCTSinMatPolyForwardProjNoGrid.h"

class ForwardProjection
{
public:
	void forwardPolyProjGrid();
	void forwardSinMatPolyProjGrid();
	void forwardSinMatPolyProjGridNoResponse();

	void forwardPolyProjNoGrid();
	void forwardSinMatPolyProjNoGrid();
	void forwardSinMatPolyProjNoGridNoResponse();

private:
	void inCTScanInfo();

	void inCTScanInfoGrid();
	void inCTScanInfoGridNoResponse();
	void inCTScanInfoSinMat();
	void inCTScanInfoNoGrid();
	void inCTScanInfoNoGridNoResponse();

private:

	CTScanParas mCTScanParas;
	GridInfo mGridInfo;   // ’§–≈œ¢
	FilePath mFilePath;


};

