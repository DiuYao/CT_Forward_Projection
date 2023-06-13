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
	void forwardSinMatPolyProjGridFoSp();
	void forwardSinMatPolyProjGridNoResponse();

	void forwardPolyProjNoGrid();
	void forwardSinMatPolyProjNoGrid();
	void forwardSinMatPolyProjNoGridNoResponse();

	void forwarProjNoGridTestI0();

private:
	void inCTScanInfo();

	void inCTScanInfoGrid();
	void inCTScanInfoGridNoResponse();
	void inCTScanInfoSinMat();
	void inCTScanInfoNoGrid();
	void inCTScanInfoNoGridNoResponse();

	void inCTScanInfoNoGridTestI0();

private:

	CTScanParas mCTScanParas;
	GridInfo mGridInfo;   // ’§–≈œ¢
	FilePath mFilePath;


};

