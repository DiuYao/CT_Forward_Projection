#pragma once

#include "DataType.h"
#include "CBCTPolyForwardProjGrid.h"
#include "CBCTPolyForwardProjNoGrid.h"
#include "CBCTSinMatPolyForwardProjGrid.h"
#include "CBCTSinMatPolyForwardProjNoGrid.h"
#include "InitFile.h"


class ForwardProjection
{
public:

	~ForwardProjection();

	void run();

private:
	void forwardPolyProjGrid();
	void forwardSinMatPolyProjGrid();
	void forwardSinMatPolyProjGridFoSp();
	void forwardSinMatPolyProjGridNoResponse();

	void forwardPolyProjNoGrid();
	void forwardSinMatPolyProjNoGrid();
	void forwardSinMatPolyProjNoGridFoSp();
	void forwardSinMatPolyProjNoGridNoResponse();

	void forwarProjNoGridTestI0();

private:
	void _getCTScanConfig();

	void _getCTScanGridConfig();
	void _getCTScanDetResponseConfig();

	void _getOutputFolderConfig();

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

	CInitFile mCInitFile;

	int _isTestI0;
	//int _isPolyenergetic;
	int _isDetResponse;
	int _isGrid;
	int _isSingleMaterial;

	PhantomMaterial _mPhantomMaterial;
};

