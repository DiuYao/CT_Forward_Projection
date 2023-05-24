#include "ForwardProjection.h"



#define PNUM 256
#define DNUM 600
#define PHANTOMLEN 51.2   // mm


void ForwardProjection::forwardPolyProjGrid()
{
    inCTScanInfo();
    inCTScanInfoGrid();

    CBCTPolyForwardProjGrid mCBCTPolyForwardProjGrid = CBCTPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTPolyForwardProjGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjGrid()
{
    inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjGridNoResponse()
{
    inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();
    inCTScanInfoGridNoResponse();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProjNoResponse();
}

void ForwardProjection::forwardPolyProjNoGrid()
{
    inCTScanInfo();
    inCTScanInfoNoGrid();

    CBCTPolyForwardProjNoGrid mCBCTPolyForwardProjNoGrid = CBCTPolyForwardProjNoGrid(mCTScanParas, mFilePath);
    mCBCTPolyForwardProjNoGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjNoGrid()
{
    inCTScanInfo();
    inCTScanInfoNoGrid();
    inCTScanInfoSinMat();

    CBCTSinMatPolyForwardProjNoGrid mCBCTSinMatPolyForwardProjNoGrid = CBCTSinMatPolyForwardProjNoGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProj();

}

void ForwardProjection::forwardSinMatPolyProjNoGridNoResponse()
{
    inCTScanInfo();
    inCTScanInfoNoGrid();
    inCTScanInfoSinMat();
    inCTScanInfoNoGridNoResponse();

    CBCTSinMatPolyForwardProjNoGrid mCBCTSinMatPolyForwardProjNoGrid = CBCTSinMatPolyForwardProjNoGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProj();

}

void ForwardProjection::inCTScanInfo()
{
    // ɨ��ϵͳ����
    mCTScanParas.pNumX = PNUM;
    mCTScanParas.pNumY = PNUM;
    mCTScanParas.pNumZ = PNUM;
	mCTScanParas.pLengthX = PHANTOMLEN;
	mCTScanParas.pLengthY = PHANTOMLEN;
	mCTScanParas.pLengthZ = PHANTOMLEN;
    mCTScanParas.dNumU = DNUM;
    mCTScanParas.dNumV = DNUM;
    mCTScanParas.projNum = 360;
    mCTScanParas.rotatedDirection = 1;      // -1 ���� ˳ʱ��, 1 ���� ��ʱ��

    mCTScanParas.dSize = 0.278;  //mm
    mCTScanParas.sdd = 1050;
    mCTScanParas.sod = 630;
    mCTScanParas.I0Val = 1 * 1e6;

    mCTScanParas.specEnergyNum = 80;
    mCTScanParas.spectrumStep = 1;

    // ��˸����Ϣ
    mCTScanParas.mScintilltorInfo.scintillatorDensity = 4.510;      // �ܶ�   g/cm^3
    mCTScanParas.mScintilltorInfo.scintillatorThickness = 0.4;      // ���   mm
    mCTScanParas.mScintilltorInfo.scintillatorThicknessErr = 0.05;   // ������� mm
    mCTScanParas.mScintilltorInfo.detResponseFactor = 1.0f;         // ��Ӧ����

    // ·��
    mFilePath.spectrumPath = "InputData/Spectrum/Spectrum_80kVp_1mmAl_GEMaxiray_1keV_Normal.txt";
    mFilePath.phantomPath = "InputData/Phantom/Al";     // ������ʱ�ɲ����룬���������ܶ�·��
    mFilePath.gridPath = "InputData/Grid/Pb_[1_80]keV_1keV.txt";
    mFilePath.scintillatorPath = "InputData/Scintillator/CsI_[1_80]keV_1keV.txt";
    mFilePath.convKernelGridPath = "InputData/ConvKernel/ConvKernelGrid_601.raw";
    
}

void ForwardProjection::inCTScanInfoGrid()
{
    // ���·��
    mFilePath.IEPath = "OutputResult/ForwardProjection/Grid";

    // դ��Ϣ
    mGridInfo.FD                        = 1400;         // mm
    mGridInfo.h                         = 2;            // mm
    mGridInfo.leadStripsWidth           = 2780;          // um
    mGridInfo.leadStripsDistance        = 2780;           // um
    mGridInfo.materialGridInterspacer   = "Air";
    mGridInfo.materialGridStrip         = "Pb";
    mGridInfo.rhoGS                     = 1.135E+01;    // g/cm^3
    mGridInfo.uGridStrip                = 0;            // �����и��������Զ���ֵ
    mGridInfo.uInterspacer              = 0;            //
}

void ForwardProjection::inCTScanInfoGridNoResponse()
{
    mFilePath.IEPath = "OutputResult/ForwardProjection/Grid/NoResponse";
}

void ForwardProjection::inCTScanInfoSinMat()
{
    mFilePath.phantomMAttenPath = "InputData/Phantom/Al_[1_80]keV_1keV.dat";
    mFilePath.phantomPath = "InputData/Phantom/Al_Density_256x256x256_float.raw";
}

void ForwardProjection::inCTScanInfoNoGrid()
{
    //mFilePath.convKernelGridPath    = "InputData/ConvKernel/ConvKernelGrid_601.raw";
    mFilePath.IEPath                = "OutputResult/ForwardProjection/NoGrid";
}

void ForwardProjection::inCTScanInfoNoGridNoResponse()
{
    mFilePath.IEPath = "OutputResult/ForwardProjection/NoGrid/NoResponse";
}
