#include "ForwardProjection.h"


//#define PHOTONS 60764.1
#define PHOTONS 8*1e4
//#define PHOTONS 6.9*1e4
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

void ForwardProjection::forwardSinMatPolyProjGridFoSp()
{
    inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProjFoSp();
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
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProjNoResponse();

}

void ForwardProjection::forwarProjNoGridTestI0()
{
    inCTScanInfo();
    inCTScanInfoNoGrid();

    CBCTPolyForwardProjNoGrid mCBCTPolyForwardProjNoGrid = CBCTPolyForwardProjNoGrid(mCTScanParas, mFilePath);
    mCBCTPolyForwardProjNoGrid.testI0();
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

    mCTScanParas.focalSpotSize = 0.6;   // mm

    mCTScanParas.I0Val = PHOTONS;

    mCTScanParas.specEnergyNum = 1;         // ������ɢ����
    mCTScanParas.spectrumStep = 73;

    // ��˸����Ϣ
    mCTScanParas.mScintilltorInfo.scintillatorDensity = 4.510;      // �ܶ�   g/cm^3
    mCTScanParas.mScintilltorInfo.scintillatorThickness = 0.4;      // ���   mm
    mCTScanParas.mScintilltorInfo.scintillatorThicknessErr = 0.05;   // ������� mm
    mCTScanParas.mScintilltorInfo.detResponseFactor = 1.0f;         // ��Ӧ����

    // ·��
    mFilePath.spectrumPath = "InputData/Spectrum/Spectrum_60keV_1mmAl.txt";
    mFilePath.phantomPath = "InputData/Phantom/Al";     // ������ʱ�ɲ����룬���������ܶ�·��
    mFilePath.gridPath = "InputData/Grid/Pb_60keV.txt";
    mFilePath.scintillatorPath = "InputData/Scintillator/CsI_60keV.txt";
    mFilePath.convKernelGridPath = "InputData/ConvKernel/ConvKernelGrid_601.raw";
    
}

void ForwardProjection::inCTScanInfoGrid()
{
    // ���·��
    mFilePath.IEPath = "OutputResult/ForwardProjection/Grid";

    // դ��Ϣ
    mGridInfo.FD                        = 1400;         // mm
    mGridInfo.h                         = 1.5f;          // mm
    mGridInfo.leadStripsWidth           = 139;          // դ����� um
    mGridInfo.leadStripsDistance        = 139;           // դ������ um
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
    mFilePath.phantomMAttenPath = "InputData/Phantom/Al_60keV.dat";
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

void ForwardProjection::inCTScanInfoNoGridTestI0()
{
    mFilePath.IEPath = "OutputResult/ForwardProjection/NoGrid";
}
