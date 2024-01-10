#include "ForwardProjection.h"



//#define PHOTONS 60764.1
//#define PHOTONS 8*1e4
#define PHOTONS 6.9*1e4
#define PNUM 512
#define DNUM 600
#define PHANTOMLEN 102.4//51.2   // mm



ForwardProjection::~ForwardProjection()
{
    DELETEARR(_mPhantomMaterial.density);
    DELETEARR(_mPhantomMaterial.massAttenPath);
}

void ForwardProjection::run()
{
    _getCTScanConfig();
    _getOutputFolderConfig();

    if (!_isTestI0)
    {
        if (_isDetResponse)
        {
            if (_isGrid)
            {
                if (_isSingleMaterial)
                {
                    forwardSinMatPolyProjGrid();
                }
                else
                {
                    forwardPolyProjGrid();
                }
            }
            else
            {
                if (_isSingleMaterial)
                {
                    forwardSinMatPolyProjNoGrid();
                }
                else
                {
                    forwardPolyProjNoGrid();
                }
            }

        }
        else
        {
            if (_isGrid)
            {
                if (_isSingleMaterial)
                {
                    forwardSinMatPolyProjGridNoResponse();
                }
                else
                {  //����Ӧ����դ�������
                
                }
            }
            else
            {
                if (_isSingleMaterial)
                {
                    forwardSinMatPolyProjNoGridNoResponse();
                }
                else
                {  //����Ӧ����դ�������
 
                }
            }
        }
        
    }
    else
    {
        forwarProjNoGridTestI0();
    }
}

void ForwardProjection::forwardPolyProjGrid()
{
    //inCTScanInfo();
    //inCTScanInfoGrid();

    _getCTScanGridConfig();
    _getCTScanDetResponseConfig();

    CBCTPolyForwardProjGrid mCBCTPolyForwardProjGrid = CBCTPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath, _mPhantomMaterial);
    mCBCTPolyForwardProjGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjGrid()
{
    /*inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();*/
    
    _getCTScanGridConfig();
    _getCTScanDetResponseConfig();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjGridFoSp()
{
    /*inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();*/

    _getCTScanGridConfig();
    _getCTScanDetResponseConfig();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProjFoSp();
}

void ForwardProjection::forwardSinMatPolyProjGridNoResponse()
{
    /*inCTScanInfo();
    inCTScanInfoGrid();
    inCTScanInfoSinMat();
    inCTScanInfoGridNoResponse();*/

    _getCTScanGridConfig();

    CBCTSinMatPolyForwardProjGrid mCBCTSinMatPolyForwardProjGrid = CBCTSinMatPolyForwardProjGrid(mCTScanParas, mGridInfo, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjGrid.computePolyForwProjNoResponse();
}

void ForwardProjection::forwardPolyProjNoGrid()
{
    /*inCTScanInfo();
    inCTScanInfoNoGrid();*/

    _getCTScanDetResponseConfig();

    CBCTPolyForwardProjNoGrid mCBCTPolyForwardProjNoGrid = CBCTPolyForwardProjNoGrid(mCTScanParas, mFilePath, _mPhantomMaterial);
    mCBCTPolyForwardProjNoGrid.computePolyForwProj();
}

void ForwardProjection::forwardSinMatPolyProjNoGrid()
{
    /*inCTScanInfo();
    inCTScanInfoNoGrid();
    inCTScanInfoSinMat();*/

    _getCTScanDetResponseConfig();

    CBCTSinMatPolyForwardProjNoGrid mCBCTSinMatPolyForwardProjNoGrid = CBCTSinMatPolyForwardProjNoGrid(mCTScanParas, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProj();

}

void ForwardProjection::forwardSinMatPolyProjNoGridFoSp()
{
    /*inCTScanInfo();
    inCTScanInfoNoGrid();
    inCTScanInfoSinMat();*/

    _getCTScanDetResponseConfig();

    CBCTSinMatPolyForwardProjNoGrid mCBCTSinMatPolyForwardProjNoGrid = CBCTSinMatPolyForwardProjNoGrid(mCTScanParas, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProjFoSp();
}

void ForwardProjection::forwardSinMatPolyProjNoGridNoResponse()
{
    /*inCTScanInfo();
    inCTScanInfoNoGrid();
    inCTScanInfoSinMat();
    inCTScanInfoNoGridNoResponse();*/

    CBCTSinMatPolyForwardProjNoGrid mCBCTSinMatPolyForwardProjNoGrid = CBCTSinMatPolyForwardProjNoGrid(mCTScanParas, mFilePath, _mPhantomMaterial);
    mCBCTSinMatPolyForwardProjNoGrid.computePolyForwProjNoResponse();

}

void ForwardProjection::forwarProjNoGridTestI0()
{

    _getCTScanDetResponseConfig();

    /*inCTScanInfo();
    inCTScanInfoNoGrid();*/

    CBCTPolyForwardProjNoGrid mCBCTPolyForwardProjNoGrid = CBCTPolyForwardProjNoGrid(mCTScanParas, mFilePath, _mPhantomMaterial);
    mCBCTPolyForwardProjNoGrid.testI0();
}

void ForwardProjection::_getCTScanConfig()
{
    

    string fileName = "Config/ScanSetting.ini";
    mCInitFile.GetFileName(fileName.c_str());

    char section[32];

    /* ******* ɨ��������� ******* */
    sprintf_s(section, sizeof(section),"%s", "ScanPara");

    mCInitFile.GetEntryValue(section, "IsTestI0", -100, _isTestI0);
    //mCInitFile.GetEntryValue(section, "IsPolyenergetic", -100, _isPolyenergetic);
    mCInitFile.GetEntryValue(section, "IsDetResponse", -100, _isDetResponse);
    mCInitFile.GetEntryValue(section, "IsGrid", -100, _isGrid);
    mCInitFile.GetEntryValue(section, "IsSingleMaterial", -100, _isSingleMaterial);
    mCInitFile.GetEntryValue(section, "SDD", -100, mCTScanParas.sdd);
    mCInitFile.GetEntryValue(section, "SOD", -100, mCTScanParas.sod);
    mCInitFile.GetEntryValue(section, "ScanAngleNum", -100, mCTScanParas.projNum);
    mCInitFile.GetEntryValue(section, "RotatedDirection", -100, mCTScanParas.rotatedDirection);

    if (_isTestI0 == -100 || /*_isPolyenergetic == -100 ||*/ _isDetResponse == -100 || _isGrid == -100 || _isSingleMaterial == -100 || mCTScanParas.sdd < 0 || mCTScanParas.sod < 0 || mCTScanParas.projNum < 0 || mCTScanParas.rotatedDirection == -100)
    {
        cout << "ɨ�������ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }

    /* ******* ����Դ���� ******* */ 
    sprintf_s(section, sizeof(section), "%s", "RaySource");

    mCInitFile.GetEntryValue(section, "Photons", -100, mCTScanParas.I0Val);
    mCInitFile.GetEntryValue(section, "SpecEnergyNum", -100, mCTScanParas.specEnergyNum);
    mCInitFile.GetEntryValue(section, "SpectrumPath", "null", mFilePath.spectrumPath);


    if (mCTScanParas.I0Val == -100 || mCTScanParas.specEnergyNum == -100 || mFilePath.spectrumPath == "null")
    {
        cout << "����Դ�����ò�����ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }


    /* ******* ģ���С�����أ����� ******* */
    sprintf_s(section, sizeof(section), "%s", "PhantomVox");

    mCInitFile.GetEntryValue(section, "PNumX", -100, mCTScanParas.pNumX);
    mCInitFile.GetEntryValue(section, "PNumY", -100, mCTScanParas.pNumY);
    mCInitFile.GetEntryValue(section, "PNumZ", -100, mCTScanParas.pNumZ);

    if (mCTScanParas.pNumX < 0 || mCTScanParas.pNumY < 0 || mCTScanParas.pNumZ < 0)
    {
        cout << "ģ�����ش�С��ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }

    mCInitFile.GetEntryValue(section, "PLengthX", -100, mCTScanParas.pLengthX);
    mCInitFile.GetEntryValue(section, "PLengthY", -100, mCTScanParas.pLengthY);
    mCInitFile.GetEntryValue(section, "PLengthZ", -100, mCTScanParas.pLengthZ);

    if (mCTScanParas.pLengthX < 0 || mCTScanParas.pLengthY < 0 || mCTScanParas.pLengthZ < 0)
    {
        cout << "ͼ��ߴ��С��ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }


    /* ******* ̽�������� ******* */
    sprintf_s(section, sizeof(section), "%s", "Detector");
    mCInitFile.GetEntryValue(section, "DNumU", -100, mCTScanParas.dNumU);
    mCInitFile.GetEntryValue(section, "DNumV", -100, mCTScanParas.dNumV);
    mCInitFile.GetEntryValue(section, "DSize", -100, mCTScanParas.dSize);


    if (mCTScanParas.dNumU == -100 || mCTScanParas.dNumV == -100 || mCTScanParas.dSize < 0)
    {
        cout << "̽���������ò�����ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }

    /* ******* ģ��������� ******* */
    sprintf_s(section, sizeof(section), "%s", "PhantomMaterial");

    mCInitFile.GetEntryValue(section, "MaterialNum", -100, _mPhantomMaterial.materialNum);
    mCInitFile.GetEntryValue(section, "StructureMarkPath", "null", _mPhantomMaterial.structureMarkPath);
    
    if (_mPhantomMaterial.materialNum < 0 || _mPhantomMaterial.structureMarkPath == "null")
    {
        cout << "ģ����ʸ�����ṹ����ļ�·�������ò�����ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }

    _mPhantomMaterial.density = new float[_mPhantomMaterial.materialNum];
    _mPhantomMaterial.massAttenPath = new string[_mPhantomMaterial.materialNum];

    char massAttenKey[64], densityKey[32];

    for (int i = 0; i < _mPhantomMaterial.materialNum; i++)
    {
        sprintf_s(densityKey, sizeof(densityKey), "Density%d", i + 1);
        sprintf_s(massAttenKey, sizeof(massAttenKey), "MassAttenPath%d", i + 1);
        
        mCInitFile.GetEntryValue(section, densityKey, -100, _mPhantomMaterial.density[i]);
        mCInitFile.GetEntryValue(section, massAttenKey, "null", _mPhantomMaterial.massAttenPath[i]);
        
        if (_mPhantomMaterial.density[i] < 0 || _mPhantomMaterial.massAttenPath[i] == "null")
        {
            cout << "��" << i+1 << "�����ʵ��ܶȻ�����˥��ϵ���ļ�·�������ò�����ȡ����" << endl;
            system("pause");
            exit(0);
            return;
        }
    }

}

/* ******* ��������դ���� ******* */
void ForwardProjection::_getCTScanGridConfig()
{
    char section[32];
    sprintf_s(section, sizeof(section), "%s", "PostGrid");

    mCInitFile.GetEntryValue(section, "FD", -100, mGridInfo.FD);
    mCInitFile.GetEntryValue(section, "Height", -100, mGridInfo.h);
    mCInitFile.GetEntryValue(section, "GridStripsWidth", -100, mGridInfo.gridStripsWidth);
    mCInitFile.GetEntryValue(section, "InterspcaceWidth", -100, mGridInfo.interspcaceWidth);
    mCInitFile.GetEntryValue(section, "GridStripMaterial", "null", mGridInfo.gridStripMaterial);
    mCInitFile.GetEntryValue(section, "InterspaceMaterial", "null", mGridInfo.interspaceMaterial);
    mCInitFile.GetEntryValue(section, "GridStripDensity", -100, mGridInfo.gridStripDensity);
    mCInitFile.GetEntryValue(section, "InterspaceDensity", -100, mGridInfo.interspaceDensity);
    mCInitFile.GetEntryValue(section, "GridStripMAPath", "null", mGridInfo.gridStripMAPath);
    mCInitFile.GetEntryValue(section, "InterspaceMAPath", "null", mGridInfo.interspaceMAPath);

    if (mGridInfo.FD < 0 || mGridInfo.h < 0 || mGridInfo.gridStripsWidth < 0 || mGridInfo.interspcaceWidth < 0 || mGridInfo.interspaceMaterial == "null" || mGridInfo.gridStripMaterial == "null" || mGridInfo.gridStripDensity < 0 || mGridInfo.interspaceDensity < 0 || mGridInfo.gridStripMAPath == "null" || mGridInfo.interspaceMAPath == "null")
    {
        cout << "��������դ�����ò�����ȡ����" << endl;
        system("pause");
        exit(0);
    }
}


/* ******* ��˸������ ******* */
void ForwardProjection::_getCTScanDetResponseConfig()
{
    char section[32];
    sprintf_s(section, sizeof(section), "%s", "Scintillator");
    mCInitFile.GetEntryValue(section, "Density", -100, mCTScanParas.mScintilltorInfo.scintillatorDensity);
    mCInitFile.GetEntryValue(section, "Thickness", -100, mCTScanParas.mScintilltorInfo.scintillatorThickness);
    mCInitFile.GetEntryValue(section, "ThicknessErr", -100, mCTScanParas.mScintilltorInfo.scintillatorThicknessErr);
    mCInitFile.GetEntryValue(section, "PhotoelectConverFactor", -100, mCTScanParas.mScintilltorInfo.detResponseFactor);
    mCInitFile.GetEntryValue(section, "ScintillatorMAPath", "null", mFilePath.scintillatorPath);

    if (mCTScanParas.mScintilltorInfo.scintillatorDensity < 0 || mCTScanParas.mScintilltorInfo.scintillatorThickness < 0 || mCTScanParas.mScintilltorInfo.scintillatorThicknessErr < 0 || mCTScanParas.mScintilltorInfo.detResponseFactor < 0 || mCTScanParas.mScintilltorInfo.detResponseFactor >1 || mFilePath.scintillatorPath == "null")
    {
        cout << "��˸������ò�����ȡ����" << endl;
        system("pause");
        exit(0);
        return;
    }
}


void ForwardProjection::_getOutputFolderConfig()
{
    CInitFile tempCInitFile;

    string fileName = "Config/OutputFolder.ini";
    tempCInitFile.GetFileName(fileName.c_str());

    char tempSection[32];

    if (!_isTestI0)
    {
        if (_isGrid)
        {
            if (_isDetResponse)
            {
                sprintf_s(tempSection, sizeof(tempSection), "%s", "GridResponse");
            }
            else
            {
                sprintf_s(tempSection, sizeof(tempSection), "%s", "GridNoResponse");
            }
        }
        else
        {
            if (_isDetResponse)
            {
                sprintf_s(tempSection, sizeof(tempSection), "%s", "NoGridResponse");
            }
            else
            {
                sprintf_s(tempSection, sizeof(tempSection), "%s", "NoGridNoResponse");
            }
        }
    }
    else
    {
        sprintf_s(tempSection, sizeof(tempSection), "%s", "NoGridResponseTestI0");
    }

    tempCInitFile.GetEntryValue(tempSection, "OutFolder", "null", mFilePath.outputFolder);

    if (mFilePath.outputFolder == "null")
    {
        cout << "������·�������ò�����ȡ����" << endl;
        system("pause");
        exit(0);
    }
}

void ForwardProjection::inCTScanInfo()
{
    // ɨ��ϵͳ����
    // ģ����Ϣ
    mCTScanParas.pNumX = PNUM;
    mCTScanParas.pNumY = PNUM;
    mCTScanParas.pNumZ = PNUM;
	mCTScanParas.pLengthX = PHANTOMLEN;
	mCTScanParas.pLengthY = PHANTOMLEN;
	mCTScanParas.pLengthZ = PHANTOMLEN;

    // ̽�������
    mCTScanParas.dNumU = DNUM;
    mCTScanParas.dNumV = DNUM; 
    mCTScanParas.dSize = 0.278;  //mm

    // ɨ�����
    mCTScanParas.projNum = 360;
    mCTScanParas.rotatedDirection = 1;      // -1 ���� ˳ʱ��, 1 ���� ��ʱ��
    mCTScanParas.sdd = 1050;
    mCTScanParas.sod = 630;

    // ���߽��㣬Ŀǰ��դ�Ľ���ģ������
    mCTScanParas.focalSpotSize = 0.4;   // mm

    mCTScanParas.I0Val = PHOTONS;

    // ������ɢ
    mCTScanParas.specEnergyNum = 80;         // ������ɢ������������� 1 ʱ������ֵ�ϵ��ڹܵ�ѹ
    mCTScanParas.spectrumStep = 1;

    // ��˸����Ϣ
    mCTScanParas.mScintilltorInfo.scintillatorDensity = 4.510;      // �ܶ�   g/cm^3
    mCTScanParas.mScintilltorInfo.scintillatorThickness = 0.4;      // ���   mm
    mCTScanParas.mScintilltorInfo.scintillatorThicknessErr = 0.02;   // ������� mm
    mCTScanParas.mScintilltorInfo.detResponseFactor = 1.0f;         // ��Ӧ����

    // ·��
    mFilePath.spectrumPath = "InputData/Spectrum/Spectrum_80kVp_1mAs_1mmAl_GEMaxiray_1keV_Normal.dat";
    mFilePath.phantomPath = "InputData/Phantom/Al";     // ������ʱ�ɲ����룬���������ܶ�·��
    mFilePath.gridPath = "InputData/Grid/Pb_[1_80]keV_1keV.txt";
    mFilePath.scintillatorPath = "InputData/Scintillator/CsI_[1_80]keV_1keV.txt";
    mFilePath.convKernelGridPath = "InputData/ConvKernel/ConvKernelGrid_601.raw";
    
}

void ForwardProjection::inCTScanInfoGrid()
{
    // ���·��
    mFilePath.outputFolder = "OutputResult/ForwardProjection/Grid";

    // դ��Ϣ
    mGridInfo.FD                        = 1400;         // mm
    mGridInfo.h                         = 1.0f;          // mm
    mGridInfo.gridStripsWidth           = 42;          // դ����� um
    mGridInfo.interspcaceWidth        = 208;          // դ������ um
    mGridInfo.interspaceMaterial    = "Air";
    mGridInfo.gridStripMaterial         = "Pb";
    mGridInfo.gridStripDensity                     = 1.135E+01;    // g/cm^3
    mGridInfo.gridStripLinearAtten                = 0;            // �����и��������Զ���ֵ
    mGridInfo.interspaceLinearAtten              = 0;            //
}

void ForwardProjection::inCTScanInfoGridNoResponse()
{
    mFilePath.outputFolder = "OutputResult/ForwardProjection/Grid/NoResponse";
}

void ForwardProjection::inCTScanInfoSinMat()
{
    mFilePath.phantomMAttenPath = "InputData/Phantom/Water_[1_150]keV_1keV.dat";
    mFilePath.phantomPath = "InputData/Phantom/WaterSphere_Density_512x512x512_float.raw";
}

void ForwardProjection::inCTScanInfoNoGrid()
{
    //mFilePath.convKernelGridPath    = "InputData/ConvKernel/ConvKernelGrid_601.raw";
    mFilePath.outputFolder                = "OutputResult/ForwardProjection/NoGrid";
}

void ForwardProjection::inCTScanInfoNoGridNoResponse()
{
    mFilePath.outputFolder = "OutputResult/ForwardProjection/NoGrid/NoResponse";
}

void ForwardProjection::inCTScanInfoNoGridTestI0()
{
    mFilePath.outputFolder = "OutputResult/ForwardProjection/NoGrid";
}
