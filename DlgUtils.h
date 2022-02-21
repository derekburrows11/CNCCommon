//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DLGUTILS_H__EC661548_4908_4F18_BC45_AC239034653A__INCLUDED_)
#define AFX_DLGUTILS_H__EC661548_4908_4F18_BC45_AC239034653A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000



// Definitions in "vc\mfc\src\dlgfloat.cpp"

//static BOOL AFXAPI AfxSimpleFloatParse(LPCTSTR lpszText, double& d);

void AFXAPI AfxTextFloatFormat(CDataExchange* pDX, int nIDC,
	void* pData, double value, int nSizeGcvt);

// Definitions in "DlgUtils.cpp"
void AFXAPI AfxTextFloatFormat(CDataExchange* pDX, int nIDC,
	void* pData, double value, int nSizeGcvt);


/*
void AFXAPI DDX_Text(CDataExchange* pDX, int nIDC, double& value, int nSizeGcvt)
{
	AfxTextFloatFormat(pDX, nIDC, &value, value, nSizeGcvt);
}
*/

void AFXAPI DDX_TextFix(CDataExchange* pDX, int nIDC, double& value, int nDigits);
void AFXAPI DDX_TextFix(CDataExchange* pDX, int nIDC, float& value, int nDigits);


#endif // !defined(AFX_DLGUTILS_H__EC661548_4908_4F18_BC45_AC239034653A__INCLUDED_)
