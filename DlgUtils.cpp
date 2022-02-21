//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <float.h>              // floating point precision

#include "DlgUtils.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif





static BOOL AFXAPI AfxSimpleFloatParse(LPCTSTR lpszText, double& d)
{
	ASSERT(lpszText != NULL);
	while (*lpszText == ' ' || *lpszText == '\t')
		lpszText++;

	TCHAR chFirst = lpszText[0];
	d = _tcstod(lpszText, (LPTSTR*)&lpszText);
	if (d == 0.0 && chFirst != '0')
		return FALSE;   // could not convert
	while (*lpszText == ' ' || *lpszText == '\t')
		lpszText++;

	if (*lpszText != '\0')
		return FALSE;   // not terminated properly

	return TRUE;
}


void AFXAPI DDX_TextFix(CDataExchange* pDX, int nIDC,
	double& value, int nDigits)
{
	HWND hWndCtrl = pDX->PrepareEditCtrl(nIDC);
	TCHAR szBuffer[32];
	if (pDX->m_bSaveAndValidate)
	{
		::GetWindowText(hWndCtrl, szBuffer, sizeof(szBuffer));
		double d;
		if (!AfxSimpleFloatParse(szBuffer, d))
		{
			AfxMessageBox(AFX_IDP_PARSE_REAL);
			pDX->Fail();            // throws exception
		}
		value = d;
	}
	else
	{
		_stprintf(szBuffer, _T("%.*f"), nDigits, value);
		::SetWindowText(hWndCtrl, szBuffer);
	}
}

void AFXAPI DDX_TextFix(CDataExchange* pDX, int nIDC,
	float& value, int nDigits)
{
	HWND hWndCtrl = pDX->PrepareEditCtrl(nIDC);
	TCHAR szBuffer[32];
	if (pDX->m_bSaveAndValidate)
	{
		::GetWindowText(hWndCtrl, szBuffer, sizeof(szBuffer));
		double d;
		if (!AfxSimpleFloatParse(szBuffer, d))
		{
			AfxMessageBox(AFX_IDP_PARSE_REAL);
			pDX->Fail();            // throws exception
		}
		value = (float)d;
	}
	else
	{
		_stprintf(szBuffer, _T("%.*f"), nDigits, value);
		::SetWindowText(hWndCtrl, szBuffer);
	}
}

