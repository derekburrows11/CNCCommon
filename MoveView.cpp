// MoveView.cpp : implementation file
//

#include "stdafx.h"
#include <afxtempl.h>		// just for Z-order CList

#include "CommonResource.h"
//#include <math.h>

//#include "MainFrm.h"		// just to use CMainFrame::ShowScale() - sends message instead!
#include "Colors.h"

#include "MoveView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif






/////////////////////////////////////////////////////////////////////////////
// Utility functions


template<class T>
int ToCPoints(const TMatrix<T>& mx, CPoint* arPt, int iMax);

template<class T>
int ToCPoints(const TMatrix<T>& mx, CPoint* arPt, int iMax)
{
	ASSERT(mx.w>=1 && mx.w<=3);
	if (iMax == -1 || iMax > mx.h)		// if iMax == -1 (default) does all rows
		iMax = mx.h;
	int mxIdx = 0;
	for (int i = 0; i < iMax; i++)
	{
		arPt[i].x = NEARINT(mx[mxIdx + 0]);
		arPt[i].y = NEARINT(mx[mxIdx + 1]);
		mxIdx += mx.w;
	}
	return iMax;
}


/////////////////////////////////////////////////////////////////////////////
// DrawNode Utility function

void DrawNodes(CDC* pDC, const SDrawNodes& dn)
{
	int r = dn.iRadius;		// radius of nodes
	int rp1 = r + 1;
	int iNodeLoc = dn.iStart % dn.iPeriod;		// a 0 is end node, non zero are control
	int nNodeType;
	CPoint arPt[2];
	for (int i = 0; i < dn.numPts; i++)
	{
		if (iNodeLoc == 0)		// end node
			nNodeType = dn.nEndType;
		else
			nNodeType = dn.nControlType;
		if (++iNodeLoc == dn.iPeriod)
			iNodeLoc = 0;
		if (nNodeType == 0)
			continue;
//		const CPoint& pt = dn.arPts[i];
		int x = dn.arPts[i].x;
		int y = dn.arPts[i].y;
		switch (nNodeType)
		{
		case 0: break;
		case 1:
			arPt[0].x = arPt[1].x = x;		// draws a point
			arPt[0].y = y;
			arPt[1].y = y+1;
			pDC->Polyline(arPt, 2);
			break;
		case 2:
			pDC->Ellipse(x-r, y-r, x+rp1, y+rp1);
			break;
		case 6:
			arPt[0].x = x-r;		// draws a x
			arPt[0].y = y-r;
			arPt[1].x = x+rp1;
			arPt[1].y = y+rp1;
			pDC->Polyline(arPt, 2);
			arPt[0].y = y+r;
			arPt[1].y = y-rp1;
			pDC->Polyline(arPt, 2);
			break;
		case 7:
			arPt[0].x = x-r;		// draws a +
			arPt[0].y = y;
			arPt[1].x = x+rp1;
			arPt[1].y = y;
			pDC->Polyline(arPt, 2);
			arPt[0].x = x;
			arPt[0].y = y-r;
			arPt[1].x = x;
			arPt[1].y = y+rp1;
			pDC->Polyline(arPt, 2);
			break;
		}
	}
}





/////////////////////////////////////////////////////////////////////////////
// CMoveView

IMPLEMENT_DYNCREATE(CMoveView, CScrollView)

BEGIN_MESSAGE_MAP(CMoveView, CScrollView)
	//{{AFX_MSG_MAP(CMoveView)
	ON_WM_ERASEBKGND()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_MBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONUP()
	ON_WM_MBUTTONUP()
	ON_WM_MOUSEWHEEL()
	ON_WM_CREATE()
	ON_WM_HSCROLL()
	ON_WM_VSCROLL()
	ON_WM_SETCURSOR()
	ON_WM_SETFOCUS()
	ON_WM_KILLFOCUS()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


CVector CMoveView::ms_ptFocusRealStored = 0;
CVect2 CMoveView::ms_ScaleStored = 0;

/////////////////////////////////////////////////////////////////////////////
// CMoveView construction/destruction

CMoveView::CMoveView()
{
	m_hCursor = NULL;
	m_hCursorMove = NULL;
	m_hCursorDir = NULL;

	// constant scaling factors
	m_mm2dev = 112.0 / 25.4;		// 112dpi / 25.4mmpi			19" screen @1600x1200 is 112dpi
	m_mm2dev.y = -m_mm2dev.y;		// device y is pos down
	m_dev2mm = 25.4 / 112.0;		// 25.4mmpi / 112dpi
	m_dev2mm.y = -m_dev2mm.y;		// device y is pos down

	// if using logical units of 0.1mm		MM_LOMETRIC
/*
	m_log2dev = 96.0 / (25.4 * 10);	// 96dpi / (10*25.4mmpi) pixel/log		windows uses 96dpi for pixel size!
	m_log2dev.y = -m_log2dev.y;		// device y is pos down
	m_dev2log = (25.4 * 10) / 96.0;	// (10*25.4mmpi) / 96dpi log/pixel
	m_dev2log.y = -m_dev2log.y;		// device y is pos down
*/
	// if using logical units of 1 pixel	MM_TEXT
	m_log2dev = 1;			// for logical units of MM_TEXT
	m_dev2log = 1;

	m_mm2log = m_mm2dev * m_dev2log;
	m_log2mm = m_log2dev * m_dev2mm;

	m_maxLogicalCoord = 0x7000;					// max logical coord with 16 bit values for win98 line functions

	m_bIsoScaling = false;
	m_b3DView = false;
	m_bShowRegionBox = true;

	// dynamic movement variables
	m_nDynMoveType = DYNA_NONE;
	m_DynaTranslateAmp = 1.0;		// for amplified translations (default = 1.0)

	m_bDevSizeOverflow = false;

	m_nShowAxes = 0;

}

CMoveView::~CMoveView()
{
}

/////////////////////////////////////////////////////////////////////////////
// CMoveView diagnostics

#ifdef _DEBUG
void CMoveView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CMoveView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMoveView drawing

void CMoveView::OnInitialUpdate()
{
/*
	DEVMODE dm;
	DWORD flags;

	dm.dmSize = sizeof(dm);
	dm.dmLogPixels = 112;
	dm.dmFields = DM_LOGPIXELS;
	flags = CDS_UPDATEREGISTRY;
//	int retVal = ChangeDisplaySettings(&dm, flags);
*/	


	
	// make sure view client size reports correctly
	GetParentFrame()->RecalcLayout();	// recalc after any new or changed toolbars etc.
	CRect rectWnd;
	GetClientRect(rectWnd);		// without scroll bars, in device coords (pixels)
	CScrollView::SetScrollSizes(MM_TEXT, rectWnd.Size() + CSize(1,1));		// make sure scroll bars are on for future client size reporting
//	m_nMapMode = MM_LOMETRIC;	// using text mode as logical coords


///////////
	m_ptFocusReal = 0;
	m_ptMinRegion = 0;
	m_ptMaxRegion = 0;

	m_ViewAzimuth = 0;	// gives default x->x, y->y, z->z (above screen)
	m_ViewElevation = 90;
	m_ViewOrientation = 0;
	SetTransform();
	SetScale(1);			// default scale

	m_ptViewBorder.x = 10;	// mm
	m_ptViewBorder.y = 10;	// mm

	UpdateTotalViewSize();

/////////
	m_hCursorMove = AfxGetApp()->LoadStandardCursor(IDC_SIZEALL);
	m_hCursorDir = AfxGetApp()->LoadStandardCursor(IDC_CROSS);


	CScrollView::OnInitialUpdate();
}

void CMoveView::OnPrepareDC(CDC* pDC, CPrintInfo* pInfo) 
{
	// TODO: Add your specialized code here and/or call the base class
	CScrollView::OnPrepareDC(pDC, pInfo);		// sets SetViewportOrg(...) from scroll position
//	pDC->SetViewportOrg(CPoint(0, 0));	// reset device coords of Viewport which are used for scrolling, m_ptWindowULReal used instead for window reference
	ASSERT(pDC->GetWindowOrg() == CPoint(0,0));		// logical coords
	SetViewportOrgToFocus(pDC);
}

/*
struct SPointArrayInfo
{
	CPoint* arPts;
	double* arElems;
	int sizeElem;
	int numPts;

	int idxPtsInit;
	int numPtsDraw;
};
*/

//void CMoveView::SetLogicalXValues(SPointArrayInfo& ai)
void CMoveView::SetLogicalXValues(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 int& idxPtsInit, int& numPtsToDraw)
{
	const double xReal2Log = m_real2log.x;		// a local copy - maybe faster?
	const double xFocusReal = m_ptFocusReal.x;
	const double maxRealCoordX = m_log2real.x * m_maxLogicalCoord;		// 8192 pixels (~5 display widths) either side
	const char* pElem = (const char*)arElems;			// char* for 1 byte size! don't want 8 byte double

	int nSect = -1;				// -1 indicates before range, 0 in range, 1 after range
	idxPtsInit = 0;
	numPtsToDraw = 0;
	for (int i = 0; i < numPts; i++, pElem += sizeElem)
	{
		double xReal = *(double*)pElem - xFocusReal;
		if (fabs(xReal) <= maxRealCoordX)				// point in range	
		{
			if (nSect == -1)		// first point in range
			{
				nSect = 0;
				idxPtsInit = i;
			}
			arPts[i].x = NEARINT(xReal2Log * xReal);
		}
		else					// don't use point - out of range!
			if (nSect == 0)		// section was in range, now out of range
				break;
	}
	if (nSect == 0)
		numPtsToDraw = i - idxPtsInit;
}

void CMoveView::SetLogicalXValuesInRange(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 double minX, double maxX, int& idxPtsInit, int& numPtsToDraw)
{
	const double xReal2Log = m_real2log.x;		// a local copy - maybe faster?
	const double xFocusReal = m_ptFocusReal.x;
	const double maxRealCoordX = m_log2real.x * m_maxLogicalCoord;		// 8192 pixels (~5 display widths) either side
	const char* pElem = (const char*)arElems;			// char* for 1 byte size! don't want 8 byte double

	const double minXReal = max(minX, xFocusReal - maxRealCoordX);
	const double maxXReal = min(maxX, xFocusReal + maxRealCoordX);

	int nSect = -1;				// -1 indicates before range, 0 in range, 1 after range
	idxPtsInit = 0;
	numPtsToDraw = 0;
	for (int i = 0; i < numPts; i++, pElem += sizeElem)
	{
		double xElem = *(double*)pElem;
		if (xElem >= minXReal && xElem <= maxXReal)				// point in range	
		{
			if (nSect == -1)		// first point in range
			{
				nSect = 0;
				idxPtsInit = i;
			}
			arPts[i].x = NEARINT(xReal2Log * (xElem - xFocusReal));
		}
		else					// don't use point - out of range!
			if (nSect == 0)		// section was in range, now out of range
				break;
	}
	if (nSect == 0)
		numPtsToDraw = i - idxPtsInit;
}

void CMoveView::SetLogicalXBezValues(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 int& idxPtsInit, int& numPtsToDraw)
{
	const double xReal2Log = m_real2log.x;		// a local copy - maybe faster?
	const double xFocusReal = m_ptFocusReal.x;
	const double maxRealCoordX = m_log2real.x * m_maxLogicalCoord;		// 8192 pixels (~5 display widths) either side
	const char* pElem = (const char*)arElems;			// char* for 1 byte size! don't want 8 byte double

	int nSect = -1;				// -1 indicates before range, 0 in range, 1 after range
	idxPtsInit = 0;
	numPtsToDraw = 0;
	for (int i = 0; i < numPts; i++, pElem += sizeElem)
	{
		double xReal = *(double*)pElem - xFocusReal;
		if (fabs(xReal) <= maxRealCoordX)				// point in range	
		{
			if (nSect == -1)		// first point in range
			{
				nSect = 0;
				idxPtsInit = i;
			}
			arPts[i].x = NEARINT(xReal2Log * xReal);
		}
		else					// don't use point - out of range!
		{
			if (nSect == -1)
			{
				i += 2;			// skip to next end node of bezier for test
				pElem += 2 * sizeElem;
			}
			else if (nSect == 0)		// section was in range, now out of range
				break;
		}
	}
	if (nSect == 0)
		numPtsToDraw = ((i - idxPtsInit - 1) / 3) * 3 + 1;		// move end point back to nearest bezier end node
}

void CMoveView::SetLogicalXBezValuesInRange(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 double minX, double maxX, int& idxPtsInit, int& numPtsToDraw)
{
	const double xReal2Log = m_real2log.x;		// a local copy - maybe faster?
	const double xFocusReal = m_ptFocusReal.x;
	const double maxRealCoordX = m_log2real.x * m_maxLogicalCoord;		// 8192 pixels (~5 display widths) either side
	const char* pElem = (const char*)arElems;			// char* for 1 byte size! don't want 8 byte double

	const double minXReal = max(minX, xFocusReal - maxRealCoordX);
	const double maxXReal = min(maxX, xFocusReal + maxRealCoordX);

	int nSect = -1;				// -1 indicates before range, 0 in range, 1 after range
	idxPtsInit = 0;
	numPtsToDraw = 0;
	for (int i = 0; i < numPts; i++, pElem += sizeElem)
	{
		double xElem = *(double*)pElem;
		if (xElem >= minXReal && xElem <= maxXReal)				// point in range	
		{
			if (nSect == -1)		// first point in range
			{
				nSect = 0;
				idxPtsInit = i;
			}
			arPts[i].x = NEARINT(xReal2Log * (xElem - xFocusReal));
		}
		else					// don't use point - out of range!
		{
			if (nSect == -1)
			{
				i += 2;			// skip to next end node of bezier for test
				pElem += 2 * sizeElem;
			}
			else if (nSect == 0)		// section was in range, now out of range
				break;
		}
	}
	if (nSect == 0)
		numPtsToDraw = ((i - idxPtsInit - 1) / 3) * 3 + 1;		// move end point back to nearest bezier end node
}

/*
void CMoveView::SetLogicalPoints(CPoint* arPts, CVect2* arRealPts, int numPts)
{
	

}

void CMoveView::TransformPoint(CVector& pt)
{
	TrxNode = m_mxViewTransform * (pt - m_ptFocusReal);

	// check screen location
#define XLOC 0x03
#define YLOC 0x0c

}
*/

#if 0
void CMoveView::CheckLocation(CVect2& pt)	// line

	int loc = 0;
	if (pt.x < xMinWnd)
		loc |= 1;
	else if (pt.x > xMaxWnd)
		loc |= 2;
	if (pt.y < yMinWnd)
		loc |= 0x40;
	else if (pt.y > yMaxWnd)
		loc |= 0x80;

	if (loc1 == 0)			// a centre
	{
		prevSegIn = true;
		nextSegIn = true;
	}


	if ((loc1 & XLOC) == 0)		// must be a mid side - top/bottom
	{

	}
	else if ((loc1 & YLOC) == 0)	// must be a mid side - left/right
	{

	}
	else		// is a corner
	{

	}
-----
//	group of 2 points for line, 4 points for bezier
	if any in centre then IN (definately) if all in centre don't need to check magnitude
	all same cell then OUT (definately)
	all same row or col then OUT (definately)
	unless common row or col is centre then IN (definately - across centre)
	else IN (possibly - can check further)




}
#endif

void CMoveView::OnDraw(CDC* pDC)
{
/*
	screen size is assumed to be 1600x1200 @ 96dpi
	but for E790 monitor screen display size 360mm x 270mm (363mm x 272mm @ 112dpi)
	windows screen data is in registry under:
	HKEY_CURRENT_CONFIG - Display - Settings
*/
/*
	LPtoDP => DP = ViewportOrgDP + (LP - WindowOrgLP) * log2dev
	DPtoLP => LP = WindowOrgLP + (DP - ViewportOrgDP) * dev2log
*/

	if (m_nShowAxes == 1)
		DrawAxes(pDC);
}

void CMoveView::DrawAxes(CDC* pDC) 
{
/*
	int hsize = pDC->GetDeviceCaps(HORZSIZE);		// calc'd from res/DPIPhysical in registry
	int vsize = pDC->GetDeviceCaps(VERTSIZE);
	int hres = pDC->GetDeviceCaps(HORZRES);
	int vres = pDC->GetDeviceCaps(VERTRES);
	int logpixX = pDC->GetDeviceCaps(LOGPIXELSX);	// DPILogical in registry
	int logpixY = pDC->GetDeviceCaps(LOGPIXELSY);
*/

	pDC->SetWindowOrg(0,0);			// logical coords
	pDC->SetViewportOrg(0,0);		// to draw axis on window edge, device coords - changes when scrolling

	CVect2 ptLog16Orig;
	GetWindowULReal(ptLog16Orig);		// real coords of upper left view window - pixel(0,0)

	// Find visible axes region
	CRect rectClientLog;
	GetClientRect(rectClientLog);	// returns device units
	pDC->DPtoLP(rectClientLog);	// convert device points to logical points. Does NEARINT(dev2log * DP) with origions 0
											// upper left should still be (0,0) in logical coords!
	int markLog, markTextLog;		// logical coords
	double markLoc, markInc;		// real coords
	int iMarkCount, iMarkCountUpper;
	CPoint ptsLine[2];


	CFont font;
	CFont* pFontOrig = pDC->GetCurrentFont();
	LOGFONT lf;
	pFontOrig->GetLogFont(&lf);
	lf.lfHeight = 14;
	lf.lfWidth = 0;
	lf.lfWeight = 400;
	strcpy(lf.lfFaceName, "");
//	font.CreateFontIndirect(&lf);
	VERIFY(font.CreatePointFont(80, "MS Sans Serif"));
//	VERIFY(font.CreatePointFont(80, "Courier"));
	pDC->SelectObject(&font);


	char txt[20];
	pDC->SetTextColor(YELLOW);
	pDC->SetBkMode(TRANSPARENT);
	CPen penMark(PS_SOLID, 0, GREEN);
	CPen* pPenOrig = pDC->SelectObject(&penMark);

	// do X axis
	markTextLog = rectClientLog.bottom + (int)(4 * m_mm2log.y);	// 4mm up from bottom
	ptsLine[0].y = rectClientLog.bottom + (int)(6 * m_mm2log.y); 
	ptsLine[1].y = rectClientLog.bottom;
	markInc = FindMarkerInterval(m_Scale.x);
	iMarkCount = (int)floor(ptLog16Orig.x / markInc);			// lower count
	iMarkCountUpper = (int)ceil((ptLog16Orig.x + m_log2real.x * rectClientLog.Width()) / markInc);		// real units
	do
	{
		markLoc = iMarkCount++ * markInc;
		markLog = NEARINT((markLoc - ptLog16Orig.x) * m_real2log.x);
		sprintf(txt, "%g", markLoc);
		pDC->TextOut(markLog + 3, markTextLog, txt, strlen(txt));
		ptsLine[0].x = ptsLine[1].x = markLog;
		pDC->Polyline(ptsLine, 2);
	} while (iMarkCount <= iMarkCountUpper);		// do one outside window edge incase part of text is inside

	// do Y axis
	markTextLog = rectClientLog.left + (int)(1 * m_mm2log.x);	// 1mm in from left
	ptsLine[0].x = rectClientLog.left + (int)(6 * m_mm2log.x);
	ptsLine[1].x = rectClientLog.left;
	markInc = FindMarkerInterval(m_Scale.y);
	iMarkCount = (int)floor((ptLog16Orig.y - fabs(m_log2real.y * rectClientLog.Height())) / markInc);		// lower count
	iMarkCountUpper = (int)ceil(ptLog16Orig.y / markInc);
	do
	{
		markLoc = iMarkCount++ * markInc;
		markLog = NEARINT((markLoc - ptLog16Orig.y) * m_real2log.y);
		sprintf(txt, "%g", markLoc);
		
		ptsLine[0].y = ptsLine[1].y = markLog;
		pDC->Polyline(ptsLine, 2);			// if line is drawn last there can be a delay in actual painting
		pDC->TextOut(markTextLog, markLog, txt, strlen(txt));
	} while (iMarkCount <= iMarkCountUpper);		// do one outside window edge incase part of text is inside
	
	pDC->SelectObject(pPenOrig);
	pDC->SelectObject(pFontOrig);
}


double CMoveView::FindMarkerInterval(double scale)
{
// find marker interval at least 14mm apart. Will be 1, 2 or 5 to a power of 10
	double minReal = 14 / scale;		// 14mm is minimum marker distance
	ASSERT(scale != 0 && minReal != 0);
	double markInc = 1;
	while (minReal <= 1)
	{
		minReal *= 10;
		markInc *= 0.1;
	}
	while (minReal > 10)
	{
		minReal *= 0.1;
		markInc *= 10;
	}

	if (minReal <= 2)		//  1 < minReal <= 10
		markInc *= 2;
	else if (minReal <= 5)
		markInc *= 5;
	else
		markInc *= 10;
	return markInc;
}


void CMoveView::DrawAxesTriple(CDC* pDC)
{

	CPen penAxis[3];
	CBrush brAxesBall;
	penAxis[0].CreatePen(PS_SOLID, 10, RGB(255,255,0));
	penAxis[1].CreatePen(PS_SOLID, 10, RGB(255,160,0));
	penAxis[2].CreatePen(PS_SOLID, 10, RGB(255,0,0));
	brAxesBall.CreateSolidBrush(RGB(50,180,50));

	int radBall = 4;
	int lenAxes = 15;
	m_mxAxes.SetSize(6,3);
	m_mxAxes = 0;
	m_mxAxes.elem(0,0) = radBall;
	m_mxAxes.elem(1,0) = lenAxes;
	m_mxAxes.elem(2,1) = radBall;
	m_mxAxes.elem(3,1) = lenAxes;
	m_mxAxes.elem(4,2) = radBall;
	m_mxAxes.elem(5,2) = lenAxes;
	radBall = NEARINT(radBall * m_mm2log.x);		// convert to logical

	m_mxAxesX.Prod(m_mxAxes, m_mxViewRotateT);
	for (int i = 0; i < 6; i++)
	{
		m_mxAxesX.elem(i,0) *= m_mm2log.x;
		m_mxAxesX.elem(i,1) *= m_mm2log.y;
	}
	CPoint arAxesPts[6];
	ToCPoints(m_mxAxesX, arAxesPts, 6);

	CList<SzRef, SzRef> Zorder(4);
	SzRef zRef = {-1, 0};
	Zorder.AddHead(zRef);		// add -1 for ball

	for (int ax = 0; ax < 3; ax++)		// is axis is same as ball(0), put it behind!
	{
		zRef.ref = ax;
		zRef.z = m_mxAxesX.elem(2*ax+1, 2);

		POSITION ip = Zorder.GetHeadPosition();	// head is lowest z (further away)
		while (1)
			if (ip == NULL)		// at tail
			{
				Zorder.InsertAfter(ip, zRef);		// if null -> tail
				break;
			}
			else if (zRef.z <= Zorder.GetAt(ip).z)
			{
				Zorder.InsertBefore(ip, zRef);	// if null -> head
				break;
			}
			else
				Zorder.GetNext(ip);
	}

	CRect rectC;
	GetClientRect(&rectC);
	CPoint ptBL(rectC.left, rectC.bottom);
	pDC->SetViewportOrg(ptBL);
	CVect2 ptWO = CVect2(-15,-15) * m_mm2log;	// -15,-15 point in mm is in bottom left corner
	pDC->SetWindowOrg((CPoint)ptWO);

	// draw axis parts by Zorder list
	CPen* pPenOrig = pDC->GetCurrentPen();
	CBrush* pBrushOrig = pDC->GetCurrentBrush();
	while (!Zorder.IsEmpty())
	{
		int iAxisOb = Zorder.RemoveHead().ref;
		if (iAxisOb == -1)
		{
			pDC->SelectStockObject(NULL_PEN);
			pDC->SelectObject(brAxesBall);
			pDC->Ellipse(-radBall, -radBall, radBall, radBall);
		}
		else
		{
			pDC->SelectObject(penAxis[iAxisOb]);
			pDC->Polyline(arAxesPts + 2*iAxisOb, 2);
		}
	}
	pDC->SelectObject(pPenOrig);
	pDC->SelectObject(pBrushOrig);
}

void CMoveView::DrawRegionBox(CDC* pDC)
{
	if (!m_bShowRegionBox)
		return;

	CVector vt0, vtx, vty, vtz;
	vt0 = vtx = vty = vtz = m_ptMinRegion - m_ptFocusReal;
	vtx.x = m_ptMaxRegion.x - m_ptFocusReal.x;
	vty.y = m_ptMaxRegion.y - m_ptFocusReal.y;
	vtz.z = m_ptMaxRegion.z - m_ptFocusReal.z;

	CVect2 pt0, ptx, pty, ptz;
	pt0.ProdPart(m_mxViewTransform2Log, vt0);
	ptx.ProdPart(m_mxViewTransform2Log, vtx);
	pty.ProdPart(m_mxViewTransform2Log, vty);
	ptz.ProdPart(m_mxViewTransform2Log, vtz);

	// got 4 transformed points, calc other 4
	CPoint arPts[10];
	CVect2 ptyz = pty + ptz - pt0;
	arPts[0] = pt0;
	arPts[1] = pty;
	arPts[2] = ptyz;
	arPts[3] = ptz;
	arPts[4] = arPts[0]		;		// to draw back to start

	CVect2 ptxm0 = ptx - pt0;
	arPts[5] = ptx;
	arPts[6] = pty + ptxm0;
	arPts[7] = ptyz + ptxm0;
	arPts[8] = ptz + ptxm0;
	arPts[9] = arPts[5];				// to draw back to start

	SetViewportOrgToFocus(pDC);
	pDC->SetWindowOrg(0,0);

	CPen pen(PS_DASH, 0, LIGHTGRAY);
	CPen* pPenOrig = pDC->SelectObject(&pen);
	pDC->SetBkMode(TRANSPARENT);

	// draw 2 end rectangles
	pDC->Polyline(arPts + 0, 5);
	pDC->Polyline(arPts + 5, 5);
	// swap points 1&4 and 3&6 then draw 4 lines joining rectangles
	CPoint ptTmp;
	ptTmp = arPts[1]; arPts[1] = arPts[5]; arPts[5] = ptTmp;
	ptTmp = arPts[3]; arPts[3] = arPts[7]; arPts[7] = ptTmp;
	pDC->Polyline(arPts + 0, 2);
	pDC->Polyline(arPts + 2, 2);
	pDC->Polyline(arPts + 5, 2);
	pDC->Polyline(arPts + 7, 2);

	pDC->SelectObject(pPenOrig);
}

void CMoveView::SetTransformedRegion()		// call after transform is set
{
	CVector ptMinRegionRel = m_ptMinRegion - m_ptFocusReal;
	CVector ptMaxRegionRel = m_ptMaxRegion - m_ptFocusReal;
	// find min and max coords of transformed region corners
	for (int r = 0; r < 3; r++)		// to get x,y,z of min/max of transform
	{
		double min = 0, max = 0;
		for (int c = 0; c < 3; c++)		// to scan through x,y,z of min/max region
		{
			double coef = m_mxViewTransform.elem(r,c);
			if (coef < 0)		// select min or max depending on sign of m_mxViewTransform coefficient
			{
				min += coef * ptMaxRegionRel[c];
				max += coef * ptMinRegionRel[c];
			}
			else
			{
				min += coef * ptMinRegionRel[c];
				max += coef * ptMaxRegionRel[c];
			}
		}
		m_ptMinTxRegion[r] = min;
		m_ptMaxTxRegion[r] = max;
	}
}


/////////////////////////////////////////////////////////////////////////////
// CMoveView functions

void CMoveView::SetScale(CVect2& scale)
{
	m_Scale = scale;
	AdjustScalingFactors();
}

void CMoveView::SetScale(double xScale, double yScale)
{
	m_Scale.x = xScale;
	m_Scale.y = yScale;
	AdjustScalingFactors();
}

void CMoveView::SetScale(double scale)
{
	m_Scale.x = scale;
	m_Scale.y = scale;
	AdjustScalingFactors();
}

void CMoveView::ZoomBy(double factor)
{
	m_Scale = m_Scale * factor;
	AdjustScalingFactors();
}

void CMoveView::AdjustScalingFactors()		// does UpdateTotalViewSize()
{
	if (m_bIsoScaling)
		if (m_Scale.x > m_Scale.y)
			m_Scale.x = m_Scale.y;
		else
			m_Scale.y = m_Scale.x;
	// sets all other changable scale factors from m_Scale & m_mm2log
	m_real2log = m_Scale * m_mm2log;
	m_log2real = m_log2mm / m_Scale;
	m_real2dev = m_Scale * m_mm2dev;
	m_dev2real = m_dev2mm / m_Scale;

	m_mxViewTransform2Log = m_mxViewTransform;
	m_mxViewTransform2Log.ScaleRowsBy(m_real2log);

//	((CMainFrame*)AfxGetMainWnd())->ShowScale(m_Scale.x, m_Scale.y);
	AfxGetMainWnd()->SendMessage(WM_USER_SETSCALE, (WPARAM)&m_Scale.x, (LPARAM)&m_Scale.y);
	UpdateTotalViewSize();
}

void CMoveView::SetTransform()
{
	// make sure view angles are within +/-360 deg
	if (m_ViewAzimuth > 360) m_ViewAzimuth -= 360;
	else if (m_ViewAzimuth < -360) m_ViewAzimuth += 360;
	if (m_ViewElevation > 360) m_ViewElevation -= 360;
	else if (m_ViewElevation < -360) m_ViewElevation += 360;
	if (m_ViewOrientation > 360) m_ViewOrientation -= 360;
	else if (m_ViewOrientation < -360) m_ViewOrientation += 360;

	m_mxViewRotate.RotateAzEl(m_ViewAzimuth, m_ViewElevation);
	m_mxRotZ.RotateZ(m_ViewOrientation);
	m_mxViewRotate.Prod(m_mxRotZ, m_mxViewRotate);
	m_mxViewRotate.Transpose(m_mxViewRotateT);	// m_mxViewRotateT used by DrawAxesTriple()
	m_mxViewTransform = m_mxViewRotate;
	m_mxViewTransform.Inverse(m_mxViewTransformInv);

	m_mxViewTransform2Log = m_mxViewTransform;
	m_mxViewTransform2Log.ScaleRowsBy(m_real2log);

//	m_vtViewTxl = CVector(m_mxViewRotate * m_vtViewPreTxl) + m_vtViewPostTxl;

//	afxDump << "m_mxViewRotate:\n" << m_mxViewRotate;
//	afxDump << "m_vtViewTxl:\n" << m_vtViewTxl;
	SetTransformedRegion();
}

void CMoveView::SetViewRegion(const CVector& cnr1, const CVector& cnr2)
{
//	ptCnr1, ptCnr2 are the two real world coords of the region to be viewed
	m_ptMinRegion.Min(cnr1, cnr2);
	m_ptMaxRegion.Max(cnr1, cnr2);
	SetTransformedRegion();
}

void CMoveView::SetTotalViewSizes(CSize sizeTotal,
	const CSize& sizePage, const CSize& sizeLine)		// override of CScrollView function
{
	ASSERT(sizeTotal.cx >= 0 && sizeTotal.cy >= 0);
	m_totalLog = sizeTotal;

	// convert logical coordinate space to device coordinates (limit to max positive 32 bit signed int)
	double fSize;
//	TRACE2("CMoveView::SetTotalViewSizes to 0x%x, 0x%x\n", sizeTotal.cx, sizeTotal.cy);
	fSize = m_log2dev.x * m_totalLog.cx;
	m_totalDev.cx = (fSize >= INT_MAX) ? INT_MAX : NEARINT(fSize);
	fSize = fabs(m_log2dev.y) * m_totalLog.cy;
	m_totalDev.cy = (fSize >= INT_MAX) ? INT_MAX : NEARINT(fSize);

	m_pageDev.cx  = NEARINT(m_log2dev.x * sizePage.cx);
	m_pageDev.cy  = abs(NEARINT(m_log2dev.y * sizePage.cy));
	m_lineDev.cx  = NEARINT(m_log2dev.x * sizeLine.cx);
	m_lineDev.cy  = abs(NEARINT(m_log2dev.y * sizeLine.cy));

	m_bDevSizeOverflow = (m_totalDev.cx == INT_MAX || m_totalDev.cy == INT_MAX);

	// now adjust device specific sizes
	ASSERT(m_totalDev.cx >= 0 && m_totalDev.cy >= 0);
	CRect rectWnd;
	GetClientRect(rectWnd);
	if (m_pageDev.cx == 0)
		m_pageDev.cx = rectWnd.Width();		// page scroll is one window
	if (m_pageDev.cy == 0)
		m_pageDev.cy = rectWnd.Height();
	if (m_lineDev.cx == 0)
		m_lineDev.cx = m_pageDev.cx / 10;
	if (m_lineDev.cy == 0)
		m_lineDev.cy = m_pageDev.cy / 10;

//m_totalDev.cx = INT_MAX;
//m_totalDev.cy = INT_MAX;

	if (m_hWnd != NULL)
		UpdateBars();
}

void CMoveView::SetTotalViewSizes(const CVect2& sizeTotalDev)		// override of CScrollView function
{
	ASSERT(sizeTotalDev.x >= 0 && sizeTotalDev.y >= 0);	// ptSizeDev is Device units

	m_totalDevReal = sizeTotalDev;
	// limit device sizes to max positive 32 bit signed int
	m_totalDev.cx = (sizeTotalDev.x >= INT_MAX) ? INT_MAX : NEARINT(sizeTotalDev.x);
	m_totalDev.cy = (sizeTotalDev.y >= INT_MAX) ? INT_MAX : NEARINT(sizeTotalDev.y);
	m_bDevSizeOverflow = (m_totalDev.cx == INT_MAX || m_totalDev.cy == INT_MAX);

	m_totalLog.cx = NEARINT(m_dev2log.x * sizeTotalDev.x);
	m_totalLog.cy = NEARINT(fabs(m_dev2log.y) * sizeTotalDev.y);

	// now adjust device specific sizes
	ASSERT(m_totalDev.cx >= 0 && m_totalDev.cy >= 0);
	CRect rectWnd;
	GetClientRect(rectWnd);
	if (m_pageDev.cx == 0)
		m_pageDev.cx = rectWnd.Width();		// page scroll is one window
	if (m_pageDev.cy == 0)
		m_pageDev.cy = rectWnd.Height();
	if (m_lineDev.cx == 0)
		m_lineDev.cx = m_pageDev.cx / 10;
	if (m_lineDev.cy == 0)
		m_lineDev.cy = m_pageDev.cy / 10;

	if (m_hWnd != NULL)
		UpdateBars();
}

void CMoveView::SetFocusReal(const CVector& ptFocus)
{
	m_ptFocusReal = ptFocus;
	// contain focus point so that window is within the total view
	if (m_ptFocusReal.x < m_ptMinFocus.x) m_ptFocusReal.x = m_ptMinFocus.x;
	else if (m_ptFocusReal.x > m_ptMaxFocus.x) m_ptFocusReal.x = m_ptMaxFocus.x;
	if (m_ptFocusReal.y < m_ptMinFocus.y) m_ptFocusReal.y = m_ptMinFocus.y;
	else if (m_ptFocusReal.y > m_ptMaxFocus.y) m_ptFocusReal.y = m_ptMaxFocus.y;
// added for z 22/7/04
	if (m_ptFocusReal.z < m_ptMinFocus.z) m_ptFocusReal.z = m_ptMinFocus.z;
	else if (m_ptFocusReal.z > m_ptMaxFocus.z) m_ptFocusReal.z = m_ptMaxFocus.z;
}

//void CMoveView::SetFocusRealContain(const CVector& ptFocus)
//{
/* if focus point is limited by a focus min/max compensate by moving
	in screen z (depth) direction to achieve transverse focus movement.
	Until limited by all axis!
*/



/*
	CVector vtOutside = 0;
	if (ptFocus.x < m_ptMinFocus.x) vtOutside.x = ptFocus.x - m_ptMinFocus.x;
	else if (ptFocus.x > m_ptMaxFocus.x) vtOutside.x = ptFocus.x - m_ptMaxFocus.x;

	// bring ptFocus.x back to focus bounding box plane
	dzScreen = vtOutside[ax] / m_mxScreen2Real.elem(ax,3);
	vtCngFocus = dzScreen * m_mxScreen2Real(col(3));
*/



/*	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	m_ptWindowULReal.x = m_ptFocusReal.x - 0.5 * m_dev2real.x * rectWnd.Width();
	m_ptWindowULReal.y = m_ptFocusReal.y + 0.5 * fabs(m_dev2real.y) * rectWnd.Height();
*/
//}

/*void CMoveView::SetWindowULReal(const CVect2& ptWindowUL)
{
	m_ptWindowULReal = ptWindowUL;
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	m_ptFocusReal.x = m_ptWindowULReal.x + 0.5 * m_dev2real.x * rectWnd.Width();
	m_ptFocusReal.y = m_ptWindowULReal.y - 0.5 * fabs(m_dev2real.y) * rectWnd.Height();
}
*/

void CMoveView::GetWindowULReal(CVect2& ptWindowUL)
{
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	ptWindowUL.x = m_ptFocusReal.x - 0.5 * m_dev2real.x * rectWnd.Width();
	ptWindowUL.y = m_ptFocusReal.y + 0.5 * fabs(m_dev2real.y) * rectWnd.Height();
}

void CMoveView::GetWindowLRReal(CVect2& ptWindowLR)
{
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	ptWindowLR.x = m_ptFocusReal.x + 0.5 * m_dev2real.x * rectWnd.Width();
	ptWindowLR.y = m_ptFocusReal.y - 0.5 * fabs(m_dev2real.y) * rectWnd.Height();
}

void CMoveView::GetWindowSizeReal(CVect2& vtWindowSize)
{
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	vtWindowSize.x = m_dev2real.x * rectWnd.Width();
	vtWindowSize.y = fabs(m_dev2real.y) * rectWnd.Height();
}

CPoint CMoveView::GetFocusPointDev()
{
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	return rectWnd.CenterPoint();
}

void CMoveView::GetFocusPointDev(CPoint& ptFocusDev, CRect& rectWnd)
{
	GetClientRect(rectWnd);
	ptFocusDev = rectWnd.CenterPoint();
}

void CMoveView::SetViewportOrgToFocus(CDC* pDC)
{
	pDC->SetViewportOrg(GetFocusPointDev());
}

/*void CMoveView::GetWindowULLog(CPoint& ptWindowUL)
{
	ptWindowUL.x = NEARINT(m_real2log.x * (m_ptWindowULReal.x - m_ptViewULReal.x));
	ptWindowUL.y = NEARINT(m_real2log.y * (m_ptWindowULReal.y - m_ptViewULReal.y));
}
*/

/////////////////////////////////////////////////////
// Coordinate transform functions
/////////////////////////////////////////////////////

void CMoveView::WndToRealCoords(CVect2& ptReal, CPoint& ptWnd)
{
	CVect2 ptWndOrig;
	GetWindowULReal(ptWndOrig);		// or GetWindowOrg??
	ptReal = (m_log2real * ptWnd) + ptWndOrig;
}

void CMoveView::WndToRealCoords(CVect2* arPtReal, CPoint* arPtWnd, int num)
{
	CVect2 ptWndOrig;
	GetWindowULReal(ptWndOrig);		// or GetWindowOrg??
	for (int i = 0; i < num; i++)
		arPtReal[i] = (m_log2real * arPtWnd[i]) + ptWndOrig;
}

void CMoveView::RealToWndCoords(CPoint& ptWnd, CVect2& ptReal)
{
	CVect2 ptWndOrig;
	GetWindowULReal(ptWndOrig);		// or GetWindowOrg??
	ptWnd.x = NEARINT(m_real2dev.x * (ptReal.x - ptWndOrig.x));
	ptWnd.y = NEARINT(m_real2dev.y * (ptReal.y - ptWndOrig.y));
}

void CMoveView::RealToWndCoords(CPoint* arPtWnd, CVect2* arPtReal, int num)
{
	CVect2 ptWndOrig;
	GetWindowULReal(ptWndOrig);		// or GetWindowOrg??
	for (int i = 0; i < num; i++)
	{
		arPtWnd[i].x = NEARINT(m_real2log.x * (arPtReal[i].x - ptWndOrig.x));
		arPtWnd[i].y = NEARINT(m_real2log.y * (arPtReal[i].y - ptWndOrig.y));
	}
}








void CMoveView::FitRectToWindow(CRect& rectFitDev)	// device coords of rectangle (rel. to window) to fit to view
{
	CSize sizeFitDev = rectFitDev.Size();
	if (sizeFitDev.cx == 0)				// correct for any zero dimensions
	{	rectFitDev.InflateRect(1, 0); sizeFitDev.cx = 2; }
	if (sizeFitDev.cy == 0)
	{	rectFitDev.InflateRect(0, 1); sizeFitDev.cy = 2; }

	// find real coords of centre of rectangle for new focus
	CPoint ptFitCenterDev = rectFitDev.CenterPoint();
	CPoint ptFocusDev;
	CRect rectWnd;
	GetFocusPointDev(ptFocusDev, rectWnd);
	m_ptFocusReal += m_dev2real * (ptFitCenterDev - ptFocusDev);	// offset from focus centre
//	SetFocusReal(m_ptFocusReal);	// will adjust m_ptWindowULReal and any other reference points

	CSize sizeWnd = rectWnd.Size();
	// calc scaling by fitting rectDevice inside window (no border spacing)
	CVect2 scaleChange = (CVect2)sizeWnd / (CVect2)sizeFitDev;	// scale to fit window

	if (sizeFitDev.cx + sizeFitDev.cy < 10)		// small accidental rectangle!
	{
		double scaleMax = 2;
		if (scaleChange.x > scaleMax) scaleChange.x = scaleMax;
		if (scaleChange.y > scaleMax) scaleChange.y = scaleMax;
	}
	m_Scale = m_Scale * scaleChange;
	AdjustScalingFactors();
}

void CMoveView::FitViewRegionToWindow()
{
	CVector spanReal = m_ptMaxTxRegion - m_ptMinTxRegion;
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	CSize sizeWnd(rectWnd.Size());

	// calc scaling by fitting plots inside window with a border spacing
	if (spanReal.x != 0)
		m_Scale.x = (m_dev2mm.x * sizeWnd.cx - 2*m_ptViewBorder.x) / spanReal.x;		// scale to fit window minus two borders
	else
		m_Scale.x = 1;
	if (spanReal.y != 0)
		m_Scale.y = (fabs(m_dev2mm.y) * sizeWnd.cy - 2*m_ptViewBorder.y) / spanReal.y;		// scale to fit window minus two borders
	else
		m_Scale.y = 1;
	if (m_bIsoScaling)
		if (m_Scale.x > m_Scale.y)
			m_Scale.x = m_Scale.y;
		else
			m_Scale.y = m_Scale.x;

	AdjustScalingFactors();

	// centre focus on plot region
	m_ptFocusReal.x = (m_ptMinRegion.x + m_ptMaxRegion.x) / 2;
	m_ptFocusReal.y = (m_ptMinRegion.y + m_ptMaxRegion.y) / 2;
	m_ptFocusReal.z = (m_ptMinRegion.z + m_ptMaxRegion.z) / 2;
	
	UpdateTotalViewSize();
}

void CMoveView::UpdateTotalViewSize()
{		// called when scale (or window size) changes, changing scroll area
	CVector spanReal = m_ptMaxTxRegion - m_ptMinTxRegion;
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	CSize sizeWnd(rectWnd.Size());

	// set total view size to the graph region plus 1 window border all around
	// limit size values to max positive 32 bit signed int
	CSize sizeTotalViewLog;			// size in logical units
	double fSize;
	fSize = m_real2log.x * spanReal.x + m_dev2log.x * 2*sizeWnd.cx;	// = plot span plus 1 window border around plot
	sizeTotalViewLog.cx = (fSize >= INT_MAX) ? INT_MAX : NEARINT(fSize);
	fSize = fabs(m_real2log.y) * spanReal.y + fabs(m_dev2log.y) * 2*sizeWnd.cy;
	sizeTotalViewLog.cy = (fSize >= INT_MAX) ? INT_MAX : NEARINT(fSize);

	SetTotalViewSizes(sizeTotalViewLog);

/*	// calc real coords of top left of total view with current view size
	m_ptViewULReal.x = m_ptMinRegion.x - m_dev2real.x * sizeWnd.cx;	// for a full window border
	m_ptViewULReal.y = m_ptMaxRegion.y + fabs(m_dev2real.y) * sizeWnd.cy;
*/

	// calc real coords of focus limits with current scaling and window size
	CVector sizeFocusBorder;
	sizeFocusBorder.x = 0.5 * m_dev2real.x * sizeWnd.cx;			// for a half window border
	sizeFocusBorder.y = 0.5 * fabs(m_dev2real.y) * sizeWnd.cy;
	sizeFocusBorder.z = 0.0;
//	changed 22/7/04
//	m_ptMinFocus = m_ptMinTxRegion - sizeFocusBorder;
//	m_ptMaxFocus = m_ptMaxTxRegion + sizeFocusBorder;
	m_ptMinFocus = m_ptMinRegion;
	m_ptMaxFocus = m_ptMaxRegion;

	SetFocusReal(m_ptFocusReal);	// Checks focus is within focus limits
	ScrollFocusToWndCentre();		// move focus point to centre of window
	Invalidate();
}


void CMoveView::ScrollFocusToWndCentre()
{
	// calc scroll location to have focus point in centre of window
	CPoint ptWndULDev;
	if (!m_bDevSizeOverflow)	// Scroll size can be correctly represented  by a 32 bit int
	{
		CVect2 ptULLimitFocus;
		GetULLimitFocus(ptULLimitFocus);
		ptWndULDev.x = NEARINT(m_real2dev.x * (m_ptFocusReal.x - ptULLimitFocus.x));	// scroll offset in device units
		ptWndULDev.y = NEARINT(m_real2dev.y * (m_ptFocusReal.y - ptULLimitFocus.y));
	}
	else		// make device pos the correct ratio of the device scroll size
	{
		int iHmin, iHmax;
		int iVmin, iVmax;
		GetScrollRange(SB_HORZ, &iHmin, &iHmax);
		GetScrollRange(SB_VERT, &iVmin, &iVmax);
		int iHlim = GetScrollLimit(SB_HORZ);
		int iVlim = GetScrollLimit(SB_VERT);
		ptWndULDev.x = NEARINT((m_ptFocusReal.x - m_ptMinFocus.x) / (m_ptMaxFocus.x - m_ptMinFocus.x) * GetScrollLimit(SB_HORZ));
		ptWndULDev.y = NEARINT((m_ptMaxFocus.y - m_ptFocusReal.y) / (m_ptMaxFocus.y - m_ptMinFocus.y) * GetScrollLimit(SB_VERT));
	}
	ScrollToDevPosition(ptWndULDev);

/*	// Set m_ptFocusReal from actual scrolled position in case scrolling got to a limit
	ptWndULdev = GetDeviceScrollPosition();
	SetWindowULReal(m_dev2real * ptWndULDev + m_ptViewULReal);
*/
}

void CMoveView::ScrollToDevPosition(POINT pt)    // device coordinates
{
	// in device coordinates - limit if out of range
	int xMax = GetScrollLimit(SB_HORZ);
	int yMax = GetScrollLimit(SB_VERT);
//	ASSERT(pt.x >= 0 && pt.x <= xMax);
//	ASSERT(pt.y >= 0 && pt.y <= yMax);
	if (pt.x < 0) pt.x = 0;
	else if (pt.x > xMax) pt.x = xMax;
	if (pt.y < 0) pt.y = 0;
	else if (pt.y > yMax) pt.y = yMax;

//	CScrollView::ScrollToDevicePosition(pt);	// following code is from CScrollView::ScrollToDevicePosition, less ScrollWindow()
	SetScrollPos(SB_HORZ, pt.x);
	SetScrollPos(SB_VERT, pt.y);
}









/////////////////////////////////////////////////////////////////////////////
// CMoveView message handlers

BOOL CMoveView::OnEraseBkgnd(CDC* pDC) 
{
	CRect rect;
	GetClientRect(&rect);
	pDC->FillSolidRect(&rect, BLACK);
	return true;           // Erased
//	return CScrollView::OnEraseBkgnd(pDC);
}


void CMoveView::OnButtonDown(UINT nFlags, CPoint point) 
{
	switch (nFlags & (MK_LBUTTON | MK_MBUTTON | MK_RBUTTON))
	{
		case 0:			  m_nDynMoveType = DYNA_NONE; break;
		case MK_LBUTTON: m_nDynMoveType = DYNA_TRANSLATE; m_hCursor = m_hCursorMove; break;
		case MK_RBUTTON: m_nDynMoveType = DYNA_ZOOMWINDOW; break;
		case MK_MBUTTON: m_nDynMoveType = DYNA_DIRECTION; m_hCursor = m_hCursorDir; break;
//		case MK_MBUTTON: m_nDynMoveType = DYNA_ZOOMROTATE; break;
		case MK_LBUTTON | MK_RBUTTON: m_nDynMoveType = DYNA_ZOOMWINDOW; break;
	}
	if (m_nDynMoveType != DYNA_NONE)
		SetCapture();
	if (m_hCursor)
		SetCursor(m_hCursor);
	m_ptButtonDown = point;
	m_ptPrevMouseMove = point;
	m_ptFocusRealInit = m_ptFocusReal;
	m_bFirstMouseMove = true;

}

void CMoveView::OnButtonUp(UINT nFlags, CPoint point) 
{
	CPoint ptChange = point - m_ptPrevMouseMove;
	CPoint ptRelative = point - m_ptButtonDown;
	if (nFlags & (MK_LBUTTON | MK_MBUTTON | MK_RBUTTON) != 0)	// if some button still down
		return;
	switch (m_nDynMoveType)
	{
		case DYNA_NONE: return;
		case DYNA_TRANSLATE: break;
		case DYNA_DIRECTION: break;
		case DYNA_ZOOMROTATE: break;
		case DYNA_ZOOMWINDOW:
		{	// if drag rectangle needs deleting
			CClientDC dc(this);
			CSize sizeDrag(3,3);
			CRect rectDrag(0,0,0,0);	// draws none
			CRect rectLast(m_ptButtonDown, m_ptPrevMouseMove);
			rectLast.NormalizeRect();
			dc.DrawDragRect(&rectDrag, sizeDrag, rectLast, sizeDrag);
			FitRectToWindow(rectLast);		// does UpdateTotalViewSize()
			break;
		}
	}
	ReleaseCapture();
	m_hCursor = NULL;
	SetCursor(m_hCursor);
	m_nDynMoveType = DYNA_NONE;
	m_ptButtonUp = point;
}

void CMoveView::OnMouseMove(UINT /*nFlags*/, CPoint point) 
{
	CPoint ptChange = point - m_ptPrevMouseMove;
	CPoint ptRelative = point - m_ptButtonDown;
	switch (m_nDynMoveType)
	{
		case DYNA_NONE:		// send location message to status bar if enabled
		{
			CVect2 ptLog16Orig;
			GetWindowULReal(ptLog16Orig);
			CVect2 realLoc = (m_log2real * point) + ptLog16Orig;
			AfxGetMainWnd()->SendMessage(WMU_SETPOSITIONINDICATOR, 2, (LPARAM)&realLoc);
		}	break;
		case DYNA_TRANSLATE: // move view in device units of point - m_ptButtonDown
			if (m_b3DView)
			{
				CVector ptOffset;		// can't use CVector constructor but can use =
				ptOffset.ProdPart(m_mxViewTransformInv, CVect2(ptRelative) * m_dev2real * m_DynaTranslateAmp);
				SetFocusReal(m_ptFocusRealInit - ptOffset);
			}
			else
				SetFocusReal(m_ptFocusRealInit - CVector(CVect2(ptRelative) * m_dev2real * m_DynaTranslateAmp));
			ScrollFocusToWndCentre();
			Invalidate();
			break;
		case DYNA_DIRECTION:
			m_ViewAzimuth -= ptChange.x / 8.0;		// rotate about window plane Y axis - rotates about world z
			m_ViewElevation += ptChange.y / 8.0;	// rotate about window plane X axis - does already
			SetTransform();
			UpdateTotalViewSize();		// Invalidate()'s - changes m_ptMinTxRegion, m_ptMaxTxRegion
			break;
		case DYNA_ZOOMROTATE:
			m_ViewOrientation += ptChange.x / 8.0;		// rotate about window plane Z axis
			SetTransform();
			ZoomBy(exp(-ptChange.y / 100.0));	// does UpdateTotalViewSize()
			break;
		case DYNA_ZOOMWINDOW:
		{
			CClientDC dc(this);
			CSize sizeDrag(3,3);
			CRect rectDrag(m_ptButtonDown, point);
			CRect rectLast(m_ptButtonDown, m_ptPrevMouseMove);
			rectDrag.NormalizeRect();
			rectLast.NormalizeRect();
			CRect *pRectLast = m_bFirstMouseMove ? NULL : &rectLast;
			dc.DrawDragRect(&rectDrag, sizeDrag, pRectLast, sizeDrag);
		}
			break;
//			HCURSOR hCursor = AfxGetApp()->LoadCursor(IDC_CURSORZOOMIN);
//			HCURSOR hCursorOrig = SetCursor(hCursor);
		default:
			break;
	}
	m_ptPrevMouseMove = point;
	m_bFirstMouseMove = false;
}

BOOL CMoveView::OnMouseWheel(UINT /*nFlags*/, short zDelta, CPoint /*pt*/) 
{
//	return CScrollView::OnMouseWheel(nFlags, zDelta, pt);
	// Use mouse wheel to zoom in or out.  WHEEL_DELTA = 120
	ZoomBy(exp(zDelta * 2e-3));		// does UpdateTotalViewSize()
	return true;
}


int CMoveView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CScrollView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	// TODO: Add your specialized creation code here
	
	return 0;
}


BOOL CMoveView::OnScrollBy(CSize sizeScroll, BOOL bDoScroll) 
{
	int bScrolled = CScrollView::OnScrollBy(sizeScroll, bDoScroll);
	if (bDoScroll && bScrolled)		// move centre focus
	{
		if (m_b3DView)
		{
			CVector ptOffset;		// can't use CVector constructor but can use =
			ptOffset.ProdPart(m_mxViewTransformInv, CVect2(sizeScroll) * m_dev2real);
			SetFocusReal(m_ptFocusReal + ptOffset);
		}
		else
			SetFocusReal(m_ptFocusReal + CVector(CVect2(sizeScroll) * m_dev2real));
		Invalidate();
	}
	return bScrolled;
}

void CMoveView::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	if (nSBCode == SB_THUMBTRACK || nSBCode == SB_THUMBPOSITION)
	{
		SCROLLINFO si;
		GetScrollInfo(SB_HORZ, &si, SIF_TRACKPOS);
		nPos = si.nTrackPos;		// nPos is only 16 bits worth in WM_HSCROLL
	}
	CScrollView::OnHScroll(nSBCode, nPos, pScrollBar);
}

void CMoveView::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	if (nSBCode == SB_THUMBTRACK || nSBCode == SB_THUMBPOSITION)
	{
		SCROLLINFO si;
		GetScrollInfo(SB_VERT, &si, SIF_TRACKPOS);
		nPos = si.nTrackPos;		// nPos is only 16 bits worth in WM_VSCROLL
	}
	CScrollView::OnVScroll(nSBCode, nPos, pScrollBar);
}

	
void CMoveView::OnViewZoomAll() 
{
	FitViewRegionToWindow();
}

void CMoveView::OnViewZoom11() 
{
	SetScale(1);
}

void CMoveView::OnViewZoomIn() 
{
	ZoomBy(2);
}

void CMoveView::OnViewZoomOut() 
{
	ZoomBy(0.5);
}

void CMoveView::OnViewIsoScale()
{
	m_bIsoScaling = !m_bIsoScaling;
	AdjustScalingFactors();
}

void CMoveView::OnUpdateViewIsoScale(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(m_bIsoScaling);	
}

BOOL CMoveView::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message) 
{
	if (nHitTest == HTCLIENT && m_hCursor)
	{
			SetCursor(m_hCursor);
			return true;
	}
	return CScrollView::OnSetCursor(pWnd, nHitTest, message);
}

void CMoveView::OnSetFocus(CWnd* pOldWnd)
{
	CScrollView::OnSetFocus(pOldWnd);
	
	// TODO: Add your message handler code here
	
}

void CMoveView::OnKillFocus(CWnd* pNewWnd)
{
	CScrollView::OnKillFocus(pNewWnd);
	
	// TODO: Add your message handler code here
	AfxGetMainWnd()->SendMessage(WMU_SETPOSITIONINDICATOR, 0, 0);	// to blank it!	
}

void CMoveView::OnViewStoreLocation()
{
	// store view focus and scale to static members
	GetFocusReal(ms_ptFocusRealStored);
	GetScale(ms_ScaleStored);
}

void CMoveView::OnViewSetLocation()
{
	if (ms_ScaleStored.x == 0)		// view location not stored yet!
		return;
	// get view focus and scale from static members
	SetFocusReal(ms_ptFocusRealStored);
	SetScale(ms_ScaleStored);
	UpdateTotalViewSize();
}

void CMoveView::OnUpdateViewSetLocation(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ms_ScaleStored.x != 0);		// view location not stored yet!
}

void CMoveView::OnViewPrevLocation()
{
}
void CMoveView::OnViewNextLocation()
{
}

void CMoveView::OnTieView()		// add this view to group of tied views that move together
{

}

void CMoveView::OnViewDirSetToX()
{
	m_ViewAzimuth = 90;
	m_ViewElevation = 0;
	m_ViewOrientation = 0;
	SetTransform();
	Invalidate();

/*	m_mxViewTransform = 0;
	m_mxViewTransform.elem(0,1) = 1;
	m_mxViewTransform.elem(1,2) = 1;
	m_mxViewTransform.elem(2,0) = 1;
	m_mxViewRotate = m_mxViewTransform;
	m_mxViewRotate.Transpose(m_mxViewRotateT);
	SetViewFrom(CVector(1,0,0));
*/
}
void CMoveView::OnViewDirSetToY()
{
	m_ViewAzimuth = 0;
	m_ViewElevation = 0;
	m_ViewOrientation = 0;
	SetTransform();
	Invalidate();

/*	m_mxViewTransform = 0;
	m_mxViewTransform.elem(0,0) = 1;
	m_mxViewTransform.elem(1,2) = 1;
	m_mxViewTransform.elem(2,1) = -1;
	m_mxViewRotate = m_mxViewTransform;
	m_mxViewRotate.Transpose(m_mxViewRotateT);
	SetViewFrom(CVector(0,1,0));
*/
}
void CMoveView::OnViewDirSetToZ()
{
	m_ViewAzimuth = 0;
	m_ViewElevation = 90;
	m_ViewOrientation = 0;
	SetTransform();
	Invalidate();

/*	m_mxViewTransform.Identity();
	m_mxViewRotate = m_mxViewTransform;
	m_mxViewRotate.Transpose(m_mxViewRotateT);
	SetViewFrom(CVector(0,0,1));
*/
}
void CMoveView::OnViewDirSetToIso()
{
	m_ViewAzimuth = 30;
	m_ViewElevation = 30;
	m_ViewOrientation = 0;
	SetTransform();
	Invalidate();

//	SetViewFrom(CVector(1,-1,1));
}
void CMoveView::OnViewDirFlip()		// views from other side - rotate 180 about screen Y axis
{
	m_ViewAzimuth += 180;
	m_ViewElevation = -m_ViewElevation;
	SetTransform();
	Invalidate();
}
void CMoveView::SetViewFrom(const CVector& /*vt*/)
{

	Invalidate();
}

void CMoveView::OnViewRegionBox() 
{
	m_bShowRegionBox = !m_bShowRegionBox;
	Invalidate();
}

void CMoveView::OnUpdateViewRegionBox(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck(m_bShowRegionBox);	
}

void CMoveView::OnOptionsScrollBars() 
{
	if (GetStyle() & (WS_HSCROLL | WS_VSCROLL))
		ModifyStyle(WS_HSCROLL | WS_VSCROLL, 0);		// remove scroll bars
	else
		ModifyStyle(0, WS_HSCROLL | WS_VSCROLL);		// remove scroll bars
}

void CMoveView::OnUpdateOptionsScrollBars(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck(GetStyle() & (WS_HSCROLL | WS_VSCROLL));	
}
