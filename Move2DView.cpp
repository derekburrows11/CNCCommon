// Move2DView.cpp : implementation file
//

#include "stdafx.h"

#include "Move2DView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMove2DView

IMPLEMENT_DYNCREATE(CMove2DView, CView)

CMove2DView::CMove2DView()
{
	m_hCursor = NULL;
	m_hCursorMove = NULL;
	m_hCursorDir = NULL;

	// constant scaling factors
	m_mm2dev = 112.0 / 25.4;		// 112dpi / 25.4mmpi			19" screen @1600x1200 is 112dpi
	m_mm2dev.y = -m_mm2dev.y;		// device y is pos down
	m_dev2mm = 25.4 / 112.0;		// 25.4mmpi / 112dpi
	m_dev2mm.y = -m_dev2mm.y;		// device y is pos down


	m_bIsoScaling = true;
	m_b3DView = false;
	m_bShowRegionBox = false;

	// dynamic movement variables
	m_nDynMoveType = DYNA_NONE;
	m_DynaTranslateAmp = 1.0;		// for amplified translations (default = 1.0)


}

CMove2DView::~CMove2DView()
{
}


BEGIN_MESSAGE_MAP(CMove2DView, CView)
	//{{AFX_MSG_MAP(CMove2DView)
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_LBUTTONDOWN()
	ON_WM_MBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MBUTTONUP()
	ON_WM_RBUTTONUP()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////

void CMove2DView::OnInitialUpdate() 
{

	m_ptFocusReal = 0;
//	m_ptMinRegion = 0;
//	m_ptMaxRegion = 0;
	SetScale(1);			// default scale

//	m_ptViewBorder.x = 10;	// mm
//	m_ptViewBorder.y = 10;	// mm

	
	CView::OnInitialUpdate();
}

/////////////////////////////////////////////////////////////////////////////
// CMove2DView drawing

void CMove2DView::OnDraw(CDC* /*pDC*/)
{
//	CDocument* pDoc = GetDocument();

}

/////////////////////////////////////////////////////////////////////////////
// CMove2DView diagnostics

#ifdef _DEBUG
void CMove2DView::AssertValid() const
{
	CView::AssertValid();
}

void CMove2DView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMove2DView message handlers

void CMove2DView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	m_nDynMoveType = DYNA_TRANSLATE;
	m_hCursor = m_hCursorMove;
	OnButtonDown(nFlags, point);
}
void CMove2DView::OnMButtonDown(UINT nFlags, CPoint point) 
{
	m_nDynMoveType = DYNA_NONE;
	m_hCursor = NULL;
	OnButtonDown(nFlags, point);
}
void CMove2DView::OnRButtonDown(UINT nFlags, CPoint point) 
{
	m_nDynMoveType = DYNA_NONE;
	m_hCursor = NULL;
	OnButtonDown(nFlags, point);
}

void CMove2DView::OnButtonDown(UINT nFlags, CPoint point) 
{
	if (m_nDynMoveType != DYNA_NONE)
		SetCapture();
	if (m_hCursor)
		SetCursor(m_hCursor);
	m_ptButtonDown = point;
	m_ptPrevMouseMove = point;
	m_ptFocusRealInit = m_ptFocusReal;
	m_bFirstMouseMove = true;
}

void CMove2DView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	OnButtonUp(nFlags, point);
}
void CMove2DView::OnMButtonUp(UINT nFlags, CPoint point) 
{
	OnButtonUp(nFlags, point);
}
void CMove2DView::OnRButtonUp(UINT nFlags, CPoint point) 
{
	OnButtonUp(nFlags, point);
}

void CMove2DView::OnButtonUp(UINT nFlags, CPoint point) 
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

void CMove2DView::OnMouseMove(UINT /*nFlags*/, CPoint point) 
{
	CPoint ptChange = point - m_ptPrevMouseMove;
	CPoint ptRelative = point - m_ptButtonDown;
	switch (m_nDynMoveType)
	{
		case DYNA_NONE:		// send location message to status bar if enabled
		{
/*			CVect2 ptLog16Orig;
			GetWindowULReal(ptLog16Orig);
			CVect2 realLoc = (m_log2real * point) + ptLog16Orig;
			AfxGetMainWnd()->SendMessage(WMU_SETPOSITIONINDICATOR, 2, (LPARAM)&realLoc);
*/		}	break;
		case DYNA_TRANSLATE: // move view in device units of point - m_ptButtonDown
			SetFocusReal(m_ptFocusRealInit - CVector2(ptRelative) * m_dev2real * m_DynaTranslateAmp);
			//ScrollFocusToWndCentre();
			Invalidate();
			break;
		case DYNA_DIRECTION:
/*			m_ViewAzimuth -= ptChange.x / 8.0;		// rotate about window plane Y axis - rotates about world z
			m_ViewElevation += ptChange.y / 8.0;	// rotate about window plane X axis - does already
			SetTransform();
			UpdateTotalViewSize();		// Invalidate()'s - changes m_ptMinTxRegion, m_ptMaxTxRegion
*/			break;
		case DYNA_ZOOMROTATE:
			//m_ViewOrientation += ptChange.x / 8.0;		// rotate about window plane Z axis
			//SetTransform();
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

BOOL CMove2DView::OnMouseWheel(UINT /*nFlags*/, short zDelta, CPoint /*pt*/) 
{
//	return CScrollView::OnMouseWheel(nFlags, zDelta, pt);
	// Use mouse wheel to zoom in or out.  WHEEL_DELTA = 120
	ZoomBy(exp(zDelta * 2e-3));		// does UpdateTotalViewSize()
	return true;		// return nonzero if mouse wheel scrolling is enabled; otherwise 0.
}

//////////////////////////////////////////////////////////////////////////
// Focus point functions

void CMove2DView::SetFocusReal(const CVector2& ptFocus)
{
	m_ptFocusReal = ptFocus;
}

CPoint CMove2DView::GetFocusPointDev()
{
	CRect rectWnd;
	GetClientRect(rectWnd);		// in device coords (pixels)
	return rectWnd.CenterPoint();
}

void CMove2DView::GetFocusPointDev(CPoint& ptFocusDev, CRect& rectWnd)
{
	GetClientRect(rectWnd);
	ptFocusDev = rectWnd.CenterPoint();
}


//////////////////////////////////////////////////////////////////////////
// Scale functions

void CMove2DView::SetScale(double xScale, double yScale)
{
	m_Scale.x = xScale;
	m_Scale.y = yScale;
	m_Scale.Set(xScale, yScale);
	AdjustScalingFactors();
}


void CMove2DView::ZoomBy(double factor)
{
	m_Scale = m_Scale * factor;
	AdjustScalingFactors();
}

void CMove2DView::AdjustScalingFactors()		// does UpdateTotalViewSize()
{
	if (m_bIsoScaling)
		if (m_Scale.x > m_Scale.y)
			m_Scale.x = m_Scale.y;
		else
			m_Scale.y = m_Scale.x;
	// sets all other changable scale factors from m_Scale & m_mm2log
//	m_real2log = m_Scale * m_mm2log;
//	m_log2real = m_log2mm / m_Scale;
	m_real2dev = m_Scale * m_mm2dev;
	m_dev2real = m_dev2mm / m_Scale;

	//UpdateTotalViewSize();
	Invalidate();
}

void CMove2DView::FitRectToWindow(CRect& rectFitDev)	// device coords of rectangle (rel. to window) to fit to view
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


