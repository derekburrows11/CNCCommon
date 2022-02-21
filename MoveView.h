// MoveView.h : header file
//


#if !defined(AFX_MOVEVIEW_H__DA6639E5_42FD_11D6_86C3_9B97A5B31D24__INCLUDED_)
#define AFX_MOVEVIEW_H__DA6639E5_42FD_11D6_86C3_9B97A5B31D24__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


#include "Matrix.h"

#define WM_USER_SETSCALE	WM_USER+2


struct SDrawNodes
{
	int iRadius;
	int nEndType;		// 0 none, 1 point, 2 circle, 3 open circle, 4 closed square, 5 open square, 6 x, 7 +,
	int nControlType;
	int iPeriod;
	int iStart;
	CPoint* arPts;
	int numPts;
};

void DrawNodes(CDC* pDC, const SDrawNodes& dn);


/////////////////////////////////////////////////////////////////////////////
// CMoveView view

class CMoveView : public CScrollView
{
protected:
	struct SzRef { int ref; double z; };
	
	CMoveView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CMoveView)


// Operations
public:

	void SetViewRegion(const CVector& cnr1, const CVector& cnr2);
	void SetTotalViewSizes(CSize sizeTotal,
				const CSize& sizePage = sizeDefault,
				const CSize& sizeLine = sizeDefault);			// override of CScrollView function
	void SetTotalViewSizes(const CVect2& sizeTotalDev);	// override of CScrollView function
	void SetScale(CVect2& scale);
	void SetScale(double xScale, double yScale);
	void SetScale(double scale);		// sets x & y to same scale
	void GetScale(CVect2& scale) { scale = m_Scale; }
	void ZoomBy(double factor);
	void SetTransform();
	void SetTransformedRegion();		// call after transform is set
	void SetViewFrom(const CVector& vt);

	void SetFocusReal(const CVector& ptFocus);
	void GetFocusReal(CVector& ptFocus) { ptFocus = m_ptFocusReal; }
	void GetFocusReal(CVect2& ptFocus) { ptFocus = m_ptFocusReal; }
	CPoint GetFocusPointDev();
	void GetFocusPointDev(CPoint& ptFocusDev, CRect& rectWnd);
	void SetViewportOrgToFocus(CDC* pDC);
	void GetULLimitFocus(CVect2& pt) { pt.x = m_ptMinFocus.x; pt.y = m_ptMaxFocus.y; }
	void GetWindowULReal(CVect2& ptWindowUL);
	void GetWindowLRReal(CVect2& ptWindowLR);
	void GetWindowSizeReal(CVect2& vtWindowSize);
//	void SetWindowULReal(const CVect2& ptWindowUL);
//	void GetWindowULLog(CPoint& ptWindowUL);

	void ScrollToDevPosition(POINT pt);		// overridden - device coordinates


// Coordinate conversion
	void WndToRealCoords(CVect2& ptReal, CPoint& ptWnd);
	void WndToRealCoords(CVect2* arPtReal, CPoint* arPtWnd, int num);
	void RealToWndCoords(CPoint& ptWnd, CVect2& ptReal);
	void RealToWndCoords(CPoint* arPtWnd, CVect2* arPtReal, int num);


	void SetLogicalXValues(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 int& idxPtsInit, int& numPtsDraw);
	void SetLogicalXValuesInRange(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 double minX, double maxX, int& idxPtsInit, int& numPtsToDraw);
	void SetLogicalXBezValues(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 int& idxPtsInit, int& numPtsDraw);
	void SetLogicalXBezValuesInRange(CPoint* arPts, const double* arElems, int sizeElem, int numPts,
											 double minX, double maxX, int& idxPtsInit, int& numPtsToDraw);



// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMoveView)
	public:
	virtual void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = NULL);
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	virtual void OnInitialUpdate();     // first time after construct
	virtual BOOL OnScrollBy(CSize sizeScroll, BOOL bDoScroll = TRUE);
	//}}AFX_VIRTUAL


// Attributes
protected:

	// Scaling x & y components
	CVect2 m_Scale;			// could be called m_real2mm
	CVect2 m_real2log;		// to convert graph data to logical units. m_Scale * m_mm2log
	CVect2 m_log2real;
	CVect2 m_real2dev;
	CVect2 m_dev2real;

	// fixed scaling variables
	CVect2 m_mm2log;		// to convert screen mm to logical units (0.1mm) * screen correction factor 112/96 (112dpi instead of 96dpi as windows expects!)
	CVect2 m_log2mm;
	CVect2 m_mm2dev;		// 112dpi / 25.4mmpi
	CVect2 m_dev2mm;		// 25.4mmpi / 112dpi
	CVect2 m_log2dev;		// log2mm * mm2dev
	CVect2 m_dev2log;		// dev2mm * mm2log

	double m_maxLogicalCoord;	// max logical coord with 16 bit values for win98 line functions

	// Reference points
	CVector m_ptFocusReal;			// focus point in real coords
//	CVect2 m_ptWindowULReal;		// previously 'm_ptLog16OrigReal' real point of Upper Left of client window (Log16 origion for DC functions)
//	CVect2 m_ptViewULReal;		// Coord's of upper left of total view - used for scrolling location

	// 3D viewing rotations
	CMatrix m_mxRotZ;					// used during calculations
	CMatrix m_mxViewRotate;			// 3 x 3 to rotate real point to view pos, scale = 1
	CMatrix m_mxViewRotateT;		// transpose of Rotate
	CMatrix m_mxViewTransform;			// = Rotate scaling
	CMatrix m_mxViewTransformInv;		// = inverse of Transform
	CMatrix m_mxViewTransform2Log;	// mxTransform x real2log - convert real point to logical screen pos
	double m_ViewAzimuth;
	double m_ViewElevation;
	double m_ViewOrientation;

	// Used during calculations for axis triple
	CMatrix m_mxAxes;
	CMatrix m_mxAxesX;

	// View location memorys
	static CVector ms_ptFocusRealStored;
	static CVect2 ms_ScaleStored;

	// View region bounds
	CVector m_ptMinRegion;		// Actual plot region bound
	CVector m_ptMaxRegion;
	CVector m_ptMinFocus;		// Plot Region plus a border which will limit focus region
	CVector m_ptMaxFocus;
	CVector m_ptMinTxRegion;	// Regions limits after transform
	CVector m_ptMaxTxRegion;

	CPoint m_ptViewBorder;		// mm plot border spacing for fit to window

	// Scrolling flag for large view in device coords (when large zoom is used)
	bool m_bDevSizeOverflow;
	CVect2 m_totalDevReal;	// total size in device units but stored as double

	// dynamic mouse movement variables
	CPoint m_ptButtonDown, m_ptButtonUp, m_ptPrevMouseMove;
	bool m_bFirstMouseMove;
	enum
	{	DYNA_NONE = 0,
		DYNA_TRANSLATE,
		DYNA_DIRECTION,
		DYNA_ZOOMROTATE,
		DYNA_ZOOMWINDOW,
	} m_nDynMoveType;
	CVector m_ptFocusRealInit;

	double m_DynaTranslateAmp;		// for amplified translations (default = 1.0)


	HCURSOR m_hCursor;	// to use different cursors in different modes
	HCURSOR m_hCursorMove;
	HCURSOR m_hCursorDir;

	// View options
	int m_nShowAxes;
	bool m_bIsoScaling;
	bool m_b3DView;					// true to use 3D view transforming
	bool m_bShowRegionBox;

// Implementation
protected:
	virtual ~CMoveView();

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif


	void DrawAxes(CDC* pDC);
	void DrawAxesTriple(CDC* pDC);
	void DrawRegionBox(CDC* pDC);
	double FindMarkerInterval(double scale);

	void AdjustScalingFactors();
	void FitRectToWindow(CRect& rectFitDev);
	void FitViewRegionToWindow();
	void UpdateTotalViewSize();
	void ScrollFocusToWndCentre();


	void OnButtonDown(UINT nFlags, CPoint point);
	void OnButtonUp(UINT nFlags, CPoint point);




	// Generated message map functions
protected:
	//{{AFX_MSG(CMoveView)
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point) { OnButtonDown(nFlags, point); }
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point) { OnButtonDown(nFlags, point); }
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point) { OnButtonDown(nFlags, point); }
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point) { OnButtonUp(nFlags, point); }
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point) { OnButtonUp(nFlags, point); }
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point) { OnButtonUp(nFlags, point); }
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnViewZoom11();
	afx_msg void OnViewZoomIn();
	afx_msg void OnViewZoomOut();
	afx_msg void OnViewZoomAll();
	afx_msg void OnViewIsoScale();
	afx_msg void OnUpdateViewIsoScale(CCmdUI* pCmdUI);
	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
	afx_msg void OnSetFocus(CWnd* pOldWnd);
	afx_msg void OnKillFocus(CWnd* pNewWnd);
	afx_msg void OnViewStoreLocation();
	afx_msg void OnViewSetLocation();
	afx_msg void OnUpdateViewSetLocation(CCmdUI* pCmdUI);
	afx_msg void OnViewDirSetToX();
	afx_msg void OnViewDirSetToY();
	afx_msg void OnViewDirSetToZ();
	afx_msg void OnViewDirSetToIso();
	afx_msg void OnViewDirFlip();
	afx_msg void OnViewRegionBox();
	afx_msg void OnUpdateViewRegionBox(CCmdUI* pCmdUI);
	afx_msg void OnOptionsScrollBars();
	afx_msg void OnUpdateOptionsScrollBars(CCmdUI* pCmdUI);
	//}}AFX_MSG
	void OnViewPrevLocation();
	void OnViewNextLocation();
	void OnTieView();
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MOVEVIEW_H__DA6639E5_42FD_11D6_86C3_9B97A5B31D24__INCLUDED_)
