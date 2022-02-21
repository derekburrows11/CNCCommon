#if !defined(AFX_MOVE2DVIEW_H__80D04CC3_6589_4373_9D3F_0B9DB8424200__INCLUDED_)
#define AFX_MOVE2DVIEW_H__80D04CC3_6589_4373_9D3F_0B9DB8424200__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Move2DView.h : header file
//


#include "Vector.h"


/////////////////////////////////////////////////////////////////////////////
// CMove2DView view

class CMove2DView : public CView
{
protected:
	CMove2DView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CMove2DView)

// Attributes
protected:

	// Scaling x & y components
	CVector2 m_Scale;			// could be called m_real2mm
//	CVector2 m_real2log;		// to convert graph data to logical units. m_Scale * m_mm2log
//	CVector2 m_log2real;
	CVector2 m_real2dev;
	CVector2 m_dev2real;

	// fixed scaling variables
//	CVect2 m_mm2log;		// to convert screen mm to logical units (0.1mm) * screen correction factor 112/96 (112dpi instead of 96dpi as windows expects!)
//	CVect2 m_log2mm;
	CVect2 m_mm2dev;		// 112dpi / 25.4mmpi
	CVect2 m_dev2mm;		// 25.4mmpi / 112dpi
//	CVect2 m_log2dev;		// log2mm * mm2dev
//	CVect2 m_dev2log;		// dev2mm * mm2log

	
	// Reference points
	CVector2 m_ptFocusReal;			// focus point in real coords


	// dynamic mouse movement variables
	CPoint m_ptButtonDown, m_ptButtonUp, m_ptPrevMouseMove;
	bool m_bFirstMouseMove;
	enum
	{
		DYNA_NONE = 0,
		DYNA_TRANSLATE,
		DYNA_DIRECTION,
		DYNA_ZOOMROTATE,
		DYNA_ZOOMWINDOW,
	} m_nDynMoveType;
	CVector2 m_ptFocusRealInit;

	double m_DynaTranslateAmp;		// for amplified translations (default = 1.0)

	HCURSOR m_hCursor;	// to use different cursors in different modes
	HCURSOR m_hCursorMove;
	HCURSOR m_hCursorDir;

	// View options
	int m_nShowAxes;
	bool m_bIsoScaling;
	bool m_b3DView;					// true to use 3D view transforming
	bool m_bShowRegionBox;


// Operations
public:

	void SetScale(double xScale, double yScale);
	void SetScale(double scale) { SetScale(scale, scale); }
	void SetScale(CVector2& vtScale) { SetScale(vtScale.x, vtScale.y); }
	void GetScale(CVector2& vtScale) { vtScale = m_Scale; }
	void ZoomBy(double factor);

	void SetFocusReal(const CVector2& ptFocus);
	void GetFocusReal(CVector2& ptFocus) { ptFocus = m_ptFocusReal; }
	CPoint GetFocusPointDev();
	void GetFocusPointDev(CPoint& ptFocusDev, CRect& rectWnd);


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMove2DView)
	public:
	virtual void OnInitialUpdate();
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CMove2DView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	void AdjustScalingFactors();
	void FitRectToWindow(CRect& rectFitDev);
	
	void OnButtonDown(UINT nFlags, CPoint point);
	void OnButtonUp(UINT nFlags, CPoint point);

	// Generated message map functions
protected:
	//{{AFX_MSG(CMove2DView)
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MOVE2DVIEW_H__80D04CC3_6589_4373_9D3F_0B9DB8424200__INCLUDED_)
