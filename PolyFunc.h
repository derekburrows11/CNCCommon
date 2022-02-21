
#if !defined(POLYFUNC_H__INCLUDED)
#define POLYFUNC_H__INCLUDED



#include "Vector.h"
//#include "Matrix.h"



// segment representation transforms

template<class T>
void Bezier2Cubic(T* poly, const T* bez)
{
	poly[0] = bez[0];
	poly[1] = 3 * (bez[1] - bez[0]);
	poly[2] = 3 * (bez[2] - 2 * bez[1] + bez[0]);
	poly[3] = 3 * (bez[1] - bez[2]) + (bez[3] - bez[0]);
}
template<class T>
void Cubic2Bezier(T* bez, const T* poly)
{
	bez[0] = poly[0];
	bez[1] = poly[0] + poly[1] / 3;
	bez[2] = poly[0] + (poly[1] * 2 + poly[2]) / 3.0;
	bez[3] = poly[0] + poly[1] + poly[2] + poly[3];
}

template<class T>
void Bezier2Thirds(T* thirds, const T* bez)
{
	thirds[0] = bez[0];
	thirds[1] = (bez[0]*8 + bez[1]*12 + bez[2]*6  + bez[3]  ) / 27.0;
	thirds[2] = (bez[0]   + bez[1]*6  + bez[2]*12 + bez[3]*8) / 27.0;
	thirds[3] = bez[3];
}
template<class T>
void Thirds2Bezier(T* bez, const T* thirds)
{
	bez[0] = thirds[0];
	bez[1] = (thirds[1]*18 - thirds[0]*5 - thirds[2]*9  + thirds[3]*2) / 6.0;
	bez[2] = (thirds[0]*2  - thirds[1]*9 + thirds[2]*18 - thirds[3]*5) / 6.0;
	bez[3] = thirds[3];
}

template<class T>
void Line2Bezier(T* bez, const T* pts)		// bez and pts can be equal!
{
	bez[0] = pts[0];
	bez[3] = pts[1];
	bez[1] = 2.0/3*bez[0] + 1.0/3*bez[3];
	bez[2] = 1.0/3*bez[0] + 2.0/3*bez[3];
}
template<class T>
void Line2Cubic(T* poly, const T* pts)		// poly and pts can be equal!
{
	poly[0] = pts[0];
	poly[1] = pts[1] - pts[0];
	poly[2] = 0.0;
	poly[3] = 0.0;
}


/*
void Bezier2Cubic(CVector2* poly, const CVector2* bez);
void Bezier2Cubic(CVector* poly, const CVector* bez);
void Cubic2Bezier(CVector2* bez, const CVector2* poly);
void Bezier2Thirds(CVector2* thirds, const CVector2* bez);
void Thirds2Bezier(CVector2* bez, const CVector2* thirds);
*/

// poly solving
/*
template<class T> T CubicAt(T* poly, double s);
template<class T> T CubicD1At(T* poly, double s);
template<class T> T CubicD2At(T* poly, double s);
template<class T> T CubicD3At(T* poly);
*/

template<class T>
T CubicAt(const T* poly, double s)		// Using Horner's polynomial rule
{
	return poly[0] + s*(poly[1] + s*(poly[2] + s*poly[3]));
}
template<class T>
T CubicD1At(const T* poly, double s)		// Using Horner's polynomial rule
{
	return poly[1] + 2*s*(poly[2] + 1.5*s*poly[3]);
}
template<class T>
T CubicD2At(const T* poly, double s)
{
	return 2*poly[2] + 6*s*poly[3];
}
template<class T>
T CubicD3At(const T* poly)
{
	return 6*poly[3];
}

///////////////////////////////////
// bezier solving
template<class T> T BezierAt(const T* bez, double s);
template<class T> T BezierD1At(const T* bez, double s);

/*
	val = n0*(1-s)^3 + c0*3*s*(1-s)^2 + c1*3*s^2*(1-s) + n1*s^3
*/

template<class T>
T BezierAt(const T* bez, double s)
{
	double oms = 1-s;
	double s2 = s*s;
	double oms2 = oms*oms;
	return oms2*oms*bez[0] + 3*s*oms2*bez[1] + 3*s2*oms*bez[2] + s2*s*bez[3];
}
template<class T>
T BezierD1At(const T* bez, double s)
{
	double oms = 1-s;
	double s2 = s*s;
	double oms2 = oms*oms;
	return -3*oms2*bez[0] + 3*oms*(oms-2*s)*bez[1] + 3*s*(2*oms-s)*bez[2] + 3*s2*bez[3];
}


///////////////////////////////////
// bezier length calculation
template<class T> double BezierLengthByThirds(const T* bez);

template<class T>
double BezierLengthByControl(const T* bez)
{
	double Li = (bez[1] - bez[0]).Mag();
	double Lm = (bez[2] - bez[1]).Mag();
	double Lf = (bez[3] - bez[2]).Mag();
	return Li + Lm + Lf;
}

template<class T>
double BezierLengthByThirds(const T* bez)
{
	// make vectors relative to bez[0] to reduce error
	T thirds1r0 = (bez[1]*12 + bez[2]*6  + bez[3]   - bez[0]*19) / 27.0;
	T thirds2r0 = (bez[1]*6  + bez[2]*12 + bez[3]*8 - bez[0]*26) / 27.0;
	double Li = thirds1r0.Mag();
	double Lm = (thirds2r0 - thirds1r0).Mag();
	double Lf = ((bez[3] - bez[0]) - thirds2r0).Mag();
//	return Li + Lm + Lf;
	double lenThirds = Li + Lm + Lf;
	double lenControl = BezierLengthByControl(bez);
	ASSERT(lenThirds <= lenControl + 1e-8);		// allow for calc error!
	ASSERT(lenThirds >= lenControl * 0.5);
	return lenThirds;
}



///////////////////////////////////
// poly span adjustment
template<class T> void CubicAdjustEnd(T* polyAdj, const T* poly, double s1);
template<class T> void CubicAdjustStart(T* polyAdj, const T* poly, double s0);
template<class T> void CubicAdjustStartEnd(T* polyAdj, const T* poly, double s0, double s1);
///////////////////////////////////
// bezier span adjustment
template<class T> void BezierAdjustEnd(T* bezAdj, const T* bez, double s1);
template<class T> void BezierAdjustStart(T* bezAdj, const T* bez, double s0);
template<class T> void BezierAdjustStartEnd(T* bezAdj, const T* bez, double s0, double s1);

template<class T> void BezierAdjustEndx(T* bezAdj, const T* bez, double s1);
template<class T> void BezierAdjustStartx(T* bezAdj, const T* bez, double s0);
template<class T> void BezierAdjustStartEndx(T* bezAdj, const T* bez, double s0, double s1);

template<class T>
void CubicAdjustEnd(T* polyAdj, const T* poly, double s1)
{	// Sadj = S(s1)				S=[0 1], Sadj=[0 s1]
	double s1P = s1*s1;
	polyAdj[0] = poly[0];
	polyAdj[1] = s1 * poly[1];
	polyAdj[2] = s1P * poly[2];
	polyAdj[3] = s1P*s1 * poly[3];
}

template<class T>
void CubicAdjustStart(T* polyAdj, const T* poly, double s0)
{	// Sadj = S(1-s0) + s0		S=[0 1], Sadj=[s0 1]
	double s0P = s0*s0;
	double Coef = 1 - s0;
	double CoefP = Coef*Coef;
	polyAdj[0] = poly[0] + s0*poly[1] + s0P*poly[2] + s0P*s0*poly[3];
	polyAdj[1] = Coef  * (poly[1] + 2*s0*poly[2] + 3*s0P*poly[3]);
	polyAdj[2] = CoefP * (poly[2] + 3*s0*poly[3]);
	polyAdj[3] = CoefP*Coef * poly[3];
}

template<class T>
void CubicAdjustStartEnd(T* polyAdj, const T* poly, double s0, double s1)
{	// Sadj = S(s1-s0) + s0		S=[0 1], Sadj=[s0 s1]
	double s0P = s0*s0;
	double Coef = s1 - s0;
	double CoefP = Coef*Coef;
	polyAdj[0] = poly[0] + s0*poly[1] + s0P*poly[2] + s0P*s0*poly[3];
	polyAdj[1] = Coef  * (poly[1] + 2*s0*poly[2] + 3*s0P*poly[3]);
	polyAdj[2] = CoefP * (poly[2] + 3*s0*poly[3]);
	polyAdj[3] = CoefP*Coef * poly[3];
}

///////////////////////////////////

template<class T>
void BezierAdjustEnd(T* bezAdj, const T* bez, double s1)
{	// Sadj = S(s1)				S=[0 1], Sadj=[0 s1]
	double s1P = s1*s1;
	double oms1 = 1-s1;
	double oms1P = oms1*oms1;
/*	Bezier2Cubic(poly, bez);
	CubicAdjustEnd(polyAdj, poly, s1);
	polyAdj[0] = bez[0];
	polyAdj[1] = 3*s1  * (bez[1] - bez[0]);
	polyAdj[2] = 3*s1P * (bez[2] - 2 * bez[1] + bez[0]);
	polyAdj[3] = s1P*s1 * (3*(bez[1] - bez[2]) + (bez[3] - bez[0]));
	Cubic2Bezier(bezAdj, poly);
	bezAdj[0] = bez[0];
	bezAdj[1] = bez[0] + s1*(bez[1] - bez[0]);;
	bezAdj[2] = bez[0] + 2*s1*(bez[1] - bez[0]) + s1P*(bez[2] - 2 * bez[1] + bez[0]);
	bezAdj[3] = bez[0] + 3*s1*(bez[1] - bez[0]) + 3*s1P*(bez[2] - 2 * bez[1] + bez[0]) + 3*s1P*s1 * (bez[1] - bez[2]) + (bez[3] - bez[0]);
*/
	bezAdj[0] = bez[0];
	bezAdj[1] = oms1*bez[0] + s1*bez[1];;
	bezAdj[2] = oms1P*bez[0] + 2*s1*oms1*bez[1] + s1P*bez[2];
	bezAdj[3] = oms1P*oms1*bez[0] + 3*s1*oms1P*bez[1] + 3*s1P*oms1*bez[2] + s1P*s1*bez[3];
}

template<class T>
void BezierAdjustStart(T* bezAdj, const T* bez, double s0)
{	// Sadj = S(1-s0) + s0		S=[0 1], Sadj=[s0 1]
	double s0P = s0*s0;
	double oms0 = 1 - s0;
	double oms0P = oms0*oms0;
/*
	bezAdj[0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	bezAdj[1] = oms0P*(oms0-oms0)*bez[0] + oms0*(3*s0*oms0+oms0*(oms0-2*s0))*bez[1] + s0*(3*s0*oms0+oms0*(2*oms0-s0))*bez[2] + s0P*(s0+oms0)*bez[3];
	bezAdj[2] = oms0*(oms0P - 2*oms0*oms0 + oms0P)*bez[0] + (3*s0*oms0P + 2*oms0*oms0*(oms0-2*s0) + oms0P*(s0-2*oms0))*bez[1] + (3*s0P*oms0 + 2*oms0*s0*(2*oms0-s0) + oms0P*(oms0-2*s0))*bez[2] + (s0P*s0 + 2*oms0*s0P + oms0P*s0)*bez[3];
	bezAdj[3] = (oms0P*oms0 - 3*oms0*oms0P + 3*oms0P*oms0 - oms0P*oms0)*bez[0] + 3*(s0*oms0P + oms0*oms0*(oms0-2*s0) + oms0P*(s0-2*oms0) + oms0P*oms0)*bez[1] + 3*(s0P*oms0 + oms0*s0*(2*oms0-s0) + oms0P*(oms0-2*s0) - oms0P*oms0)*bez[2] + (s0P*s0 + 3*oms0*s0P + 3*oms0P*s0 + oms0P*oms0)*bez[3];
*/
	bezAdj[0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	bezAdj[1] = oms0P*bez[1] + 2*s0*oms0*bez[2] + s0P*bez[3];
	bezAdj[2] = oms0*bez[2] + s0*bez[3];
	bezAdj[3] = bez[3];
}

template<class T>
void BezierAdjustStartEnd(T* bezAdj, const T* bez, double s0, double s1)
{	// Sadj = S(s1-s0) + s0		S=[0 1], Sadj=[s0 s1]
	double s0P = s0*s0;
	double oms0 = 1-s0;
	double oms0P = oms0*oms0;
	double s1P = s1*s1;
	double oms1 = 1-s1;
	double oms1P = oms1*oms1;

	double s0oms1 = s0 * oms1;
	double s1oms0 = s1 * oms0;
	double w2s0oms1 = 2*s0oms1 + s1oms0;
	double w2s1oms0 = 2*s1oms0 + s0oms1;

/*
//	Bezier2Cubic(poly, bez);
	poly[0] = bez[0];
	poly[1] = 3*(bez[1] - bez[0]);
	poly[2] = 3*(bez[2] - 2*bez[1] + bez[0]);
	poly[3] = (3*(bez[1] - bez[2]) + (bez[3] - bez[0]));
//	CubicAdjustStartEnd(polyAdj, poly, s0, s1);
	polyAdj[0] = bez[0] + s0*3*(bez[1] - bez[0]) + s0P*3*(bez[2] - 2*bez[1] + bez[0]) + s0P*s0*(3*(bez[1] - bez[2]) + (bez[3] - bez[0]));
	polyAdj[0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	polyAdj[1] = Coef * (3*(bez[1] - bez[0]) + 2*s0*3*(bez[2] - 2*bez[1] + bez[0]) + 3*s0P*(3*(bez[1] - bez[2]) + (bez[3] - bez[0])));
	polyAdj[1] = Coef*3 * (-oms0P*bez[0] + oms0*(oms0-2*s0)*bez[1] + s0*(2*oms0-s0)*bez[2] + s0P*bez[3]);
	polyAdj[2] = CoefP * (3*(bez[2] - 2*bez[1] + bez[0]) + 3*s0*(3*(bez[1] - bez[2]) + (bez[3] - bez[0])));
	polyAdj[2] = CoefP*3 * (oms0*bez[0] + (s0-2*oms0)*bez[1] + (oms0-2*s0)*bez[2] + s0*bez[3]);
	polyAdj[3] = CoefP*Coef*(3*(bez[1] - bez[2]) + (bez[3] - bez[0]));
//	Cubic2Bezier(bezAdj, polyAdj);
	bezAdj[0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	bezAdj[1] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3] + Coef * (-oms0P*bez[0] + oms0*(oms0-2*s0)*bez[1] + s0*(2*oms0-s0)*bez[2] + s0P*bez[3]);
	bezAdj[2] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3] + 2*Coef * (-oms0P*bez[0] + oms0*(oms0-2*s0)*bez[1] + s0*(2*oms0-s0)*bez[2] + s0P*bez[3]) + CoefP * (oms0*bez[0] + (s0-2*oms0)*bez[1] + (oms0-2*s0)*bez[2] + s0*bez[3]);
	bezAdj[3] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3] + 3*Coef * (-oms0P*bez[0] + oms0*(oms0-2*s0)*bez[1] + s0*(2*oms0-s0)*bez[2] + s0P*bez[3]) + CoefP*3 * (oms0*bez[0] + (s0-2*oms0)*bez[1] + (oms0-2*s0)*bez[2] + s0*bez[3]) + CoefP*Coef*(3*(bez[1] - bez[2]) + (bez[3] - bez[0]));
*/

	bezAdj[0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	bezAdj[1] = oms1*oms0P*bez[0] + oms0*w2s0oms1*bez[1] + s0*w2s1oms0*bez[2] + s1*s0P*bez[3];
	bezAdj[2] = oms0*oms1P*bez[0] + oms1*w2s1oms0*bez[1] + s1*w2s0oms1*bez[2] + s0*s1P*bez[3];
	bezAdj[3] = oms1P*oms1*bez[0] + 3*s1*oms1P*bez[1] + 3*s1P*oms1*bez[2] + s1P*s1*bez[3];

/*
	double Coef = s1 - s0;
	double CoefP = Coef*Coef;

	T arbezAdj[3][4];
	arbezAdj[0][0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	arbezAdj[0][1] = oms0P*(oms0-Coef)*bez[0] + oms0*(3*s0*oms0+Coef*(oms0-2*s0))*bez[1] + s0*(3*s0*oms0+Coef*(2*oms0-s0))*bez[2] + s0P*(s0+Coef)*bez[3];
	arbezAdj[0][2] = oms0*(oms0P - 2*Coef*oms0 + CoefP)*bez[0] + (3*s0*oms0P + 2*Coef*oms0*(oms0-2*s0) + CoefP*(s0-2*oms0))*bez[1] + (3*s0P*oms0 + 2*Coef*s0*(2*oms0-s0) + CoefP*(oms0-2*s0))*bez[2] + (s0P*s0 + 2*Coef*s0P + CoefP*s0)*bez[3];
	arbezAdj[0][3] = (oms0P*oms0 - 3*Coef*oms0P + 3*CoefP*oms0 - CoefP*Coef)*bez[0] + 3*(s0*oms0P + Coef*oms0*(oms0-2*s0) + CoefP*(s0-2*oms0) + CoefP*Coef)*bez[1] + 3*(s0P*oms0 + Coef*s0*(2*oms0-s0) + CoefP*(oms0-2*s0) - CoefP*Coef)*bez[2] + (s0P*s0 + 3*Coef*s0P + 3*CoefP*s0 + CoefP*Coef)*bez[3];

	arbezAdj[1][0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	arbezAdj[1][1] = oms0P*(oms0-Coef)*bez[0] + oms0*(3*s0*oms0+Coef*(oms0-2*s0))*bez[1] + s0*(3*s0*oms0+Coef*(2*oms0-s0))*bez[2] + s0P*(s0+Coef)*bez[3];
	arbezAdj[1][2] = oms0*oms1P*bez[0] + (3*s0*oms0P + 2*Coef*oms0*(oms0-2*s0) + CoefP*(s0-2*oms0))*bez[1] + (3*s0P*oms0 + 2*Coef*s0*(2*oms0-s0) + CoefP*(oms0-2*s0))*bez[2] + (s0P*s0 + 2*Coef*s0P + CoefP*s0)*bez[3];
	arbezAdj[1][3] = oms1P*oms1*bez[0] + 3*s1*oms1P*bez[1] + 3*s1P*oms1*bez[2] + s1P*s1*bez[3];

	arbezAdj[2][0] = oms0P*oms0*bez[0] + 3*s0*oms0P*bez[1] + 3*s0P*oms0*bez[2] + s0P*s0*bez[3];
	arbezAdj[2][1] = oms1*oms0P*bez[0] + oms0*(3*oms0*s1 - 2*Coef)*bez[1] + s0*(3*oms0*s1 - Coef)*bez[2] + s1*s0P*bez[3];
	arbezAdj[2][2] = oms0*oms1P*bez[0] + oms1*(3*oms1*s0 + 2*Coef)*bez[1] + s1*(3*oms1*s0 + Coef)*bez[2] + s0*s1P*bez[3];
	arbezAdj[2][3] = oms1P*oms1*bez[0] + 3*s1*oms1P*bez[1] + 3*s1P*oms1*bez[2] + s1P*s1*bez[3];

	for (int i1 = 0; i1 < 3; i1++)
		for (int i = 0; i < 4; i++)
		{
			double diff = (arbezAdj[i1][i] - bezAdj[i]).SumAbs();
			ASSERT(diff <= 1e-10);
		}
*/
}


template<class T>
void BezierAdjustEndx(T* bezAdj, const T* bez, double s1)
{	// Sadj = S(s1)				S=[0 1], Sadj=[0 s1]
	T poly[4];
	T polyAdj[4];
	Bezier2Cubic(poly, bez);
	CubicAdjustEnd(polyAdj, poly, s1);
	Cubic2Bezier(bezAdj, polyAdj);
}

template<class T>
void BezierAdjustStartx(T* bezAdj, const T* bez, double s0)
{	// Sadj = S(1-s0) + s0		S=[0 1], Sadj=[s0 1]
	T poly[4];
	T polyAdj[4];
	Bezier2Cubic(poly, bez);
	CubicAdjustStart(polyAdj, poly, s0);
	Cubic2Bezier(bezAdj, polyAdj);
}

template<class T>
void BezierAdjustStartEndx(T* bezAdj, const T* bez, double s0, double s1)
{	// Sadj = S(s1-s0) + s0		S=[0 1], Sadj=[s0 s1]
	T poly[4];
	T polyAdj[4];
	Bezier2Cubic(poly, bez);
	CubicAdjustStartEnd(polyAdj, poly, s0, s1);
	Cubic2Bezier(bezAdj, polyAdj);
}

template<class T>
void LineAdjustStartEnd(T* ptsAdj, const T* pts, double s0, double s1)
{
	ptsAdj[0] = (1-s0)*pts[0] + s0*pts[1];
	ptsAdj[1] = (1-s1)*pts[0] + s1*pts[1];
}


/*
double CubicAt(double* poly, double x);
double CubicD1At(double* poly, double x);
double CubicD2At(double* poly, double x);
double CubicD3At(double* poly);

CVector2 CubicAt(CVector2* poly, double x);
CVector2 CubicD1At(CVector2* poly, double x);
CVector2 CubicD2At(CVector2* poly, double x);
CVector2 CubicD3At(CVector2* poly);

CVector CubicAt(CVector* poly, double x);
CVector CubicD1At(CVector* poly, double x);
CVector CubicD2At(CVector* poly, double x);
CVector CubicD3At(CVector* poly);
*/

// Approx B-Spline method
void BSplinePolyFromFitPoints(CVector2* poly, const CVector2* pts);
void BSplineBezFromFitPoints(CVector2* bez, const CVector2* pts);

// Arc to bezier functions
double GetBezierArcControlNodeRadiusRatio(double angArc);
double GetBezier90ArcControlNodeRadiusRatio();
double GetBezier45ArcControlNodeRadiusRatio();
enum { CW = 1, CCW };
void FitBezierToArcCentre(int dir, const CVector& ptInit, const CVector& ptCentre, const CVector& ptFinal, CVector& ptCi, CVector& ptCf);
void FitBezierToArcDir(const CVector& ptInit, const CVector& vtInitDir, const CVector& ptFinal, CVector& ptCi, CVector& ptCf);
void FitBezierToNoCurveArcDir(const CVector& ptInit, const CVector& vtInitDir, const CVector& ptFinal, CVector& ptCi, CVector& ptCf);
void FitBezierToArc3Point(const CVector& ptInit, const CVector& ptMid, const CVector& ptFinal, CVector& ptCi, CVector& ptCf);

// Point fitting to bezier functions
void FitBezierTo4Points(const CVector& pt0, const CVector& pt1, const CVector& pt2, const CVector& pt3, CVector& ptCi, CVector& ptCf);

// Bezier segment offset functions
void OffsetBezier4Points(CVector2* nDest, const CVector2* nSrc, double dist);
void OffsetBezierEndDeriv(CVector2* nDest, const CVector2* nSrc, double dist);


// matrix solving functions
bool LUFullSymSolve(int n, const double* elemA, double* vtX, const double* vtVal);
//	Solves for x where: mxA * vtX = vtVal   where mxA[j,k] = elemA[j+k], mxA is n*n and elemA is 2n-1 long



#endif // !defined(POLYFUNC_H__INCLUDED)
