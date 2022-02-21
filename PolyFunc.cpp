
#include "stdafx.h"

#include "Matrix.h"
#include "PolyFunc.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif




/*
	poly is [a0  a1  a2  a3]'
	thirds  [f(0)  f(1/3)  f(2/3)  f(1)]'
	bezier  [n0  c0  c1  n1]'
	bezier midpoint = 1/8*[1 3 3 1] * bezier

	val = a0 + a1*s + a2*s^2 + a3*s^3
	val = n0 + 3*(c0-n0)*s + 3*(c1-2c0+n0)*s^2 + (n1-3c1+3c0-n0)*s^3		 - using bezier2poly
	val = n0(1-3s+3s^2-s^3) + 3c0(s-2s^2+s^3) + 3c1(s^2-s^3) + n1(s^3)

	val = n0*(1-s)^3 + c0*3*s*(1-s)^2 + c1*3*s^2*(1-s) + n1*s^3

poly2thirds =
     1     0     0      0
     1   1/3   1/9   1/27
     1   2/3   4/9   8/27
     1     1     1      1

thirds2poly =
    1      0      0      0
   -5.5    9     -4.5    1
    9    -22.5   18     -4.5
   -4.5   13.5  -13.5    4.5

//////////

poly2bezier =
    1    0    0    0
    1  1/3    0    0
    1  2/3  1/3    0
    1    1    1    1

bezier2poly =
    1    0    0    0
   -3    3    0    0
    3   -6    3    0
   -1    3   -3    1

thirds2bezier =
    1      0      0      0
   -5/6    3     -1.5    1/3
    1/3   -1.5    3     -5/6
    0      0      0      1

bezier2thirds =
    1      0      0      0
 8/27    4/9    2/9   1/27
 1/27    2/9    4/9   8/27
    0      0      0      1

*/


///////////////////////////////////



/*
double CubicAt(double* poly, double x)		// Using Horner's polynomial rule
{
	return ((poly[3]*x + poly[2])*x + poly[1])*x + poly[0];
}
double CubicD1At(double* poly, double x)		// Using Horner's polynomial rule
{
	return (poly[3]*(1.5*x) + poly[2])*(2*x) + poly[1];
}
double CubicD2At(double* poly, double x)
{
	return poly[3]*(6*x) + poly[2]*2;
}
double CubicD3At(double* poly)
{
	return poly[3]*6;
}

CVector2 CubicAt(CVector2* poly, double x)		// Using Horner's polynomial rule
{
	return ((poly[3]*x + poly[2])*x + poly[1])*x + poly[0];
}
CVector2 CubicD1At(CVector2* poly, double x)		// Using Horner's polynomial rule
{
	return (poly[3]*(1.5*x) + poly[2])*(2*x) + poly[1];
}
CVector2 CubicD2At(CVector2* poly, double x)
{
	return poly[3]*(6*x) + poly[2]*2;
}
CVector2 CubicD3At(CVector2* poly)
{
	return poly[3]*6;
}

CVector CubicAt(CVector* poly, double x)		// Using Horner's polynomial rule
{
	return ((poly[3]*x + poly[2])*x + poly[1])*x + poly[0];
}
CVector CubicD1At(CVector* poly, double x)		// Using Horner's polynomial rule
{
	return (poly[3]*(1.5*x) + poly[2])*(2*x) + poly[1];
}
CVector CubicD2At(CVector* poly, double x)
{
	return poly[3]*(6*x) + poly[2]*2;
}
CVector CubicD3At(CVector* poly)
{
	return poly[3]*6;
}
*/

/*
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
*/
//template CVector CubicAt<CVector>(const CVector* poly, double s);
//template CVector CubicAt(const CVector* poly, double s);
//template CVector CubicAt(const CVector* poly, double s);






///////////////////////////////////


void BSplinePolyFromFitPoints(CVector2* poly, const CVector2* pts)
{
	// poly is approx bspline between pts[0] & pts[1] for s = [0,1]
	poly[0] = (pts[1] + pts[0]*4 + pts[-1]) * (1.0/6);
	poly[1] = (pts[1] - pts[-1]) * 0.5;
	poly[2] = (pts[1] - pts[0]*2 + pts[-1]) * 0.5;
	poly[3] = (pts[2] - pts[1]*3 + pts[0]*3 - pts[-1]) * (1.0/6);
}

void BSplineBezFromFitPoints(CVector2* bez, const CVector2* pts)
{
	// bezier is approx bspline between pts[0] & pts[1]
	bez[0] = (pts[1] + pts[0]*4 + pts[-1]) * (1.0/6);
	bez[1] = (pts[1] + pts[0]*2) * (1.0/3);
	bez[2] = (pts[1]*2 + pts[0]) * (1.0/3);
	bez[3] = (pts[2] + pts[1]*4 + pts[0]) * (1.0/6);
}


///////////////////////////////////
// Arc to bezier functions
///////////////////////////////////
// work best for arcs of small angle, <90 deg



double GetBezierArcControlNodeRadiusRatio(double angArc)
{
	// control node distance from end nodes = radius * 4/3 * (1 - cos(angArc/2)) / sin(angArc/2)   -> radius * 0.5523 for 90deg
	// this gives a bezier that touches arc midway
	return (4.0/3.0) * (1 - cos(angArc/2)) / sin(angArc/2);
}

double GetBezier90ArcControlNodeRadiusRatio()		// 90 deg arc
{
	// control node distance from end nodes = radius * 4/3 * (1 - cos(angArc/2)) / sin(angArc/2)   -> radius * 0.5523 for 90deg
	// this gives a bezier that touches arc midway
	static double val = (4.0/3.0) * (sqrt(2) - 1);
	return val;
}

double GetBezier45ArcControlNodeRadiusRatio()		// 45 deg arc
{
	// control node distance from end nodes = radius * 4/3 * (1 - cos(angArc/2)) / sin(angArc/2)   -> radius * 0.5523 for 90deg
	// this gives a bezier that touches arc midway
	static double val = (4.0/3.0) * (1 - cos(45*deg2rad/2)) / sin(45*deg2rad/2);
	return val;
}

void FitBezierToArcCentre(int dir, const CVector& ptInit, const CVector& ptCentre, const CVector& ptFinal, CVector& ptCi, CVector& ptCf)
{
	// uses initial and final points, arc centre point and arc direction CW or CCW
	CVector vtInitDirUnit = ptInit - ptCentre;		// radius vector
	vtInitDirUnit.Unit();
	ASSERT(vtInitDirUnit.z == 0);
	CVector2 vtDir = (CVector2)vtInitDirUnit;
	if (dir == CW)
		vtDir.Rotate90CW();		// rotate 90deg CW
	else if (dir == CCW)
		vtDir.Rotate90CCW();	// rotate 90deg CCW
	vtInitDirUnit = vtDir;

	// The following is copy of function below
	// uses initial and final points and initial direction
	CVector vtChord = ptFinal - ptInit;
	double magChord = vtChord.Mag();
	CVector vtChordUnit = vtChord / magChord;

	// Ang is angArc/2
	double cosAng, sinAngSq;
	if (vtInitDirUnit.z == 0 && vtChord.z == 0)
	{
		double sinAng;
		cosAng = vtChordUnit.x * vtInitDirUnit.x + vtChordUnit.y * vtInitDirUnit.y;	// ok if arc is in x-y plane
		sinAng = vtChordUnit.x * vtInitDirUnit.y - vtChordUnit.y * vtInitDirUnit.x;	// ok if arc is in x-y plane
		sinAngSq = sinAng*sinAng;
	}
	else
	{
		cosAng = 0;
		sinAngSq = 0;
		ASSERT(0);	// not handled yet!
	}

	// radius = magChord/2 / sin(angArc/2)
	// cnodeDist = 4.0/3 * radius * (1-cos(angArc/2)) / sin(angArc/2);	// this gives bezier that touches arc midway
	double cnodeDist;
	if (fabs(sinAngSq) > 1e-6)
		cnodeDist = 2 * magChord * (1 - cosAng) / (3 * sinAngSq);
	else
		cnodeDist = magChord / 3;

	ptCi = vtInitDirUnit * cnodeDist + ptInit;
	// refect vtCi about perp bisector of chord
	ptCf = ptCi + vtChordUnit * (magChord - 2 * cosAng * cnodeDist);
}

void FitBezierToArcDir(const CVector& ptInit, const CVector& vtInitDir, const CVector& ptFinal, CVector& ptCi, CVector& ptCf)
{
	// uses initial and final points and initial direction
	CVector vtInitDirUnit = vtInitDir / vtInitDir.Mag();
	CVector vtChord = ptFinal - ptInit;
	double magChord = vtChord.Mag();
	CVector vtChordUnit = vtChord / magChord;

	// Ang is angArc/2
	double cosAng, sinAngSq;
	if (vtInitDir.z == 0 && vtChord.z == 0)
	{
		double sinAng;
		cosAng = vtChordUnit.x * vtInitDirUnit.x + vtChordUnit.y * vtInitDirUnit.y;	// ok if arc is in x-y plane
		sinAng = vtChordUnit.x * vtInitDirUnit.y - vtChordUnit.y * vtInitDirUnit.x;	// ok if arc is in x-y plane
		sinAngSq = sinAng*sinAng;
	}
	else		// not in x-y plane
	{
		cosAng = dot(vtChordUnit, vtInitDirUnit);
		sinAngSq = cross(vtChordUnit, vtInitDirUnit).MagSq();
		//ASSERT(0);		// check this result!!
	}

	// radius = magChord/2 / sin(angArc/2)
	// cnodeDist = 4.0/3 * radius * (1-cos(angArc/2)) / sin(angArc/2);	// this gives bezier that touches arc midway
	double cnodeDist;
	if (fabs(sinAngSq) > 1e-6)
		cnodeDist = 2 * magChord * (1 - cosAng) / (3 * sinAngSq);
	else
		cnodeDist = magChord / 3;

	ptCi = vtInitDirUnit * cnodeDist + ptInit;
	// refect vtCi about perp bisector of chord
	ptCf = ptCi + vtChordUnit * (magChord - 2 * cosAng * cnodeDist);
}

void FitBezierToNoCurveArcDir(const CVector& ptInit, const CVector& vtInitDir, const CVector& ptFinal, CVector& ptCi, CVector& ptCf)
{
	// uses initial and final points and initial direction
	// only good for small angle changes
	// instead of an arc, can do cubic with no curve at ptInit if ptCf is in line with ptInit & ptCi
	CVector vtChord = ptFinal - ptInit;
	double cnodeDist = dot(vtChord, vtInitDir) / (3 * vtInitDir.MagSq());	// (= cnodeDist/magInitDir) is Chord component in initial direction - gives exact cubic (no quadratic)
	ptCi = vtInitDir * cnodeDist + ptInit;
	ptCf = vtInitDir * (2*cnodeDist) + ptInit;
}

void FitBezierToArc3Point(const CVector& ptInit, const CVector& ptMid, const CVector& ptFinal, CVector& ptCi, CVector& ptCf)
{
	// uses initial, mid and final points to find arc from 3 points on curve
	CVector vtChord = ptFinal - ptInit;
	CVector vtSeg1 = ptMid - ptInit;
	CVector vtSeg2 = ptFinal - ptMid;
	double magChord = vtChord.Mag();
	CVector vtChordUnit = vtChord / magChord;
	double magSeg1Seg2 = sqrt(vtSeg1.MagSq() * vtSeg2.MagSq());

	// Ang is angArc/2 is angle b/w vtSeg1 & vtSeg2, and is angle b/w vtChord & end tangents
	double cosAng, sinAngSq;
	CVector vtInitDirUnit;
	if (vtSeg1.z == 0 && vtSeg2.z == 0)
	{
		double sinAng;
		cosAng = (vtSeg2.x * vtSeg1.x + vtSeg2.y * vtSeg1.y) / magSeg1Seg2;	// ok if arc is in x-y plane
		sinAng = (vtSeg2.x * vtSeg1.y - vtSeg2.y * vtSeg1.x) / magSeg1Seg2;	// ok if arc is in x-y plane
		vtInitDirUnit.x = vtChordUnit.x * cosAng - vtChordUnit.y * sinAng;
		vtInitDirUnit.y = vtChordUnit.x * sinAng + vtChordUnit.y * cosAng;
		vtInitDirUnit.z = 0;
		sinAngSq = sinAng*sinAng;
	}
	else		// not in x-y plane
	{
		cosAng = dot(vtSeg2, vtSeg1) / magSeg1Seg2;
		CVector vtSinAng = cross(vtSeg2, vtSeg1) / magSeg1Seg2;
		vtInitDirUnit = cross(vtSinAng, vtChordUnit) + vtChordUnit * cosAng;
		sinAngSq = vtSinAng.MagSq();
		//ASSERT(0);		// check this result!!
	}

	//	radius = magChord/2 / sin(angArc/2)
	// cnodeDist = 4.0/3 * radius * (1-cos(angArc/2)) / sin(angArc/2);	// this gives bezier that touches arc midway
	double cnodeDist;
	if (fabs(sinAngSq) > 1e-6)
		cnodeDist = 2 * magChord * (1 - cosAng) / (3 * sinAngSq);
	else
		cnodeDist = magChord / 3;
	ptCi = vtInitDirUnit * cnodeDist + ptInit;
	// refect vtCi about perp bisector of chord
	ptCf = ptCi + vtChordUnit * (magChord - 2 * cosAng * cnodeDist);
}

///////////////////////////////////
// Point to bezier functions
///////////////////////////////////

void FitBezierTo4Points(const CVector& pt0, const CVector& pt1, const CVector& pt2, const CVector& pt3, CVector& ptCi, CVector& ptCf)
{
	// find s at pt1 & pt2 using segment lengths
	double s1 = (pt1 - pt0).Mag();
	double s2 = s1 + (pt2 - pt1).Mag();
	double sTotInv = 1.0 / (s2 + (pt3 - pt2).Mag());
	s1 *= sTotInv;		// normalise
	s2 *= sTotInv;
	// find bezier c0, c1, given n0, n1 and pt1/2 at s values  - solve 2x2
	//	val = n0(1-3s+3s^2-s^3) + 3c0(s-2s^2+s^3) + 3c1(s^2-s^3) + n1(s^3)
	double s1p2 = s1*s1;
	double s1p3 = s1*s1p2;
	double s2p2 = s2*s2;
	double s2p3 = s2*s2p2;
	CMatrix mxEqu(2,2);
	mxEqu.elem(0,0) = 3*(s1-2*s1p2+s1p3);
	mxEqu.elem(0,1) = 3*(s1p2-s1p3);
	mxEqu.elem(1,0) = 3*(s2-2*s2p2+s2p3);
	mxEqu.elem(1,1) = 3*(s2p2-s2p3);
	mxEqu.Invert();
	double mxVal00 = 1-3*s1+3*s1p2-s1p3;
	double mxVal01 = s1p3;
	double mxVal10 = 1-3*s2+3*s2p2-s2p3;
	double mxVal11 = s2p3;

	CMatrix vtVal(2), vtSolu(2);
	// solve for x
	vtVal[0] = pt1.x - mxVal00*pt0.x - mxVal01*pt3.x;
	vtVal[1] = pt2.x - mxVal10*pt0.x - mxVal11*pt3.x;
	vtSolu = mxEqu * vtVal;
	ptCi.x = vtSolu[0];
	ptCf.x = vtSolu[1];
	// solve for y
	vtVal[0] = pt1.y - mxVal00*pt0.y - mxVal01*pt3.y;
	vtVal[1] = pt2.y - mxVal10*pt0.y - mxVal11*pt3.y;
	vtSolu = mxEqu * vtVal;
	ptCi.y = vtSolu[0];
	ptCf.y = vtSolu[1];
	// solve for z
	vtVal[0] = pt1.z - mxVal00*pt0.z - mxVal01*pt3.z;
	vtVal[1] = pt2.z - mxVal10*pt0.z - mxVal11*pt3.z;
	vtSolu = mxEqu * vtVal;
	ptCi.z = vtSolu[0];
	ptCf.z = vtSolu[1];
}



///////////////////////////////////
//
///////////////////////////////////


double SolveCubicFor()
//double SolveCubicFor(double* poly, double pos, double sInit, double sMin, double sMax, int signVel, int signAcc)
{
// will find correct s point given sign of vel and acc
	ASSERT(0);
	return 0;
}






//////////////////

void OffsetBezier4Points(CVector2* nDest, const CVector2* nSrc, double dist)
{
	// Offsets bezier segment in nSrc by dist perpendicular to line
// method 1: offset end nodes and 1/3 points
	CVector2 poly[4];
	Bezier2Cubic(poly, nSrc);

	CVector2 ptsThirds[4];		// for offset as a thirds segment
	Bezier2Thirds(ptsThirds, nSrc);
	CVector2 vtNorm;

	vtNorm = poly[1];		// = (nSrc[1] - nSrc[0]) * 3;
	vtNorm.Rotate90();
	ptsThirds[0] += vtNorm * (dist / vtNorm.Mag());

//	ptsThirds[1] = CubicAt(poly, 1.0/3);
	vtNorm = CubicD1At(poly, 1.0/3);
	vtNorm.Rotate90();
	ptsThirds[1] += vtNorm * (dist / vtNorm.Mag());

//	ptsThirds[2] = CubicAt(poly, 2.0/3);
	vtNorm = CubicD1At(poly, 2.0/3);
	vtNorm.Rotate90();
	ptsThirds[2] += vtNorm * (dist / vtNorm.Mag());

	vtNorm = nSrc[3] - nSrc[2];
	vtNorm.Rotate90();
	ptsThirds[3] += vtNorm * (dist / vtNorm.Mag());

	Thirds2Bezier(nDest, ptsThirds);
}

void OffsetBezierEndDeriv(CVector2* nDest, const CVector2* nSrc, double dist)
{
// method 2: offset end nodes and end derivatives to left by dist
/*
	Curve (=1/radius) of parametric curve =
	C = (V x A) / |V|^3		where V = dP/dt and A = dV/dt

	Voffset = V * (R + dist)/R		where R = radius of curve = 1/curve
*/
	CVector2 poly[4];
	Bezier2Cubic(poly, nSrc);

	CVector2 vtV, vtA, vtNorm;
	double VmagSq, Vmag, curve, rad, Vratio;

	// start nodes
	vtV = poly[1];		// = (nSrc[1] - nSrc[0]) * 3;
	vtA = poly[2] * 2;
	VmagSq = vtV.MagSq();
	Vmag = sqrt(VmagSq);
	vtNorm = vtV;
	vtNorm.Rotate90();
	nDest[0] = nSrc[0] + vtNorm * (dist / Vmag);

	curve = cross(vtV, vtA) / (VmagSq * Vmag);
	rad = 1 / curve;
	Vratio = 1.0 - dist / rad;
	nDest[1] = nDest[0] + vtV * (Vratio / 3);

	// end nodes
	vtV = (nSrc[3] - nSrc[2]) * 3;
	vtA += poly[3] * 6;
	VmagSq = vtV.MagSq();
	Vmag = sqrt(VmagSq);
	vtNorm = vtV;
	vtNorm.Rotate90();
	nDest[3] = nSrc[3] + vtNorm * (dist / Vmag);

	curve = cross(vtV, vtA) / (VmagSq * Vmag);
	rad = 1 / curve;
	Vratio = 1.0 - dist / rad;
	nDest[2] = nDest[3] - vtV * (Vratio / 3);
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Matrix solving functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*
	Reference: LU-Factorization, Crout's or Doolittle's method, Kreyszig p.1011
	Takes (n^3)/3 operations to get L & U
	Then   n^2    operations to solve vtX
	Takes n^2 + 2n doubles of memory
	Should work well without the need for pivoting on a symmetric, positive definite matrix!
	for matrix A, vector b and unknow vector x
		mxA * vtX = vtVal
		find L, U such that:  mxA = L*U
		where L is a lower triangular with diagonals of 1 and U is upper triangular
			mxA * vtX = LU * vtX = vtVal
		->	Ly = vtVal	where U * vtX = y
		L and U are be combined in LU to save memory!
		LU diagonal belongs to U.  L's diagonal is 1's
*/
int TMatrix<double>::LUSolve(double* vtX, const double* vtVal) const
{
	ASSERT(w == h);
	if (w != h)
		return 0;			// matrix not square!
	int n = w;				// order of matrix
	int j, k, s, smax;		// j-row index, k-column index
	double dSum;
	const int maxN = 8;
	double statLU[maxN*maxN];	// reserve space for up to 8 square matrix without having to 'new' memory
	double statUinv[maxN];
	double statY[maxN];
	double* LU = statLU;
	double* Uinv = statUinv;
	double* Y = statY;

	if (n > maxN)
	{
//		ASSERT(n <= maxN);		// if not 'new' some memory for arrays
		LU = new double[n*(n+2)];	// LU[n*n]		delete only LU
		Uinv = LU + n*n;				// Uinv[n]
		Y = LU + n*(n+1);				// Y[n]
	}

	for (j = 0; j < n; j++)			// row
	{
		bool bIsLower = true;
		for (k = 0; k < n; k++)		// column
		{
			if (k == j)		// on diagonal element
				bIsLower = false;
			smax = bIsLower ? k : j;
			dSum = elem(j, k);
			for (s = 0; s < smax; s++)
				dSum -= LU[j*n + s] * LU[s*n + k];	// L[j, s] * U[s, k]
			if (bIsLower)
				LU[j*n + k] = dSum * Uinv[k];		// to lower
			else
			{
				LU[j*n + k] = dSum;		// using upper
				if (k == j)
					if (fabs(dSum) > 1e-10)		// limit of a small number!
						Uinv[j] = 1.0 / dSum;	// do div's for diagonal upper's only once, used again in back substitution
					else
					{
						TRACE1("Very small pivot value in TMatrix::LUSolve() of: %g\n", dSum);
						if (n <= 30)
							afxDump << "Matrix needs reordering or is singular in TMatrix::LUSolve():\n" << *this;
						ASSERT(0);			// matrix needs reordering or is singular
						return 0;
					}
			}
		}
	}
//CMatrix mxLU(n,n);
//mxLU.CopyFromArrayByRows(LU);
//afxDump << "mxLU: " << mxLU;
	// now find	Ly = vtVal
	for (j = 0; j < n; j++)		// row
	{
		dSum = vtVal[j];
		for (k = 0; k < j; k++)
			dSum -= LU[j*n + k] * Y[k];	// using lower
		Y[j] = dSum;
	}
	// now find	Ux = y
	for (j = n-1; j >= 0; j--)		// row
	{
		dSum = Y[j];
		for (k = j+1; k < n; k++)
			dSum -= LU[j*n + k] * vtX[k];	// using upper
		vtX[j] = dSum * Uinv[j];			// inverse stored earlier
	}
	if (LU != statLU)		// then was allocated with new
		delete[] LU;
	return 1;
}




/*
	LUFullSymSolve uses code shortcuts for applicable to fully symmetrical matricies
	as with least squares curve fitting
	where element Ajk = elemA[j+k]
	for a square matrix n*n there are (2n-1) elements in elemA
	vtX and vtVal are of length (n)
	solve for vtX:
		mxA * vtX = vtVal
*/
bool LUFullSymSolve(int n, const double* elemA, double* vtX, const double* vtVal)
{
	int j, k, s, smax;
	double dSum;
	bool bIsLower;
	double LU[4*4], Uinv[4], y[4];	// reserve space for up to cubic polys
	if (n > 4)
	{
		ASSERT(n <= 4);		// Order for LU solve is > 4 - need to increase array sizes
		return false;
	}

	for (j = 0; j < n; j++)			// row
	{
		bIsLower = true;
		for (k = 0; k < n; k++)		// column
		{
			if (k == j)		// on diagonal element
				bIsLower = false;
			smax = bIsLower ? k : j;
			dSum = elemA[j+k];
			for (s = 0; s < smax; s++)
				dSum -= LU[j*n + s] * LU[s*n + k];	// L[j, s] * U[s, k]
			if (bIsLower)
				LU[j*n + k] = dSum * Uinv[k];		// to lower
			else
			{
				LU[j*n + k] = dSum;		// using upper
				if (k == j)
					if (fabs(dSum) > 1e-20)		// limit of a small number!
						Uinv[j] = 1.0 / dSum;	// do div's for diagonal upper's only once, used again in back substitution
					else
					{
						ASSERT(0);			// marix needs reordering or is singular
						return false;
					}
			}
		}
	}
	// now find	Ly = b
	for (j = 0; j < n; j++)		// row
	{
		dSum = vtVal[j];
		for (k = 0; k < j; k++)
			dSum -= LU[j*n + k] * y[k];	// using lower
		y[j] = dSum;
	}
	// now find	Ux = y
	for (j = n-1; j >= 0; j--)		// row
	{
		dSum = y[j];
		for (k = j+1; k < n; k++)
			dSum -= LU[j*n + k] * vtX[k];	// using upper
		vtX[j] = dSum * Uinv[j];			// inverse stored earlier
	}
	return true;
}






//////////////////////////////////////
// Matrix inverse related functions
//////////////////////////////////////



// Downloaded from internet
// Original file 'Derek\Maths\Matrix Maths.cpp'
/*
* Solve a set of linear equations. a is a square matrix of coefficients.
* b is the right hand side. b is replaced by solution. 
* Target is replaced by its LU decomposition.
*/
int TMatrix<double>::LUSolve2(double* vtX, const double* vtVal) const
{
	ASSERT(w == h);		// should be square!
	if (w != h)
		return 0;				// matrix not square!

	int* arIndx = new int[h];

	for (int i = 0; i < h; i++)
		vtX[i] = vtVal[i];
	CMatrix mxLUD(*this);			// make a copy of matrix
	mxLUD.LUDecompose(arIndx, true);
	mxLUD.LUBackSubstitute(vtX, arIndx);
	delete[] arIndx;
	return 1;
}

/*
* Invert a matrix.
*/
int TMatrix<double>::LUInvert()
{
	ASSERT(w == h);		// should be square!
	if (w != h)
		return 0;				// matrix not square!

	int size = w;			//The size of the system.

	//temporary storage
	int* arIndx = new int[h];			// An array holding the permutations used by LU decomposition
	double* arCol = new double[size];
	
	CMatrix mxLUD(*this);		// copy to mxLUD
	mxLUD.LUDecompose(arIndx);
	
	//Do backsubstitution with the b matrix being all zeros except for
	//a 1 in the row that matches the column we're in.
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
			arCol[i] = 0;
		arCol[j] = 1;
		mxLUD.LUBackSubstitute(arCol, arIndx);
		for (i = 0; i < size; i++)
			elem(i,j) = arCol[i];			//plug values into result
	}

	delete[] arIndx;
	delete[] arCol;
	return 1;
}

int TMatrix<double>::LUInverse(TMatrix<T>& matrixI) const
{
	matrixI = *this;
	return matrixI.LUInvert();
/*
	ASSERT(w == h);		// should be square!
	if (w != h)
		return 0;				// matrix not square!

 	int size = w;			//The size of the system
	matrixI.SetSize(size, size);

	//temporary storage
	int* arIndx = new int[h];			// An array holding the permutations used by LU decomposition
	double* arCol = new double[size];
	
	CMatrix mxLUD(*this);		// copy to mxLUD
	mxLUD.LUDecompose(arIndx);

	//Do backsubstitution with the b matrix being all zeros except for
	//a 1 in the row that matches the column we're in.
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
			arCol[i] = 0;
		arCol[j] = 1;
		mxLUD.LUBackSubstitute(arCol, arIndx);
		//plug values into result
		for (i = 0; i < size; i++)
			matrixI.elem(i,j) = arCol[i];
	}

	delete[] arIndx;
	delete[] arCol;
	return 1;
*/
}


// Downloaded from internet
// Original file 'Derek\Maths\Matrix Maths.cpp'
/*
* luDecomposition performs LU Decomposition on itself.
* must be given an array to mark the row permutations and a flag
* to mark whether the number of permutations was even or odd.
* Reference: Numerical Recipes in C.
*/
//template<class T>
void TMatrix<double>::LUDecompose(int indx[], bool parity)
{
	int w = this->w;				// number of columns (make a local copy)
	int h = this->h;				// number of rows
	ASSERT(w == h);		// should be square!
	if (w != h)
		return;				// matrix not square!
	double* const matrix = arVal;

	int i, j, k, imax;			//imax is position of largest element in the row. i,j,k, are counters
	double amax, dum;				// amax is value of largest element in the row. 
	double* scaling = new double[h];	// scaling factor for each row is stored here
	double tiny = 1.0e-20;		// a small number != zero

	parity = true;	// Is the number of pivots even?  parity use???

	//Loop through rows to get the scaling information
	//The largest element in the row is the inverse of the scaling factor.
	for (i = 0; i < h; i++)
	{
		amax = 0;
		for (j = 0; j < w; j++)
			if ((dum = fabs(matrix[i*w + j])) > amax)
				amax = dum;
		if (amax == 0)
			TRACE1("Singular Matrix (row %i is all 0)\n", i);
		scaling[i] = 1.0 / amax;		//Save the scaling
	}

	//Loop through columns using Crout's Method.
	imax = 0;
	for (j = 0; j < w; j++)
	{
		//lower left corner?? upper left
		for (i = 0; i < j; i++)
		{
			dum = matrix[i*w + j];
			for (k = 0; k < i; k++)
				dum -= matrix[i*w + k] * matrix[k*w + j];
			matrix[i*w + j] = dum;
		}
		//Initialize search for largest element
		amax = 0;
		//upper right corner
		for (i = j; i < h; i++)
		{
			dum = matrix[i*w + j];
			for (k = 0; k < j; k++)
				dum -= matrix[i*w + k] * matrix[k*w + j];
			matrix[i*w + j] = dum;
			if (scaling[i] * fabs(dum) > amax)
			{
				amax = scaling[i]* fabs(dum);
				imax = i;
			}
		}

		//Change rows if it is necessary
		if (j != imax)
		{
			for (k = 0; k < w; k++)
			{
				dum = matrix[imax*w + k];
				matrix[imax*w + k] = matrix[j*w + k];
				matrix[j*w + k] = dum;
			}
			//Change parity
			parity = !parity;
			scaling[imax] = scaling[j];
		}
		//Mark the column with the pivot row.
		indx[j] = imax;

		//replace zeroes on the diagonal with a small number.
		if (matrix[j*w + j] == 0.0)
		{
			matrix[j*w + j] = tiny;
			TRACE2("Diagonal (%i,%i) is zero, setting to tiny!\n", j, j);
		}
		//Divide by the pivot element
		if (j != h)			// w or h ???
		{
			dum = 1.0 / matrix[j*w + j];
			for (i = j+1; i < h; i++)
				matrix[i*w + j] *= dum;
		}
	}
	delete[] scaling;
}

// Downloaded from internet
// Original file 'Derek\Maths\Matrix Maths.cpp'
/*
* Do the backsubstitution on matrix a which is the LU decomposition
* of the original matrix. b is the right hand side vector which is NX1. b 
* is replaced by the solution. indx is the array that marks the row
* permutations.
*/
void TMatrix<double>::LUBackSubstitute(double* b, int indx[]) const
{
	int w = this->w;			// number of columns (make a local copy)
	int h = this->h;			// number of rows

	//counters
	int i, ip, j, ii = -1;
	double sum = 0;
	
	for (i = 0; i < h; i++)		// scan through rows
	{		
		ip = indx[i];
		sum = b[ip];

		b[ip] = b[i];
		if (ii != -1)
			for (j = ii; j < i; j++)
				sum -= elem(i,j) * b[j];
		else
			if ( sum != 0)
				ii = i;
		b[i] = sum;
	}
	for (i = h-1; i >= 0; i--)
	{
		sum = b[i];
		for (j = i+1; j < w; j++)
			sum -= elem(i,j) * b[j];
		b[i] = sum / elem(i,i);
	}
}


