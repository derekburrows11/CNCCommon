
#if !defined(SPLINEMATH_H__INCLUDED)
#define SPLINEMATH_H__INCLUDED




// See 'Spline Functions Approx.txt' for algorithms

template<class T>
void FindSpline1seg(T* n);		// p0 c0 c1 p1  -> 2 points, 2 control points

template<class T>
void FindSpline2segs(T* n);		// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points
template<class T>
void FindSpline2segsC0(T* n);	// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points
template<class T>
void FindSpline2segsC0C3(T* n);	// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points


template<class T>
void FindSplineSegFlatEnds(T* vtBez, const T* vtPts, int numPts);



////////////////////////////////////////////
// See 'Spline Functions Approx.txt' for algorithms



enum { p0, c0, c1, p1, c2, c3, p2};

template<class T>
void FindSpline1seg(T* n)	// p0 c0 c1 p1  -> 2 points, 2 control points
{
	T p0p1on3 = (n[p1] - n[p0]) / 3;
	n[c0] = n[p0] + p0p1on3;
	n[c1] = n[p1] - p0p1on3;
}

template<class T>
void FindSpline2segs(T* n)	// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points
{														//  0  1  2  3  4  5  6
/*	G2(0) = (L1^2*(p3 - p2) + L2^2*(p2 - p1)) / (L1*L2 + L1^2)
	G1(1) = L1/L2 * G2(0)	
	2G2(1) = 3p3 - 3p2 - G2(0)
	2G1(0) = 3p2 - 3p1 - G1(1)
*/
	T p0p1 = n[p1] - n[p0];
	T p1p2 = n[p2] - n[p1];
	double L1sq = p0p1.MagSq();
	double L2sq = p1p2.MagSq();
	double L1 = sqrt(L1sq);
	double L2 = sqrt(L2sq);
	T g[4];		// gradients / 3

	if (L2 == 0)
	{
		g[1] = g[2] = g[3] = 0;
		g[0] = p0p1 * 0.5;
	}
	else if (L1 == 0)
	{
		g[0] = g[1] = g[2] = 0;
		g[3] = p1p2 * 0.5;
	}
	else
	{
		g[2] = (p1p2 * L1sq + p0p1 * L2sq) / (3 * (L1*L2 + L1sq));
		g[1] = g[2] * (L1 / L2);
		g[3] = (p1p2 - g[2]) * 0.5;
		g[0] = (p0p1 - g[1]) * 0.5;
	}
	n[c0] = n[p0] + g[0];
	n[c1] = n[p1] - g[1];
	n[c2] = n[p1] + g[2];
	n[c3] = n[p2] - g[3];
}

template<class T>
void FindSpline2segsC0(T* n)	// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points
{		// first control point is defined
/*
	G2(0) = (L1^2 * (3p3 - 3p2) + L2^2 * (6p2 - 6p1 - 2g1)) / (4*L1*L2 + 3*L1^2)
	G1(1) = L1/L2 * G2(0)	
	2G2(1) = 3p3 - 3p2 - G2(0)
*/
	T p0p1 = n[p1] - n[p0];
	T p1p2 = n[p2] - n[p1];
	double L1sq = p0p1.MagSq();
	double L2sq = p1p2.MagSq();
	double L1 = sqrt(L1sq);
	double L2 = sqrt(L2sq);
	T g[4];		// gradients / 3

	if (L2 == 0)
		g[1] = g[2] = g[3] = 0;
	else if (L1 == 0)
	{
		g[1] = g[2] = 0;
		g[3] = p1p2 * 0.5;
	}
	else
	{
		g[0] = n[c0] - n[p0];
		g[2] = (p1p2 * L1sq + (p0p1 - g[0]) * (2*L2sq)) / (4*L1*L2 + 3*L1sq);
		g[1] = g[2] * (L1 / L2);
		g[3] = (p1p2 - g[2]) * 0.5;
	}
	n[c1] = n[p1] - g[1];
	n[c2] = n[p1] + g[2];
	n[c3] = n[p2] - g[3];
}

template<class T>
void FindSpline2segsC0C3(T* n)	// p0 c0 c1 p1 c2 c3 p2  -> 3 points, 4 control points
{		// first and last control points are defined
/*
	G2(0) = (L1^2 * (3p3 - 3p2 - g3) + L2^2 * (3p2 - 3p1 - g1)) / (2*L1*L2 + 2*L1^2)
	G1(1) = L1/L2 * G2(0)	
*/
	T p0p1 = n[p1] - n[p0];
	T p1p2 = n[p2] - n[p1];
	double L1sq = p0p1.MagSq();
	double L2sq = p1p2.MagSq();
	double L1 = sqrt(L1sq);
	double L2 = sqrt(L2sq);
	T g[4];		// gradients / 3

	g[0] = n[c0] - n[p0];
	g[3] = n[p2] - n[c3];
	g[2] = ((p1p2 - g[3]) * L1sq + (p0p1 - g[0]) * L2sq) / (2 * (L1*L2 + L1sq));
	g[1] = g[2] * (L1 / L2);

	n[c1] = n[p1] - g[1];
	n[c2] = n[p1] + g[2];
}



template<class T>
void FindSplineSegFlatEnds(T* vtBez, const T* vtPts, int numPts)
{
	// solve bezier spline sequentially
	// requires 3*numPts + 1 vtBez vectors
	if (numPts < 3)
	{
		if (numPts <= 0)
			return;
		else if (numPts == 1)
		{
			vtBez[0] = vtPts[0];
			return;
		}
		else if (numPts == 2)
		{
			vtBez[0] = vtPts[0];
			vtBez[3] = vtPts[1];
			FindSpline1seg(vtBez);
			return;
		}
	}
	ASSERT(numPts >= 3);

	vtBez[0] = vtPts[0];
	vtBez[3] = vtPts[1];
	vtBez[6] = vtPts[2];
	FindSpline2segs(vtBez);
	T* vtBezSolve = vtBez;
	for (int i = 3; i < numPts; i++)
	{
		vtBezSolve += 3;
		vtBezSolve[6] = vtPts[i];
		FindSpline2segsC0(vtBezSolve);
	}
}



#endif // !defined(SPLINEMATH_H__INCLUDED)
