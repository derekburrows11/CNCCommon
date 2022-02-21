/////////////////////////////////////////////
// Vector.h
/////////////////////////////////////////////

#if !defined(VECTOR_H)
#define VECTOR_H

#include <math.h>
#include <iostream.h>
//#include <stdlib.h>		// for __min, __max


template<class T> class TMatrix;

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000



#ifndef ASSERT
	#includ <assert.h>			// use if not MFC
	#defin ASSERT assert			// ASSERT is MFC
#endif

// math definitions and macros

const double pi = acos(-1.0);
const double pi2 = pi * 2;
const double deg2rad = pi / 180;
const double rad2deg = 180 / pi;

#define NEARINT(x) ((int)floor((x) + 0.5))

/*
inline int sign(double val)
{
	if (val > 0) return 1;
	else if (val == 0) return 0;
	else return -1;
}
*/
#define sign(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)




//////////////////////////////////////////////
// template TVector2<T>
//////////////////////////////////////////////

template<class T> class TVector;	// notification for TVector2(const TVector<T>& vt3)

template<class T> class TVector2
{
public:
	T x, y;
/*	union
	{
		T arVal[2];
		struct { T x, y; };
	}; */
// constructors
	TVector2() {}
	TVector2(T valX, T valY) { x = valX; y = valY; }
	TVector2(T val) { x = y = val; }
//	explicit TVector2(const TVector<T>& vt3) { x = vt3.x; y = vt3.y; }
	TVector2(const TVector<T>& vt3) { x = vt3.x; y = vt3.y; }
	TVector2(const TMatrix<T>& matrix);

#ifdef __AFXWIN_H__
	TVector2(const CPoint pt)	{ x = (T)pt.x; y = (T)pt.y; }
	TVector2(const CSize  pt)	{ x = (T)pt.cx; y = (T)pt.cy; }
	operator CPoint() const		{ return CPoint(NEARINT(x), NEARINT(y)); }
	CPoint NearInt() { return CPoint(NEARINT(x), NEARINT(y)); }
#endif


// member functions
	void Set(T valX, T valY) { x = valX; y = valY; }
	void Set(T val) { x = y = val; }
	void Get(T& valX, T& valY) const { valX = x; valY = y; }
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void SetFromArray(const T2* arVal) { x = (T)arVal[0]; y = (T)arVal[1]; }
	bool IsVector() const { return true; }
	int Length() const { return 2; }	
	T* GetArray() { return &x; }
	void Transpose() {}
	void Transpose(TVector2<T>& vtT) const { vtT.x = x; vtT.y = y; }
	
	T& operator[](int idx) { ASSERT(idx >= 0 && idx < 2); return (&x)[idx]; }
	T operator[](int idx) const { ASSERT(idx >= 0 && idx < 2); return (&x)[idx]; }

	void Min(const TVector2<T>& vt);
	void Max(const TVector2<T>& vt);
	void Min(const TVector2<T>& vt1, const TVector2<T>& vt2);
	void Max(const TVector2<T>& vt1, const TVector2<T>& vt2);
	T MinElem() const { return x<y ? x:y; }
	T MaxElem() const { return x>y ? x:y; }


	// logical functions
	bool Any() const { return x || y; }
	bool All() const { return x && y; }

	// vector maths
	double MagSq() const { return x*x + y*y; }
	double Mag() const { return sqrt(x*x + y*y); }
	double Sum() const { return x + y; }
	double SumAbs() const { return fabs(x) + fabs(y); }
	void Unit(TVector2<T>& vtU) const;
	void Unit();
	void Neg()			{ x = (T)-x; y = (T)-y; }
	void Rotate90()	{ T t = x; x = (T)-y; y = t; }	// 90 deg CCW
	void Rotate90CCW(){ T t = x; x = (T)-y; y = t; }	// 90 deg CCW
	void Rotate180()	{ x = (T)-x; y = (T)-y; }
	void Rotate270()	{ T t = x; x = y; y = (T)-t; }	// 90 deg CW
	void Rotate90CW()	{ T t = x; x = y; y = (T)-t; }	// 90 deg CW
	void Scale(const TVector2<T>& lhs, const TVector2<T>& rhs) { x = (T)(lhs.x*rhs.x); y = (T)(lhs.y*rhs.y); }
	void Scale(const TVector<T>& lhs, const TVector2<T>& rhs)  { x = (T)(lhs.x*rhs.x); y = (T)(lhs.y*rhs.y); }
	void Scale(const TVector2<T>& lhs, const TVector<T>& rhs)  { x = (T)(lhs.x*rhs.x); y = (T)(lhs.y*rhs.y); }
	void Prod(const TMatrix<T>& lhs, const TVector<T>& rhs);
	void ProdPart(const TMatrix<T>& lhs, const TVector<T>& rhs);	// only calc for TVector2 part of result
	double Dot(const TVector2<T>& vt2) const;
	double Cross(const TVector2<T>& rhs) const;




	// assignment operators
	TVector2& operator+=(const TVector2& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x += rhs.x;
		y += rhs.y;
		return *this;
	#pragma warning(default: 4244)
	}
	TVector2& operator-=(const TVector2& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	#pragma warning(default: 4244)
	}
//	TVector2<T>& operator+=(const T val) { x += val; y += val; return *this; }
//	TVector2<T>& operator-=(const T val) { x -= val; y -= val; return *this; }
	TVector2<T>& operator*=(const T val) { x = T(x * val); y = T(y * val); return *this; }
	TVector2<T>& operator/=(const T val) { x = T(x / val); y = T(y / val); return *this; }
	
	TVector2<T> operator-() const { return TVector2<T>((T)-x, (T)-y); }
	int operator==(const TVector2& rhs) const { return x == rhs.x && y == rhs.y; }

};


//////////////////////////////////////////////
// template TVector<T>
//////////////////////////////////////////////

template<class T> class TVector
{
public:
	T x, y, z;
/*	union
	{
		T arVal[3];
		struct { T x, y, z; };
	}; */
// constructors
	TVector() {}
	TVector(T valX, T valY, T valZ) { x = valX; y = valY; z = valZ; }
	TVector(T val) { x = y = z = val; }
	TVector(const TVector2<T>& vt2) { x = vt2.x; y = vt2.y; z = 0; }
		template<class T2>
	TVector(const TVector<T2>& vect) { x = (T)vect.x; y = (T)vect.y; z = (T)vect.z; }
	TVector(const TMatrix<T>& matrix);

#ifdef __AFXWIN_H__
	operator CPoint() const { return CPoint(NEARINT(x), NEARINT(y)); }
#endif


// member functions
	void Set(T valX, T valY, T valZ) { x = valX; y = valY; z = valZ; }
	void Set(T val) { x = y = z = val; }
	void Get(T& valX, T& valY, T& valZ) const { valX = x; valY = y; valZ = z; }
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void SetFromArray(const T2* arVal) { x = (T)arVal[0]; y = (T)arVal[1]; z = (T)arVal[2]; }
	bool IsVector() const { return true; }
	int Length() const { return 3; }	
	T* GetArray() { return &x; }
	void Transpose() {}
	void Transpose(TVector<T>& vtT) const { vtT.x = x; vtT.y = y; vtT.z = z; }

	T& operator[](int idx) { ASSERT(idx >= 0 && idx < 3); return (&x)[idx]; }
	T operator[](int idx) const { ASSERT(idx >= 0 && idx < 3); return (&x)[idx]; }

	void Abs(TVector<T>& vtAbs) const { vtAbs.x = (T)fabs(x); vtAbs.y = (T)fabs(y); vtAbs.z = (T)fabs(z); }
	void Min(const TVector<T>& vt);
	void Max(const TVector<T>& vt);
	void Min(const TVector<T>& vt1, const TVector<T>& vt2);
	void Max(const TVector<T>& vt1, const TVector<T>& vt2);
	T MinElem() const { return x<y ? (x<z? x:z):(y<z? y:z); }
	T MaxElem() const { return x>y ? (x>z? x:z):(y>z? y:z); }
	T MinAbsElem() const;
	T MaxAbsElem() const;
	int MinAbsElemAxis() const;
	int MaxAbsElemAxis() const;

	// logical functions
	bool Any() const { return x || y || z; }
	bool All() const { return x && y && z; }

	// vector maths
	double MagSq() const { return x*x + y*y + z*z; }
	double Mag() const { return sqrt(x*x + y*y + z*z); }
	double Sum() const { return x + y + z; }
	double SumAbs() const { return fabs(x) + fabs(y) + fabs(z); }
	void Unit(TVector<T>& vtU) const;
	void Unit();
	void Neg() { x = (T)-x; y = (T)-y; z = (T)-z; }
	void Neg(TVector<T>& vtN) const { vtN.x = (T)-x; vtN.y = (T)-y; vtN.z = (T)-z; }
	void Rotate90(TVector<T>& vtR)  const { vtR.x = (T)-y; vtR.y = x; vtR.z = z; }
	void Rotate180(TVector<T>& vtR) const { vtR.x = (T)-x; vtR.y = (T)-y; vtR.z = z; }
	void Rotate270(TVector<T>& vtR) const { vtR.x = y; vtR.y = (T)-x; vtR.z = z; }
	void Scale(const TVector<T>& lhs, const TVector<T>& rhs) { x = (T)(lhs.x*rhs.x); y = (T)(lhs.y*rhs.y); z = (T)(lhs.z*rhs.z); }
	void Prod(const TMatrix<T>& lhs, const TVector<T>& rhs);
	void ProdPart(const TMatrix<T>& lhs, const TVector<T>& rhs);	// only calc for TVector part of result
	void ProdPart(const TMatrix<T>& lhs, const TVector2<T>& rhs);	// only calc for TVector part of result
	double Dot(const TVector<T>& vt2) const;
	void Cross(const TVector<T>& lhs, const TVector<T>& rhs);


	// assignment operators
	TVector<T>& operator=(const T val) { x = y = z = val; return *this; }
//	TVector<T>& operator=(const TVector<T>& vt) { x = vt.x; y = vt.y; z = vt.z; return *this; }
		template<class T2>
	TVector<T>& operator=(const TVector<T2>& vt) { x = (T)vt.x; y = (T)vt.y; z = (T)vt.z; return *this; }
		template<class T2>
	TVector<T>& operator=(const TVector2<T2>& vt2) { x = vt2.x; y = vt2.y; z = 0; return *this; }
		template<class T2>
	TVector<T>& operator=(const TMatrix<T2>& matrix)
	{
		ASSERT(matrix.IsVector());
		ASSERT(matrix.Length() <= 3);
		int iMax = matrix.Length();
		for (int i = 0; i < iMax; i++)
			arVal[i] = (T)matrix[i];
		return *this;
	}
	TVector<T>& operator=(const TMatrix<T>& matrix);
	TVector<T>& operator+=(const TVector<T>& vt);
	TVector<T>& operator-=(const TVector<T>& vt);
//	TVector<T>& operator*=(const T val) { x = T(x * val); y = T(y * val); z = T(z * val); return *this; }
//	TVector<T>& operator/=(const T val) { x = T(x / val); y = T(y / val); z = T(z / val); return *this; }
	TVector<T>& operator*=(const T val);
	TVector<T>& operator/=(const T val);
	TVector<T> operator-() const { return TVector<T>((T)-x, (T)-y, (T)-z); }
	int operator==(const TVector& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z; }

	// conversion operators
//	operator TVector2<T>&() { return *(TVector2<T>*)&x; }



/*	{
		#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x += vt.x; y += vt.y; z += vt.z;
		#pragma warning(default: 4244)
		return *this;
	}
*/	


// functions same as TMatrix counterparts
//	operator T*() { return arVal; }	// gives overload ambiguity!

/*
// TPoint3<T> conversions
	TVector<T>& operator=(const TPoint2<T>& pt) { x = pt.x; y = pt.y; z = 0; return *this; }
	TVector<T>& operator=(const TPoint3<T>& pt) { x = pt.x; y = pt.y; z = pt.z; return *this; }
	operator TPoint2<T>() const { return TPoint2<T>(x, y); }
//	operator TPoint3<T>() const { return TPoint3<T>(x, y, z); }
//	operator TPoint3<T>*() { return (TPoint3<T>*)&x; }
	operator TPoint3<T>&() { return *(TPoint3<T>*)&x; }
*/
};





///////////////////////////////////////////
///////////////////////////////////////////

typedef TVector<float> CVectorf;
typedef TVector<double> CVector;

typedef TVector2<float> CVector2f;
typedef TVector2<double> CVector2;

typedef TVector2<float> CVect2f;
typedef TVector2<double> CVect2;

///////////////////////////////////////////
// TVector member functions
///////////////////////////////////////////




template<class T>
TVector2<T>::TVector2(const TMatrix<T>& matrix)
{
	int iLen = matrix.Length();
	ASSERT(matrix.IsVector());
	ASSERT(iLen <= 2);
	x = (iLen > 0) ? matrix[0] : 0;
	y = (iLen > 1) ? matrix[1] : 0;
}

template<class T>
TVector<T>::TVector(const TMatrix<T>& matrix)
{
	int iLen = matrix.Length();
	ASSERT(matrix.IsVector());
	ASSERT(iLen <= 3);
	x = (iLen > 0) ? matrix[0] : 0;
	y = (iLen > 1) ? matrix[1] : 0;
	z = (iLen > 2) ? matrix[2] : 0;
}

template<class T>
void TVector2<T>::Min(const TVector2<T>& vt)
{
	if (vt.x < x) x = vt.x;
	if (vt.y < y) y = vt.y;
}
template<class T>
void TVector2<T>::Max(const TVector2<T>& vt)
{
	if (vt.x > x) x = vt.x;
	if (vt.y > y) y = vt.y;
}
template<class T>
void TVector2<T>::Min(const TVector2<T>& vt1, const TVector2<T>& vt2)
{
	x = (vt1.x <= vt2.x) ? vt1.x : vt2.x;
	y = (vt1.y <= vt2.y) ? vt1.y : vt2.y;
}
template<class T>
void TVector2<T>::Max(const TVector2<T>& vt1, const TVector2<T>& vt2)
{
	x = (vt1.x >= vt2.x) ? vt1.x : vt2.x;
	y = (vt1.y >= vt2.y) ? vt1.y : vt2.y;
}

template<class T>
void TVector<T>::Min(const TVector<T>& vt)
{
	if (vt.x < x) x = vt.x;
	if (vt.y < y) y = vt.y;
	if (vt.z < z) z = vt.z;
}
template<class T>
void TVector<T>::Max(const TVector<T>& vt)
{
	if (vt.x > x) x = vt.x;
	if (vt.y > y) y = vt.y;
	if (vt.z > z) z = vt.z;
}
template<class T>
void TVector<T>::Min(const TVector<T>& vt1, const TVector<T>& vt2)
{
	x = (vt1.x <= vt2.x) ? vt1.x : vt2.x;
	y = (vt1.y <= vt2.y) ? vt1.y : vt2.y;
	z = (vt1.z <= vt2.z) ? vt1.z : vt2.z;
}
template<class T>
void TVector<T>::Max(const TVector<T>& vt1, const TVector<T>& vt2)
{
	x = (vt1.x >= vt2.x) ? vt1.x : vt2.x;
	y = (vt1.y >= vt2.y) ? vt1.y : vt2.y;
	z = (vt1.z >= vt2.z) ? vt1.z : vt2.z;
}

template<class T>
T TVector<T>::MinAbsElem() const
{
	T minElem = fabs(x) < fabs(y) ? x : y;
	return fabs(z) < fabs(minElem) ? z : minElem;
}
template<class T>
T TVector<T>::MaxAbsElem() const
{
	T maxElem = fabs(x) > fabs(y) ? x : y;
	return fabs(z) > fabs(maxElem) ? z : maxElem;
}

template<class T>
int TVector<T>::MinAbsElemAxis() const
{
	int minAx = fabs(y) < fabs(x) ? 1 : 0;
	return fabs(z) < fabs(arVal[minAx]) ? 2 : minAx;
}
template<class T>
int TVector<T>::MaxAbsElemAxis() const
{
	int maxAx = fabs(y) > fabs(x) ? 1 : 0;
	return fabs(z) > fabs(arVal[maxAx]) ? 2 : maxAx;
}

template<class T>
void TVector2<T>::Unit(TVector2<T>& vectU) const
{
	double invMag = 1 / sqrt(x*x + y*y);		// do division only once!
	vectU.x = x * invMag; vectU.y = y * invMag;
}
template<class T>
void TVector2<T>::Unit()
{
	double invMag = 1 / sqrt(x*x + y*y);		// do division only once!
	x *= invMag; y *= invMag;
}
template<class T>
void TVector<T>::Unit(TVector<T>& vectU) const
{
	double invMag = 1 / sqrt(x*x + y*y + z*z);		// do division only once!
	vectU.x = x * invMag; vectU.y = y * invMag; vectU.z = z * invMag;
}
template<class T>
void TVector<T>::Unit()
{
	double invMag = 1 / sqrt(x*x + y*y + z*z);		// do division only once!
	x *= invMag; y *= invMag; z *= invMag;
}


template<class T>
void TVector2<T>::Prod(const TMatrix<T>& lhs, const TVector<T>& rhs)
{
	const int inner = rhs.Length();
	ASSERT(lhs.w == inner);
	ASSERT(lhs.h == Length());
	int mxIdx = 0;
	T elem;
	for (int row = 0; row < Length(); row++)
	{
		elem = 0;
		for (int i = 0; i < inner; i++)
			elem += lhs[mxIdx++] * rhs[i];
		(&x)[row] = elem;
	}
}

template<class T>
void TVector<T>::Prod(const TMatrix<T>& lhs, const TVector<T>& rhs)
{
	const int inner = rhs.Length();
	ASSERT(lhs.w == inner);
	ASSERT(lhs.h == Length());
	int mxIdx = 0;
	T elem;
	for (int row = 0; row < Length(); row++)
	{
		elem = 0;
		for (int i = 0; i < inner; i++)
			elem += lhs[mxIdx++] * rhs[i];
		(&x)[row] = elem;
	}
}

template<class T>
void TVector2<T>::ProdPart(const TMatrix<T>& lhs, const TVector<T>& rhs)
{
	const int inner = rhs.Length();
	ASSERT(lhs.w == inner);
	ASSERT(lhs.h >= Length());
	int mxIdx = 0;
	T elem;
	for (int row = 0; row < Length(); row++)
	{
		elem = 0;
		for (int i = 0; i < inner; i++)
			elem += lhs[mxIdx++] * rhs[i];
		(&x)[row] = elem;
	}
}

template<class T>
void TVector<T>::ProdPart(const TMatrix<T>& lhs, const TVector<T>& rhs)
{
	const int inner = rhs.Length();
	ASSERT(lhs.w == inner);
	ASSERT(lhs.h >= Length());
	int mxIdx = 0;
	T elem;
	for (int row = 0; row < Length(); row++)
	{
		elem = 0;
		for (int i = 0; i < inner; i++)
			elem += lhs[mxIdx++] * rhs[i];
		(&x)[row] = elem;
	}
}

template<class T>
void TVector<T>::ProdPart(const TMatrix<T>& lhs, const TVector2<T>& rhs)
{
	const int inner = rhs.Length();
	ASSERT(lhs.w >= inner);
	ASSERT(lhs.h >= Length());
	ASSERT(inner == 2);
	for (int row = 0; row < Length(); row++)
	{
		int mxIdx = lhs.w * row;
		(&x)[row] = lhs[mxIdx++] * rhs[0] + lhs[mxIdx] * rhs[1];
	}
}


template<class T>
double TVector2<T>::Dot(const TVector2<T>& vt2) const
{
	return (x * vt2.x + y * vt2.y);
}
template<class T>
double TVector<T>::Dot(const TVector<T>& vt2) const
{
	return (x * vt2.x + y * vt2.y + z * vt2.z);
}

template<class T>
double TVector2<T>::Cross(const TVector2<T>& rhs) const
{
	return (x*rhs.y - y*rhs.x);
}
template<class T>
void TVector<T>::Cross(const TVector<T>& lhs, const TVector<T>& rhs)
{
	x = lhs.y*rhs.z - lhs.z*rhs.y;
	y = lhs.z*rhs.x - lhs.x*rhs.z;
	z = lhs.x*rhs.y - lhs.y*rhs.x;
}




///////////////////////////////////////////
// TVector member operators
///////////////////////////////////////////

template<class T> inline
TVector<T>& TVector<T>::operator+=(const TVector<T>& vt)
{
	x += vt.x; y += vt.y; z += vt.z;
	return *this;
}

template<class T> inline
TVector<T>& TVector<T>::operator-=(const TVector<T>& vt)
{
	x -= vt.x; y -= vt.y; z -= vt.z;
	return *this;
}

template<class T> inline
TVector<T>& TVector<T>::operator*=(const T val)
{
	x *= val; y *= val; z *= val;
	return *this;
}

template<class T> inline
TVector<T>& TVector<T>::operator/=(const T val)
{
	x /= val; y /= val; z /= val;
	return *this;
}


template<class T>
TVector<T>& TVector<T>::operator=(const TMatrix<T>& matrix)
{
	ASSERT(matrix.IsVector());
	ASSERT(matrix.Length() <= 3);
	int iMax = matrix.Length();
	for (int i = 0; i < iMax; i++)
		arVal[i] = matrix[i];
	return *this;
}


/////////////////////////////////////
// TVector associated functions
/////////////////////////////////////

template<class T>
double dot(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	return (lhs.x * rhs.x) + (lhs.y * rhs.y);
}
template<class T>
double cross(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	return (lhs.x * rhs.y) - (lhs.y * rhs.x);
}

template<class T> inline
double dot(const TVector<T>& lhs, const TVector<T>& rhs)
{
	return (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
}
template<class T>
TVector<T> cross(const TVector<T>& lhs, const TVector<T>& rhs)
{
	return TVector<T>(lhs.y*rhs.z - lhs.z*rhs.y, lhs.z*rhs.x - lhs.x*rhs.z, lhs.x*rhs.y - lhs.y*rhs.x);
}


///////////////////////////////////////////
// TVector associated operators
///////////////////////////////////////////


template<class T> inline
TVector2<T> operator+(const TVector2<T>& lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs.x + rhs.x, lhs.y + rhs.y); }
template<class T> inline
TVector2<T> operator+(const TVector2<T>& lhs, T rhs)
	{	return TVector2<T>(lhs.x + rhs, lhs.y + rhs); }
template<class T> inline
TVector2<T> operator+(T lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs + rhs.x, lhs + rhs.y); }

template<class T> inline
TVector<T> operator+(const TVector<T>& lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z); }
template<class T> inline
TVector<T> operator+(const TVector<T>& lhs, T rhs)
	{	return TVector<T>(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs); }
template<class T> inline
TVector<T> operator+(T lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z); }

template<class T> inline
TVector2<T> operator-(const TVector2<T>& lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs.x - rhs.x, lhs.y - rhs.y); }
template<class T> inline
TVector2<T> operator-(const TVector2<T>& lhs, T rhs)
	{	return TVector2<T>(lhs.x - rhs, lhs.y - rhs); }
template<class T> inline
TVector2<T> operator-(T lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs - rhs.x, lhs - rhs.y); }

template<class T> inline
TVector<T> operator-(const TVector<T>& lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z); }
template<class T> inline
TVector<T> operator-(const TVector<T>& lhs, T rhs)
	{	return TVector<T>(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs); }
template<class T> inline
TVector<T> operator-(T lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z); }


template<class T> inline
TVector2<T> operator*(const TVector2<T>& lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs.x * rhs.x, lhs.y * rhs.y); }
template<class T> inline
TVector2<T> operator*(const TVector2<T>& lhs, double rhs)
	{	return TVector2<T>(lhs.x * rhs, lhs.y * rhs); }
template<class T> inline
TVector2<T> operator*(double lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs * rhs.x, lhs * rhs.y); }

template<class T> inline
TVector2<T> operator*(const TVector2<T>& lhs, const TVector<T>& rhs)
	{	return TVector2<T>(lhs.x * rhs.x, lhs.y * rhs.y); }
template<class T> inline
TVector2<T> operator*(const TVector<T>& lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs.x * rhs.x, lhs.y * rhs.y); }

#ifdef __AFXWIN_H__
template<class T>		// use POINT not CPoint - CPoint(DWORD pt) enables this to be used with int's!
TVector2<T> operator*(const TVector2<T>& lhs, const POINT rhs)
	{	return TVector2<T>(lhs.x * rhs.x, lhs.y * rhs.y); }
template<class T>		// use SIZE not CSize - CSize(DWORD pt) enables this to be used with int's!
TVector2<T> operator*(const TVector2<T>& lhs, const SIZE rhs)
	{	return TVector2<T>(lhs.x * rhs.cx, lhs.y * rhs.cy); }
#endif



template<class T> inline
TVector<T> operator*(const TVector<T>& lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z); }
template<class T> inline
TVector<T> operator*(const TVector<T>& lhs, T rhs)
	{	return TVector<T>(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs); }
template<class T> inline
TVector<T> operator*(T lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z); }
template<class T> inline
TVector<T> operator*(const TVector<T>& lhs, int rhs)
	{	return TVector<T>(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs); }
template<class T> inline
TVector<T> operator*(int lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z); }




template<class T> inline
TVector2<T> operator/(const TVector2<T>& lhs, const TVector2<T>& rhs)
	{	return TVector2<T>(lhs.x / rhs.x, lhs.y / rhs.y); }
template<class T> inline
TVector2<T> operator/(const TVector2<T>& lhs, double rhs)
	{	return TVector2<T>(lhs.x / rhs, lhs.y / rhs); }

template<class T> inline
TVector<T> operator/(const TVector<T>& lhs, const TVector<T>& rhs)
	{	return TVector<T>(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z); }
template<class T> inline
TVector<T> operator/(const TVector<T>& lhs, double rhs)
	{	return TVector<T>(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs); }





#define BOOL_TYPE  char

template<class T>
TVector2<BOOL_TYPE> operator>(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 2; i++)
		res[i] = lhs[i] > rhs[i];
	return res;
}

template<class T>
TVector2<BOOL_TYPE> operator<(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 2; i++)
		res[i] = lhs[i] < rhs[i];
	return res;
}

template<class T>
TVector2<BOOL_TYPE> operator||(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 2; i++)
		res[i] = lhs[i] || rhs[i];
	return res;
}

template<class T>
TVector2<BOOL_TYPE> operator&&(const TVector2<T>& lhs, const TVector2<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 2; i++)
		res[i] = lhs[i] && rhs[i];
	return res;
}

template<class T>
TVector<BOOL_TYPE> operator>(const TVector<T>& lhs, const TVector<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 3; i++)
		res[i] = lhs[i] > rhs[i];
	return res;
}

template<class T>
TVector<BOOL_TYPE> operator<(const TVector<T>& lhs, const TVector<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 3; i++)
		res[i] = lhs[i] < rhs[i];
	return res;
}

template<class T>
TVector<BOOL_TYPE> operator||(const TVector<T>& lhs, const TVector<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 3; i++)
		res[i] = lhs[i] || rhs[i];
	return res;
}

template<class T>
TVector<BOOL_TYPE> operator&&(const TVector<T>& lhs, const TVector<T>& rhs)
{
	TVector<BOOL_TYPE> res;
	for (int i = 0; i < 3; i++)
		res[i] = lhs[i] && rhs[i];
	return res;
}



template<class T>
ostream& operator<<(ostream& os, const TVector<T>& vect)
{
	os.width(10);
	os << vect[0] << ' ';
	os.width(10);
	os << vect[1] << ' ';
	os.width(10);
	os << vect[2];
	return os;
}



#endif // !defined(VECTOR_H)
