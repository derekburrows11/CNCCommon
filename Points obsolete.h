// Points.h
//
//////////////////////////////////////////////////////////////////////

#if !defined(POINTS_H__INCLUDED_)
#define POINTS_H__INCLUDED_


#include <math.h>


#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


#ifndef ASSERT
	#include <assert.h>			// use if not MFC
	#define ASSERT assert			// ASSERT is MFC
#endif



#define NEARINT(x) ((int)floor(x + 0.5))



//#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
//#pragma warning(disable: 4800)	// forcing value to bool 'true' or 'false' (performance warning)
//#pragma warning(disable: 4804)	// unsafe use of type 'bool' in operation


template<class T> class TPoint3;	// notification for TPoint2(const TPoint3<T>& pt3)

template<class T> class TPoint2
{
public:
	T x, y;
	TPoint2() {}
	TPoint2(const T ix, const T iy)	{ x = ix; y = iy; }
	TPoint2(const T val)					{ x = y = val; }
	TPoint2(const TPoint3<T>& pt3)	{ x = pt3.x; y = pt3.y; }
#ifdef __AFXWIN_H__
	TPoint2(const CPoint pt)	{ x = (T)pt.x; y = (T)pt.y; }
	TPoint2(const CSize  pt)	{ x = (T)pt.cx; y = (T)pt.cy; }
	operator CPoint() const		{ return CPoint(NEARINT(x), NEARINT(y)); }
#endif

	void Set(const T ix, const T iy) { x = ix; y = iy; }
	
	T& operator[](int idx) { ASSERT(idx >= 0 && idx < 2); return (&x)[idx]; }
	T operator[](int idx) const { ASSERT(idx >= 0 && idx < 2); return (&x)[idx]; }

	void Lowest(const TPoint2& other)
	{
		if (other.x < x) x = other.x;
		if (other.y < y) y = other.y;
	}
	void Highest(const TPoint2& other)
	{
		if (other.x > x) x = other.x;
		if (other.y > y) y = other.y;
	}

	TPoint2& operator+=(const TPoint2& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x += rhs.x;
		y += rhs.y;
		return *this;
	#pragma warning(default: 4244)
	}
	TPoint2& operator-=(const TPoint2& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	#pragma warning(default: 4244)
	}

	int operator==(const TPoint2& rhs) const
	{
		return x == rhs.x && y == rhs.y;
	}

	double MagSq() const { return x*x + y*y; }
	double Mag() const { return sqrt(x*x + y*y); }
	void Unit();

	void Rotate90()	{ T t = x; x = (T)-y; y = t; }	// 90 deg CCW
	void Rotate90CCW(){ T t = x; x = (T)-y; y = t; }	// 90 deg CCW
	void Rotate180()	{ x = (T)-x; y = (T)-y; }
	void Rotate270()	{ T t = x; x = y; y = (T)-t; }	// 90 deg CW
	void Rotate90CW()	{ T t = x; x = y; y = (T)-t; }	// 90 deg CW
};


template<class T> class TPoint3
{
public:
	T x, y, z;
	TPoint3() {}
	TPoint3(const T ix, const T iy, const T iz)	{ x = ix; y = iy, z = iz; }
	TPoint3(const T val)					{ x = y = z = val; }
	TPoint3(const TPoint2<T>& pt2)	{ x = pt2.x; y = pt2.y; z = 0; }
#ifdef __AFXWIN_H__
	TPoint3(const CPoint pt)			{ x = (T)pt.x; y = (T)pt.y; z = 0; }
	TPoint3(const CSize  pt)			{ x = (T)pt.cx; y = (T)pt.cy; z = 0; }
	operator CPoint() const				{ return CPoint(NEARINT(x), NEARINT(y)); }
#endif

//	operator TPoint2<T>() const { return TPoint2<T>(x, y); }
//	operator TVector<T>() const { return TVector<T>(x, y, z); }

	void Set(const T ix, const T iy, const T iz) { x = ix; y = iy, z = iz; }
	
	T& operator[](int idx) { ASSERT(idx >= 0 && idx < 3); return (&x)[idx]; }
	T operator[](int idx) const { ASSERT(idx >= 0 && idx < 3); return (&x)[idx]; }

	void Lowest(const TPoint3& other)
	{
		if (other.x < x) x = other.x;
		if (other.y < y) y = other.y;
		if (other.z < z) z = other.z;
	}
	void Highest(const TPoint3& other)
	{
		if (other.x > x) x = other.x;
		if (other.y > y) y = other.y;
		if (other.z > z) z = other.z;
	}
	void Lowest(const TPoint3& pt1, const TPoint3& pt2)
	{
		x = (pt1.x <= pt2.x) ? pt1.x : pt2.x;
		y = (pt1.y <= pt2.y) ? pt1.y : pt2.y;
		z = (pt1.z <= pt2.z) ? pt1.z : pt2.z;
	}
	void Highest(const TPoint3& pt1, const TPoint3& pt2)
	{
		x = (pt1.x >= pt2.x) ? pt1.x : pt2.x;
		y = (pt1.y >= pt2.y) ? pt1.y : pt2.y;
		z = (pt1.z >= pt2.z) ? pt1.z : pt2.z;
	}

	TPoint3& operator+=(const TPoint3& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	#pragma warning(default: 4244)
	}
	TPoint3& operator-=(const TPoint3& rhs)
	{
	#pragma warning(disable: 4244)	// conversion from 'int' to 'char', possible loss of data
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	#pragma warning(default: 4244)
	}

	template<class T2>
	TPoint3<T>& operator=(const TPoint2<T2>& pt2) { x = pt2.x; y = pt2.y; z = 0; return *this; }

	int operator==(const TPoint3& rhs) const
	{
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}

	double MagSq() const { return x*x + y*y + z*z; }
	double Mag() const { return sqrt(x*x + y*y + z*z); }
	void Unit();


};


typedef TPoint2<float>  CPoint2f;
typedef TPoint2<double> CPoint2d;

typedef TPoint3<float>  CPoint3f;
typedef TPoint3<double> CPoint3d;

/////////////////////////////////////
// Member functions
/////////////////////////////////////

template<class T>
void TPoint2<T>::Unit()
{ double ivMag = 1/sqrt(x*x+y*y); x*=ivMag; y*=ivMag; }

template<class T>
void TPoint3<T>::Unit()
{ double ivMag = 1/sqrt(x*x+y*y+z*z); x*=ivMag; y*=ivMag; z*=ivMag; }



/////////////////////////////////////
// TPoint associated functions
/////////////////////////////////////

template<class T>
double cross(TPoint2<T>& lhs, TPoint2<T>& rhs)
{
	return (lhs.x * rhs.y) - (lhs.y * rhs.x);
}

template<class T>
double dot(TPoint2<T>& lhs, TPoint2<T>& rhs)
{
	return (lhs.x * rhs.x) + (lhs.y * rhs.y);
}

/////////////////////////////////////
// TPoint associated operators
/////////////////////////////////////

// TPoint2 operators

template<class T>
TPoint2<T> operator+(const TPoint2<T>& lhs, const TPoint2<T>& rhs)
{
	return TPoint2<T>(lhs.x + rhs.x, lhs.y + rhs.y);
}

template<class T>
TPoint2<T> operator-(const TPoint2<T>& lhs, const TPoint2<T>& rhs)
{
	return TPoint2<T>(lhs.x - rhs.x, lhs.y - rhs.y);
}

template<class T>
TPoint2<T> operator*(const TPoint2<T>& lhs, const TPoint2<T>& rhs)
{
	return TPoint2<T>(lhs.x * rhs.x, lhs.y * rhs.y);
}
template<class T>
TPoint2<T> operator*(const TPoint2<T>& lhs, double rhs)
{
	return TPoint2<T>(lhs.x * rhs, lhs.y * rhs);
}

#ifdef __AFXWIN_H__
template<class T>		// use POINT not CPoint - CPoint(DWORD pt) enables this to be used with int's!
TPoint2<T> operator*(const TPoint2<T>& lhs, const POINT rhs)
{
	return TPoint2<T>(lhs.x * rhs.x, lhs.y * rhs.y);
}
template<class T>		// use SIZE not CSize - CSize(DWORD pt) enables this to be used with int's!
TPoint2<T> operator*(const TPoint2<T>& lhs, const SIZE rhs)
{
	return TPoint2<T>(lhs.x * rhs.cx, lhs.y * rhs.cy);
}
#endif

template<class T>
TPoint2<T> operator/(const TPoint2<T>& lhs, const TPoint2<T>& rhs)
{
	return TPoint2<T>(lhs.x / rhs.x, lhs.y / rhs.y);
}
template<class T>
TPoint2<T> operator/(const TPoint2<T>& lhs, double rhs)
{
	return TPoint2<T>(lhs.x / rhs, lhs.y / rhs);
}



// TPoint3 operators

template<class T>
TPoint3<T> operator+(const TPoint3<T>& lhs, const TPoint3<T>& rhs)
{
	return TPoint3<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

template<class T>
TPoint3<T> operator-(const TPoint3<T>& lhs, const TPoint3<T>& rhs)
{
	return TPoint3<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

template<class T>
TPoint3<T> operator*(const TPoint3<T>& lhs, const TPoint3<T>& rhs)
{
	return TPoint3<T>(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}
template<class T>
TPoint3<T> operator*(const TPoint3<T>& lhs, double rhs)
{
	return TPoint3<T>(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

template<class T>
TPoint3<T> operator/(const TPoint3<T>& lhs, const TPoint3<T>& rhs)
{
	return TPoint3<T>(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);
}
template<class T>
TPoint3<T> operator/(const TPoint3<T>& lhs, double rhs)
{
	return TPoint3<T>(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}





//#pragma warning(default: 4244)
//#pragma warning(default: 4800)
//#pragma warning(default: 4804)


#endif // !defined(POINTS_H__INCLUDED_)
