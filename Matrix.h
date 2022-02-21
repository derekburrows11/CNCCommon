//////////////////////////////////////////////////////////////////////
// Matrix.h: interface for the CMatrix and CVector classes.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__48FCACC3_FE73_11D4_8C1E_9D2D1A323C2C__INCLUDED_)
#define AFX_MATRIX_H__48FCACC3_FE73_11D4_8C1E_9D2D1A323C2C__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000



#include "Vector.h"

#include <strstrea.h>




//////////////////////////////////////////////
// template TMatrix<T>
//////////////////////////////////////////////
/*
template<class T> struct TMatrixBase
{
	T* arVal;
};

template<class T> struct TVectorBase
{
	T arVal[3];
	union
	{
		T arVal[3];
		T x, y, z;
	};
};
*/


//template<class BASE, class T> class TMatrix : public BASE
template<class T> class TMatrix
{
public:
	int w, h;
	T* arVal;		// matrix array is row by row
// constructors
public:
	TMatrix()
		{ w = h = 0; arVal = NULL; }
	explicit TMatrix(int vectorLength);		// column vector
	TMatrix(int rows, int columns);
		template<class T2>
	TMatrix(const TMatrix<T2>& matrix)
	{
		w = matrix.w;	h = matrix.h;
		int elements = h * w;
		arVal = new T[elements];
		for (int i = 0; i < elements; i++)
			arVal[i] = (T)matrix.arVal[i];
	}
	TMatrix(const TMatrix<T>& matrix);
	TMatrix(const TVector<T>& vect);
	TMatrix(const TVector<T>& vect, int vectorLength);

// destructor
	virtual ~TMatrix() { if (arVal != NULL) delete[] arVal; }

// member functions
public:
	virtual void SetSize(int rows, int columns);
	virtual void SetSize(int vectorLength) { SetSize(vectorLength, 1); }		// default column vector
		template<class T2>
	void SetSize(const TMatrix<T2>& matrix) { SetSize(matrix.h, matrix.w); }
	T* GetArray() { return arVal; }
	int Width() const { return w; }
	int Height() const { return h; }
	int Length() const { return __max(w, h); }
	bool IsVector() const { return w == 1 || h == 1; }
	void Transpose();								// transpose this matrix
	void Transpose(TMatrix<T>& matrixT) const;	// matrixT = transpose
	void ResizeToVectorByRows() { h = h * w; w = 1; }
	void ResizeToVectorByColumns();
	void Neg();

	void CopyToDestLoc(TMatrix<T>& mxDest, int destRow, int destCol) const;
	void CopyFromSrcLoc(const TMatrix<T>& mxSrc, int srcRow,  int srcCol);
	void CopyToLocFromSrcLoc(const CPoint& ptDest, const TMatrix<T>& mxSrc, const CPoint& ptSrc, const CSize& size);

		template<class T2>
	void CopyToArray(T2* arDest) { CopyToArrayByRow(arDest); };
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void CopyToArrayByRow(T2* arDest)		// row by row
	{
		int elements = h * w;
		for (int i = 0; i < elements; i++)
			arDest[i] = arVal[i];
	}
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void CopyToArrayByCol(T2* arDest)		// column by column
	{
		int elements = h * w;
		int i = 0;
		for (int col = 0; col < w; col++)
			for (int mxIdx = col; mxIdx < elements; mxIdx += w)
				arDest[i++] = arVal[mxIdx];
	}

		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void SetFromArray(const T2* arSrc)		// copies into vector
	{
		if (!IsVector())
			TRACE0("TMatrix::SetFromArray() matrix is not a vector!\n");
		int elements = h * w;
		for (int i = 0; i < elements; i++)
			arVal[i] = (T)arSrc[i];
	}
	
	// where CopyFromArrayBy....()
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void SetFromArrayByRows(const T2* arSrc)		// copies into matrix row by row
	{
		int elements = h * w;
		for (int i = 0; i < elements; i++)
			arVal[i] = (T)arSrc[i];
	}
		template<class T2>			// multi class templates need to be defined in class definition apparently!
	void SetFromArrayByCols(const T2* arSrc)		// copies into matrix column by column
	{
		int elements = h * w;
		int i = 0;
		for (int col = 0; col < w; col++)
			for (int mxIdx = col; mxIdx < elements; mxIdx += w)
				arVal[mxIdx] = (T)arSrc[i++];
	}
	
	void Prod(const TMatrix<T>& lhs, const TMatrix<T>& rhs);
	void Prod(const TMatrix<T>& lhs, const T& rhs);
	void ProdTr(const TMatrix<T>& lhs, const TMatrix<T>& rhsT);
	void ProdTl(const TMatrix<T>& lhsT, const TMatrix<T>& rhs);
	bool ProdPart(const TMatrix<T>& lhs, const TMatrix<T>& rhs);

	// logical functions
	inline bool Any() const;
	bool All() const;
	int Find(TMatrix<int>& index) const;		// sets to indices of all non-zero elements

	// vector maths
	double Mag() const { return sqrt(MagSq()); }
	double MagSq() const;
	double Sum() const;
	double SumAbs() const;
	T MaxElem() const;
	T MinElem() const;
	T MaxAbsElem() const;
	T MinAbsElem() const;
	void Unit(TMatrix<T>& matrixU) const;

	// set matrix types
	void Identity();
	void RotateZ(double theta);
	void RotateY(double theta);
	void RotateX(double theta);
	void RotateAzEl(double az, double el);
//	int ToCPoints(CPoint* arPt, int iMax = -1) const;
	void SetRow(int row, const TMatrix<T>& vect);
	void SetColumn(int col, const TMatrix<T>& vect);
	void SetRow(int row, const TVector<T>& vect);
	void SetColumn(int col, const TVector<T>& vect);
	void SetRow(int row, const T& val);
	void SetColumn(int col, const T& val);
	void SetRow(int row, const T* arSrc);
	void SetColumn(int col, const T* arSrc);
	void ScaleRowsBy(const TMatrix<T>& vt);
	void ScaleRowsBy(const TVector2<T>& vt);
	void ScaleRowsBy(const TVector<T>& vt);
	void ScaleColumnsBy(const TMatrix<T>& vt);
	void ScaleColumnsBy(const TVector2<T>& vt);
	void ScaleColumnsBy(const TVector<T>& vt);
	void SwapRows(int row1, int row2);
	void SwapColumns(int col1, int col2);
	void CopyRow(int rowDest, int rowSrc, int num = 1);
	void CopyColumn(int colDest, int colSrc, int num = 1);
	void BlockScale(const TMatrix<T>& mxBlock, const TMatrix<T>& mxScaler);



// polynomial functions
	void Power3Vector(double s);
	void Power2Vector(double s);
	void Power1Vector(double s);
	void Power3Matrix(const TMatrix<T>& vtS);
	void Power2Matrix(const TMatrix<T>& vtS);
	void Power1Matrix(const TMatrix<T>& vtS);
	void SolvePolyAt(const TMatrix<T>& vectS, TMatrix<T>& vectRes) const;
	void SolvePolyAt(const TMatrix<T>& vectS, TVector<T>& vectRes) const;
	void SolvePolyAt0(TVector<T>& vectRes) const;
	void SolvePolyAt1(TVector<T>& vectRes) const;
	void Bezier2Poly(TMatrix<T>& poly) const;
	void BezierRel2Poly(TMatrix<T>& poly) const;
	void Poly2Bezier(TMatrix<T>& bezier) const;
	void Poly2BezierRel(TMatrix<T>& bezierRel) const;
	void Derivative(TMatrix<T>& deriv) const;
	void Intergral(TMatrix<T>& inter) const;
	
	T& elem(int row, int col) { ASSERT(row<h && col<w); return arVal[row*w + col]; }
	T elem(int row, int col) const { ASSERT(row<h && col<w); return arVal[row*w + col]; }

	// Inverse functions
	double Invert();											// Invert this matrix, returns determinant
	double Inverse(TMatrix<T>& matrixI) const;		// matrixI = inverse, returns determinant
	T Determinant() const;

	// the following LU definitions are in "PolyFunc.cpp" - only defined for <double>
	int LUSolve(double* vtX, const double* vtVal) const;		// original version used - from Kreyszig
	int LUSolve2(double* vtX, const double* vtVal) const;		// downloaded from internet
	int LUInvert();														// downloaded from internet
	int LUInverse(TMatrix<T>& matrixI) const;						// downloaded from internet
	void LUDecompose(int indx[], bool parity = true);			// downloaded from internet
	void LUBackSubstitute(double* b, int indx[]) const;		// downloaded from internet

	


	T& operator[](int idx) { ASSERT(idx < w*h); return arVal[idx]; }
	T operator[](int idx) const { ASSERT(idx < w*h); return arVal[idx]; }

	operator T*() { return arVal; }	// gives overload ambiguity!

		template<class T2>
	TMatrix<T>& operator=(const TMatrix<T2>& matrix)	// multi class templates need to be defined in class definition apparently!
	{
		SetSize(matrix);
		int elements = h * w;
		for (int i = 0; i < elements; i++)
			arVal[i] = (T)matrix[i];
		return *this;
	}
	TMatrix& operator=(const TMatrix<T>& matrix);
	TMatrix<T>& operator=(const T val);

	TMatrix<T>& operator+=(const TMatrix<T>& rhs);
	TMatrix<T>& operator-=(const TMatrix<T>& rhs);
	TMatrix<T>& operator+=(const T& rhs);
	TMatrix<T>& operator-=(const T& rhs);
	TMatrix<T>& operator*=(const T& rhs);
	TMatrix<T>& operator/=(const T& rhs);

	TMatrix<T> operator-() const;


};




///////////////////////////////////////////
///////////////////////////////////////////

//typedef TMatrix<TMatrixBase<float>, float> CMatrixf;
//typedef TMatrix<TMatrixBase<double>, double> CMatrixd;
//typedef TMatrix<TMatrixBase<double>, double> CMatrix;
typedef TMatrix<float> CMatrixf;
typedef TMatrix<double> CMatrix;


//////////////////////////////////////////////
// template TMatrixWrap<T>
//////////////////////////////////////////////
// enables a TVector to be used as a TMatrix
template<class T> class TMatrixWrap : public TMatrix<T>
{
public:
// constructors
	TMatrixWrap(TVector<T>& vect) { arVal = vect.GetArray(); h = 3; w = 1; }
	TMatrixWrap(TVector<T>& vect, int length)
		{ arVal = vect.GetArray(); ASSERT(length <= 3); h = length; w = 1; }
	TMatrixWrap(T* array, int rows, int columns)
		{ arVal = array; h = rows; w = columns; }

// destructor
	virtual ~TMatrixWrap() { arVal = NULL; }		// also try to stop call to base destructor
};

typedef TMatrixWrap<float> CMatrixWrapf;
typedef TMatrixWrap<double> CMatrixWrap;



///////////////////////////////////////////
// TMatrix member functions
///////////////////////////////////////////

// Constructors

template<class T>
TMatrix<T>::TMatrix(int rows, int columns)
{
	h = rows;
	w = columns;
	arVal = new T[w * h];
}

template<class T>
TMatrix<T>::TMatrix(int vectorLength)		// column vector
{
	h = vectorLength;
	w = 1;
	arVal = new T[h];
}

template<class T>
TMatrix<T>::TMatrix(const TMatrix<T>& matrix)
{
	h = matrix.h;	w = matrix.w;
	int elements = h * w;
	arVal = new T[elements];
	memcpy(arVal, matrix.arVal, sizeof(T) * elements);
}

template<class T>
TMatrix<T>::TMatrix(const TVector<T>& vect)
{
	h = 3;
	w = 1;
	arVal = new T[3];
	arVal[0] = vect[0];
	arVal[1] = vect[1];
	arVal[2] = vect[2];
}

template<class T>
TMatrix<T>::TMatrix(const TVector<T>& vect, int vectorLength)
{
	h = vectorLength;
	w = 1;
	arVal = new T[h];
	int vectMax = __min(3, vectorLength);
	for (int i = 0; i < vectMax; i++)
		arVal[i] = vect[i];
	while (i < vectorLength)
		arVal[i] = 0;				// fill rest with 0's
}


// functions

template<class T>
void TMatrix<T>::SetSize(int rows, int columns)
{
	if (rows * columns != h * w)
	{
		if (arVal) delete[] arVal;
		if (rows == 0 || columns == 0)
		{
			w = h = 0; 
			arVal = NULL;
			return;
		}
		arVal = new T[rows * columns];
	}
	h = rows;
	w = columns;
}


/*
template<class T, class T2>	// multi class templates need to be defined in class definition apparently!
TMatrix<T>& TMatrix<T>::operator=(const TMatrix<T2>& matrix)
{
	SetSize(matrix);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] = matrix[i];
	return *this;
}
*/

template<class T>
TMatrix<T>& TMatrix<T>::operator=(const TMatrix<T>& matrix)
{
	SetSize(matrix);
	memcpy(arVal, matrix.arVal, sizeof(T) * w * h);
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator=(const T val)
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] = val;
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator+=(const TMatrix<T>& rhs)
{
	ASSERT(w == rhs.w && h == rhs.h);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] += rhs[i];
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator+=(const T& rhs)
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] += rhs;
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator-=(const TMatrix<T>& rhs)
{
	ASSERT(w == rhs.w && h == rhs.h);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] -= rhs[i];
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator-=(const T& rhs)
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] -= rhs;
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator*=(const T& rhs)
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] *= rhs;
	return *this;
}

template<class T>
TMatrix<T>& TMatrix<T>::operator/=(const T& rhs)
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] /= rhs;
	return *this;
}

template<class T>
void TMatrix<T>::Identity()
{
	ASSERT(w == h);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] = 0;
	for (i = 0; i < elements; i += w+1)
		arVal[i] = 1;
}

template<class T>
bool TMatrix<T>::Any() const
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		if (arVal[i] != 0) return true;
	return false;
}


#define BOOL_TYPE  char

bool TMatrix<BOOL_TYPE>::Any() const
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		if (arVal[i] != false) return true;
	return false;
}

template<class T>
bool TMatrix<T>::All() const
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		if (arVal[i] == 0) return false;
	return true;
}

template<class T>
double TMatrix<T>::MagSq() const
{
	int elements = h * w;
	double sumSq = 0;
	for (int i = 0; i < elements; i++)
		sumSq += arVal[i] * arVal[i];
	return sumSq;
}

template<class T>
double TMatrix<T>::Sum() const
{
	int elements = h * w;
	double sum = 0;
	for (int i = 0; i < elements; i++)
		sum += arVal[i];
	return sum;
}

template<class T>
double TMatrix<T>::SumAbs() const
{
	int elements = h * w;
	double sum = 0;
	for (int i = 0; i < elements; i++)
		sum += fabs(arVal[i]);
	return sum;
}

template<class T>
T TMatrix<T>::MaxElem() const
{
	int elements = h * w;
	T maxElem = arVal[0];
	for (int i = 1; i < elements; i++)
		if (maxElem < arVal[i])
			maxElem = arVal[i];
	return maxElem;
}

template<class T>
T TMatrix<T>::MinElem() const
{
	int elements = h * w;
	T minElem = arVal[0];
	for (int i = 1; i < elements; i++)
		if (minElem > arVal[i])
			minElem = arVal[i];
	return minElem;
}

template<class T>
T TMatrix<T>::MaxAbsElem() const
{
	int elements = h * w;
	T maxElem = fabs(arVal[0]);
	for (int i = 1; i < elements; i++)
		if (maxElem < fabs(arVal[i]))
			maxElem = fabs(arVal[i]);
	return maxElem;
}

template<class T>
T TMatrix<T>::MinAbsElem() const
{
	int elements = h * w;
	T minElem = fabs(arVal[0]);
	for (int i = 1; i < elements; i++)
		if (minElem > fabs(arVal[i]))
			minElem = fabs(arVal[i]);
	return minElem;
}

template<class T>
void TMatrix<T>::Unit(TMatrix<T>& matrixU) const
{
	double oneOnMag = 1 / Mag();		// do division only once!
	matrixU.SetSize(h, w);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		matrixU[i] = arVal[i] * oneOnMag;
}

template<class T>
void TMatrix<T>::Neg()
{
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		arVal[i] = -arVal[i];
}

template<class T>
TMatrix<T> TMatrix<T>::operator-() const
{
	TMatrix<T> res(h, w);
	int elements = h * w;
	for (int i = 0; i < elements; i++)
		res[i] = -arVal[i];
	return res;
}

template<class T>
int TMatrix<T>::Find(TMatrix<int>& index) const
{
	// sets to indices of all non-zero elements
	int elements = h * w;
	int count = 0;
	for (int i = 0; i < elements; i++)
		if (arVal[i] != 0)
			count++;
	index.SetSize(count);
	count = 0;
	for (int i = 0; i < elements; i++)
		if (arVal[i] != 0)
			index.arVal[count++] = i;
	return count;
}


template<class T>
void TMatrix<T>::RotateZ(double theta)
{
	if (w<3 || h<3)
		SetSize(3, 3);
	arVal[0] = arVal[w + 1] = (T)cos(theta);
	arVal[w + 0] = -(arVal[1] = (T)sin(theta));
	arVal[2*w + 2] = 1;
	arVal[2] = arVal[w + 2] = arVal[2*w + 0] = arVal[2*w + 1] = 0;
}
template<class T>
void TMatrix<T>::RotateY(double theta)
{
	if (w<3 || h<3)
		SetSize(3, 3);
	arVal[0] = arVal[2*w + 2] = (T)cos(theta);
	arVal[2] = -(arVal[2*w + 0] = (T)sin(theta));
	arVal[w + 1] = 1;
	arVal[1] = arVal[w + 0] = arVal[w + 2] = arVal[2*w + 1] = 0;
}
template<class T>
void TMatrix<T>::RotateX(double theta)
{
	if (w<3 || h<3)
		SetSize(3, 3);
	arVal[w + 1] = arVal[2*w + 2] = (T)cos(theta);
	arVal[2*w + 1] = -(arVal[w + 2] = (T)sin(theta));
	arVal[0] = 1;
	arVal[1] = arVal[2] = arVal[w + 0] = arVal[2*w + 0] = 0;
}

/*
void TMatrix<T>::RotateAzEl(double az, double el)
% 	VIEW(AZ,EL) and VIEW([AZ,EL]) set the angle of the view from which an
%	observer sees the current 3-D plot.  AZ is the azimuth or horizontal 
%	rotation and EL is the vertical elevation (both in degrees). Azimuth 
%	revolves about the z-axis, with positive values indicating counter-
%	clockwise rotation of the viewpoint. Positive values of elevation 
%	correspond to moving above the object; negative values move below.
%	VIEW([X Y Z]) sets the view angle in cartesian coordinates. The
%	magnitude of vector X,Y,Z is ignored.
% 
% 	Here are some examples:
% 
% 	AZ = -37.5, EL = 30 is the default 3-D view.
% 	AZ = 0, EL = 90 is directly overhead and the default 2-D view.
% 	AZ = EL = 0 looks directly up the first column of the matrix.
% 	AZ = 180 is behind the matrix.

	if a unit direction vector [x y z] convert to az & el
	az = atan2(x, -y) * 180/pi;
	el = atan2(z, sqrt(x^2 + y^2)) * 180/pi;
	ATAN2	Four quadrant inverse tangent.
	ATAN2(Y,X) is the four quadrant arctangent of the real parts of the
	elements of X and Y.  -pi <= ATAN2(Y,X) <= pi.
 
	View transformation matrix formed by composing two rotations:
		1) Rotate about the z axis -AZ radians			= rotZ(AZ)
		2) Rotate about the x axis (EL-pi/2) radians	= rotX(pi/2-EL)
T = [  cos(az)           sin(az)           0       ]
    [ -sin(el)*sin(az)   sin(el)*cos(az)   cos(el) ]
    [  cos(el)*sin(az)  -cos(el)*cos(az)   sin(el) ]
*/
template<class T>
void TMatrix<T>::RotateAzEl(double az, double el)
{
	if (w<3 || h<3)
		SetSize(3, 3);
	double sinaz = sin(az * deg2rad);
	double cosaz = cos(az * deg2rad);
	double sinel = sin(el * deg2rad);
	double cosel = cos(el * deg2rad);
	arVal[0] = (T)cosaz;
	arVal[1] = (T)sinaz;
	arVal[2] = 0;
	arVal[w + 0] = (T)(-sinaz * sinel);
	arVal[w + 1] = (T)(cosaz * sinel);
	arVal[w + 2] = (T)cosel;
	arVal[2*w + 0] = (T)(sinaz * cosel);
	arVal[2*w + 1] = (T)(-cosaz * cosel);
	arVal[2*w + 2] = (T)sinel;
}

template<class T>
void TMatrix<T>::Transpose(TMatrix<T>& matrixT) const
{
	matrixT.SetSize(w, h);
	int idx, idxT = 0;
	for (int col = 0; col < w; col++)
	{
		idx = col;
		for (int row = 0; row < h; row++)
		{
			matrixT.arVal[idxT++] = arVal[idx];
			idx += w;
		}
	}
}

template<class T>
void TMatrix<T>::Transpose()
{
	int idx, idxT;
	T tmpVal;
	if (w != 1 && h != 1)		// if not a vector
		if (h == w)					// if square
			for (int col = 1; col < w; col++)
			{
				idx = col;
				idxT = col * w;
				for (int row = 0; row < col; row++)
				{
					tmpVal = arVal[idx];
					arVal[idx] = arVal[idxT];
					arVal[idxT] = tmpVal;
					idxT++;
					idx += w;
				}
			}
		else		// general rectangular matrix
		{
			CMatrix mxT;
			Transpose(mxT);
			*this = mxT;
			return;
			ASSERT(0);		// not finished
			idxT = 0;
			for (int col = 0; col < w; col++)
			{
				idx = col;
				for (int row = 0; row < h; row++)
				{
					tmpVal = arVal[idx];
					arVal[idx] = arVal[idxT];
					arVal[idxT] = tmpVal;
					idxT++;
					idx += w;
				}
			}
		}

	idx = w;
	w = h;
	h = idx;
}

template<class T>
void TMatrix<T>::ResizeToVectorByColumns()
{
	int elements = h * w;
	ASSERT(0);		// not done yet!
}


template<class T>
void TMatrix<T>::CopyToDestLoc(TMatrix<T>& mxDest, int destRow, int destCol) const
{
// Copies a part of the matrix from all/part of the given source matrix to origRow/Col of this matrix
	int rowI, rowF, colI, colF, rowS, colSI, idx, idxS;
	rowI = destRow;
	colI = destCol;
	rowF = __min(mxDest.h, rowI + h);
	colF = __min(mxDest.w, colI + w);
	rowS = 0;				// row of source matrix
	colSI = 0;				// initial column of source matrix
	for (int row = rowI; row < rowF; row++, rowS++)
	{
		idx = row*mxDest.w + colI;
		idxS = rowS*w + colSI;
		for (int col = colI; col < colF; col++)
			mxDest[idx++] = arVal[idxS++];
	}
}

template<class T>
void TMatrix<T>::CopyFromSrcLoc(const TMatrix<T>& mxSrc, int srcRow, int srcCol)
{
// Copies the matrix from the given location of the source matrix
	int rowI, rowF, colI, colF, rowS, colSI, idx, idxS;
	rowI = 0;
	colI = 0;
	rowF = __min(h, mxSrc.h-srcRow);
	colF = __min(w, mxSrc.w-srcCol);
	rowS = srcRow;				// row of source matrix
	colSI = srcCol;				// initial column of source matrix
	for (int row = rowI; row < rowF; row++, rowS++)
	{
		idx = row*w + colI;
		idxS = rowS*mxSrc.w + colSI;
		for (int col = colI; col < colF; col++)
			arVal[idx++] = mxSrc[idxS++];
	}
}

template<class T>
void TMatrix<T>::CopyToLocFromSrcLoc(const CPoint& ptDest, const TMatrix<T>& mxSrc, const CPoint& ptSrc, const CSize& size)
{
// Copies to ptDest matrix from the given location and size of the source matrix
	int rowI, rowF, colI, colF, rowS, colSI, idx, idxS;
	rowI = ptDest.y;
	colI = ptDest.x;
	rowF = __min(h, ptDest.y+mxSrc.h-ptSrc.y);
	if (rowF >= ptDest.y+size.cy)
		rowF = ptDest.y+size.cy;
	else
		TRACE("Can't fit size rows in TMatrix::CopyToLocFromSrcLoc()\n");
	colF = __min(w, ptDest.x+mxSrc.w-ptSrc.x);
	if (colF >= ptDest.x+size.cx)
		colF = ptDest.x+size.cx;
	else
		TRACE("Can't fit size columns in TMatrix::CopyToLocFromSrcLoc()\n");

	rowS = ptSrc.y;				// row of source matrix
	colSI = ptSrc.x;				// initial column of source matrix
	for (int row = rowI; row < rowF; row++, rowS++)
	{
		idx = row*w + colI;
		idxS = rowS*mxSrc.w + colSI;
		for (int col = colI; col < colF; col++)
			arVal[idx++] = mxSrc[idxS++];
	}
}


// Basic matrix operations

template<class T>
void TMatrix<T>::Prod(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.h);
	SetSize(lhs.h, rhs.w);
	int inner = lhs.w;
	T elem;
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
		{
			elem = 0;
			for (int i = 0; i < inner; i++)
				elem += lhs[row*inner + i] * rhs[i*w + col];
			arVal[row*w + col] = elem;
		}
}

template<class T>
void TMatrix<T>::Prod(const TMatrix<T>& lhs, const T& rhs)
{
	SetSize(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		arVal[i] = (T)(lhs[i] * rhs);
}

// ProdTr = lhs * Transpose(rhsT)
template<class T>
void TMatrix<T>::ProdTr(const TMatrix<T>& lhs, const TMatrix<T>& rhsT)
{
	ASSERT(lhs.w == rhsT.w);
	SetSize(lhs.h, rhsT.h);
	int inner = lhs.w;
	T elem;
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
		{
			elem = 0;
			for (int i = 0; i < inner; i++)
				elem += lhs[row*inner + i] * rhsT[col*inner + i];
			arVal[row*w + col] = elem;
		}
}
// ProdTl = Transpose(lhsT) * rhs
template<class T>
void TMatrix<T>::ProdTl(const TMatrix<T>& lhsT, const TMatrix<T>& rhs)
{
	ASSERT(lhsT.h == rhs.h);
	SetSize(lhsT.w, rhs.w);
	int inner = lhsT.h;
	T elem;
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
		{
			elem = 0;
			for (int i = 0; i < inner; i++)
				elem += lhsT[i*h + row] * rhs[i*w + col];
			arVal[row*w + col] = elem;
		}
}

template<class T>
bool TMatrix<T>::ProdPart(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
/*	Inner product.  If matrix is too small for result returns false but all
	applicable elements are stored
*/
	ASSERT(lhs.w == rhs.h);
	int inner = lhs.w;
	bool bRes = true;
	if (h < lhs.h) bRes = false;
	if (w < rhs.w) bRes = false;
	int numRows = __min(h, lhs.h);
	int numCols = __min(w, rhs.w);
	T elem;
	for (int row = 0; row < numRows; row++)
		for (int col = 0; col < numCols; col++)
		{
			elem = 0;
			for (int i = 0; i < inner; i++)
				elem += lhs[row*lhs.w + i] * rhs[i*rhs.w + col];
			arVal[row*w + col] = elem;
		}
	return bRes;
}


template<class T>
double TMatrix<T>::Inverse(TMatrix<T>& matrixI) const
{
	ASSERT(w == h);		// must be square
	if (w != h)
		return 0;			// matrix not square!
	int n = w;				// order of matrix

	matrixI.SetSize(w, h);
	const T* const a = arVal;
	T* const ai = matrixI.arVal;
	double det, invDet;
	if (n == 2)
	{
		det = a[0]*a[3] - a[1]*a[2];
		if (det == 0)
			return det;
		invDet = 1.0 / det;
		ai[0] =  a[3] * invDet;
		ai[3] =  a[0] * invDet;
		ai[1] = -a[1] * invDet;
		ai[2] = -a[2] * invDet;
	}
	else if (n == 3)		// 3 x 3 inverse from Kreysig p.411
	{
		det  = a[0] * (a[4]*a[8] - a[7]*a[5]);
		det -= a[3] * (a[1]*a[8] - a[7]*a[2]);
		det += a[6] * (a[1]*a[5] - a[4]*a[2]);
		if (det == 0)
			return det;
		invDet = 1.0 / det;
		ai[0] = invDet *  (a[4]*a[8] - a[7]*a[5]);		// A11
		ai[1] = invDet * -(a[1]*a[8] - a[7]*a[2]);		// A21
		ai[2] = invDet *  (a[1]*a[5] - a[4]*a[2]);		// A31
		ai[3] = invDet * -(a[3]*a[8] - a[6]*a[5]);		// A12
		ai[4] = invDet *  (a[0]*a[8] - a[6]*a[2]);		// A22
		ai[5] = invDet * -(a[0]*a[5] - a[3]*a[2]);		// A32
		ai[6] = invDet *  (a[3]*a[7] - a[6]*a[4]);		// A13
		ai[7] = invDet * -(a[0]*a[7] - a[6]*a[1]);		// A23
		ai[8] = invDet *  (a[0]*a[4] - a[3]*a[1]);		// A33
	}
	else
	{
		ASSERT(0);		// not handled yet!
		return 0;
	}
	return det;
}

template<class T>
double TMatrix<T>::Invert()
{
	ASSERT(w == h);		// must be square
	if (w != h)
		return 0;			// matrix not square!
	int n = w;				// order of matrix

	T* const a = arVal;
	double det, invDet;
	if (n == 2)
	{
		det = a[0]*a[3] - a[1]*a[2];
		if (det == 0)
			return det;
		invDet = 1.0 / det;
		T tmp = a[0];
		a[0] =  a[3] * invDet;
		a[3] =   tmp * invDet;
		a[1] *= -invDet;
		a[2] *= -invDet;
	}
	else if (n == 3)		// 3 x 3 inverse from Kreysig p.411
	{
		T ai[9];		// for temp inverse storage
		det  = a[0] * (a[4]*a[8] - a[7]*a[5]);
		det -= a[3] * (a[1]*a[8] - a[7]*a[2]);
		det += a[6] * (a[1]*a[5] - a[4]*a[2]);
		if (det == 0)
			return det;
		invDet = 1.0 / det;
		ai[0] = invDet *  (a[4]*a[8] - a[7]*a[5]);		// A11
		ai[1] = invDet * -(a[1]*a[8] - a[7]*a[2]);		// A21
		ai[2] = invDet *  (a[1]*a[5] - a[4]*a[2]);		// A31
		ai[3] = invDet * -(a[3]*a[8] - a[6]*a[5]);		// A12
		ai[4] = invDet *  (a[0]*a[8] - a[6]*a[2]);		// A22
		ai[5] = invDet * -(a[0]*a[5] - a[3]*a[2]);		// A32
		ai[6] = invDet *  (a[3]*a[7] - a[6]*a[4]);		// A13
		ai[7] = invDet * -(a[0]*a[7] - a[6]*a[1]);		// A23
		ai[8] = invDet *  (a[0]*a[4] - a[3]*a[1]);		// A33
		for (int i = 0; i < 9; i++)
			a[i] = ai[i];
	}
	else
	{
		ASSERT(0);		// not handled yet!
		return 0;
	}
	return det;
}

template<class T>
T TMatrix<T>::Determinant() const
{
	ASSERT(w == h);		// must be square
	if (w != h)
		return 0;			// matrix not square!
	int n = w;				// order of matrix

	const T* const a = arVal;
	double det;
	if (n == 2)
		det = a[0]*a[3] - a[1]*a[2];
	else if (n == 3)		// 3 x 3 determinant from Kreysig p.396
	{
		det  = a[0] * (a[4]*a[8] - a[7]*a[5]);
		det -= a[3] * (a[1]*a[8] - a[7]*a[2]);
		det += a[6] * (a[1]*a[5] - a[4]*a[2]);
	}
	else
	{
		ASSERT(0);		// not handled yet!
		det = 0;
	}
	return det;
}

/*
template<class T>
int TMatrix<T>::ToCPoints(CPoint* arPt, int iMax) const
{
	ASSERT(w>=1 && w<=3);
	if (iMax == -1 || iMax > h)		// if iMax == -1 (default) does all rows
		iMax = h;
	int mxIdx = 0;
	for (int i = 0; i < iMax; i++)
	{
		arPt[i].x = (int)arVal[mxIdx + 0];
		arPt[i].y = (int)arVal[mxIdx + 1];
		mxIdx += w;
	}
	return iMax;
}
*/

template<class T>
void TMatrix<T>::Power3Vector(double s)	// for cubics
{
	SetSize(1, 4);		// make it a row vector
	double sPow = s;
	arVal[0] = 1;
	arVal[1] = (T)sPow;
	sPow *= s;
	arVal[2] = (T)sPow;
	sPow *= s;
	arVal[3] = (T)sPow;
}

template<class T>
void TMatrix<T>::Power2Vector(double s)	// for quadratics
{
	SetSize(1, 3);		// make it a row vector
	arVal[0] = 1;
	arVal[1] = (T)s;
	arVal[2] = (T)(s*s);
}

template<class T>
void TMatrix<T>::Power1Vector(double s)	// for linear
{
	SetSize(1, 2);		// make it a row vector
	arVal[0] = 1;
	arVal[1] = (T)s;
}

template<class T>
void TMatrix<T>::Power3Matrix(const TMatrix<T>& vtS)	// for cubics
{
	ASSERT(vtS.IsVector());
	SetSize(vtS.Length(), 4);		// make columns of vtS to power 0, 1, 2, 3
	int i = 0;
	for (int row = 0; row < h; row++)
	{
		double s = vtS[row];
		double sPow = s;
		arVal[i++] = 1;
		arVal[i++] = (T)sPow;
		sPow *= s;
		arVal[i++] = (T)sPow;
		sPow *= s;
		arVal[i++] = (T)sPow;
	}
}

template<class T>
void TMatrix<T>::Power2Matrix(const TMatrix<T>& vtS)	// for quadratics
{
	ASSERT(vtS.IsVector());
	SetSize(vtS.Length(), 3);		// make columns of vtS to power 0, 1, 2
	int i = 0;
	for (int row = 0; row < h; row++)
	{
		double s = vtS[row];
		arVal[i++] = 1;
		arVal[i++] = (T)s;
		arVal[i++] = (T)(s*s);
	}
}

template<class T>
void TMatrix<T>::Power1Matrix(const TMatrix<T>& vtS)	// for linear
{
	ASSERT(vtS.IsVector());
	SetSize(vtS.Length(), 2);		// make columns of vtS to power 0, 1
	int i = 0;
	for (int row = 0; row < h; row++)
	{
		arVal[i++] = 1;
		arVal[i++] = (T)vtS[row];
	}
}

template<class T>
void TMatrix<T>::SolvePolyAt(const TMatrix<T>& vectS, TMatrix<T>& vectRes) const
{
//		vectRes = vectS * this		assumes vectS is a row but returns a column vector
//		vectS is [1  s  s^2  s^3 ...] and can have extra elements!

	ASSERT(vectS.IsVector());
	ASSERT(vectS.Length() >= h);
	vectRes.SetSize(w, 1);			// result is column vector
	int inner = h;
	for (int col = 0; col < w; col++)
	{
		T elem = 0;
		int mxIdx = col;
		for (int i = 0; i < inner; i++)
		{
			elem += vectS.arVal[i] * arVal[mxIdx];
			mxIdx += w;
		}
		vectRes[col] = elem;
	}
}

template<class T>
void TMatrix<T>::SolvePolyAt(const TMatrix<T>& vectS, TVector<T>& vectRes) const
{
//		vectRes = vectS * this		assumes vectS is a row but returns a column vector
//		vectS is [1  s  s^2  s^3 ...] and can have extra elements!

	ASSERT(vectS.IsVector());
	ASSERT(vectS.Length() >= h);
	ASSERT(vectRes.Length() >= w);	// result is vector
	int inner = h;
	for (int col = 0; col < w; col++)
	{
		T elem = 0;
		int mxIdx = col;
		for (int i = 0; i < inner; i++)
		{
			elem += vectS.arVal[i] * arVal[mxIdx];
			mxIdx += w;
		}
		vectRes[col] = elem;
	}
}

template<class T>
void TMatrix<T>::SolvePolyAt0(TVector<T>& vectRes) const
{
//		vectRes = vectS * this		assumes vectS is a row but returns a column vector
//		vectS is [1  0  0  0 ...] and can have extra elements!
	ASSERT(vectRes.Length() >= w);	// result is vector
	for (int col = 0; col < w; col++)
		vectRes[col] = arVal[col];
}

template<class T>
void TMatrix<T>::SolvePolyAt1(TVector<T>& vectRes) const
{
//		vectRes = vectS * this		assumes vectS is a row but returns a column vector
//		vectS is [1  1  1  1 ...] and can have extra elements!
	ASSERT(vectRes.Length() >= w);	// result is vector
	int inner = h;
	for (int col = 0; col < w; col++)
	{
		T elem = 0;
		int mxIdx = col;
		for (int i = 0; i < inner; i++)
		{
			elem += arVal[mxIdx];
			mxIdx += w;
		}
		vectRes[col] = elem;
	}
}

template<class T>
void TMatrix<T>::Bezier2Poly(TMatrix<T>& poly) const
{
/*
	pre * by
	|  1  0  0  0 |
	| -3  3  0  0 |
	|  3 -6  3  0 |
	| -1  3 -3  1 |
*/
	ASSERT(h == 4);		// must be a cubic polynomial for now!
	poly.SetSize(h, w);
	int r0, r1, r2, r3;
	r0 = 0; r1 = w; r2 = 2*w; r3 = 3*w;
	while (r0 < w)
	{
		poly.arVal[r0] = arVal[r0];
		poly.arVal[r1] = 3 * (arVal[r1] - arVal[r0]);
		poly.arVal[r2] = 3 * (arVal[r2] - (2*arVal[r1]) + arVal[r0]);
		poly.arVal[r3] = arVal[r3] - arVal[r0] + (3 * (arVal[r1] - arVal[r2]));
		r0++; r1++; r2++; r3++;
	}
}

template<class T>
void TMatrix<T>::BezierRel2Poly(TMatrix<T>& poly) const
{
/*
	pre * by
	|  1  0  0  0 |
	|  0  3  0  0 |
	| -3 -6  3  3 |
	|  2  3 -3 -2 |
*/
	ASSERT(h == 4);		// must be a cubic polynomial for now!
	poly.SetSize(h, w);
	int r0, r1, r2, r3;
	r0 = 0; r1 = w; r2 = 2*w; r3 = 3*w;
	while (r0 < w)
	{
		poly.arVal[r0] = arVal[r0];
		poly.arVal[r1] = 3*arVal[r1];
		poly.arVal[r2] = 3 * (arVal[r3] + arVal[r2] - (2*arVal[r1]) - arVal[r0]);
		poly.arVal[r3] = (2 * (arVal[r0] - arVal[r3])) + (3 * (arVal[r1] - arVal[r2]));
		r0++; r1++; r2++; r3++;
	}
}

template<class T>
void TMatrix<T>::Poly2Bezier(TMatrix<T>& bezier) const
{
/*
	pre * by
	|  1   0   0   0 |
	|  1 1/3   0   0 |
	|  1 2/3 1/3   0 |
	|  1   1   1   1 |
*/
	ASSERT(h == 4);		// must be a cubic polynomial for now!
	bezier.SetSize(h, w);
	int r0, r1, r2, r3;
	r0 = 0; r1 = w; r2 = 2*w; r3 = 3*w;
	while (r0 < w)
	{
		bezier.arVal[r0] = arVal[r0];
		bezier.arVal[r1] = arVal[r0] + (arVal[r1]/3);
		bezier.arVal[r2] = arVal[r0] + (((2*arVal[r1]) + arVal[r2]) / 3);
		bezier.arVal[r3] = arVal[r0] + arVal[r1] + arVal[r2] + arVal[r3];
		r0++; r1++; r2++; r3++;
	}
}

template<class T>
void TMatrix<T>::Poly2BezierRel(TMatrix<T>& bezierRel) const
{
/*
	pre * by
	|  1    0    0    0 |
	|  0  1/3    0    0 |
	|  0 -1/3 -2/3   -1 |
	|  1    1    1    1 |
*/
	ASSERT(h == 4);		// must be a cubic polynomial for now!
	bezierRel.SetSize(h, w);
	int r0, r1, r2, r3;
	r0 = 0; r1 = w; r2 = 2*w; r3 = 3*w;
	while (r0 < w)
	{
		bezierRel.arVal[r0] = arVal[r0];
		bezierRel.arVal[r1] = arVal[r1]/3;
		bezierRel.arVal[r2] = ((-arVal[r1] - (2*arVal[r2])) / 3) - arVal[r3];
		bezierRel.arVal[r3] = arVal[r0] + arVal[r1] + arVal[r2] + arVal[r3];
		r0++; r1++; r2++; r3++;
	}
}

template<class T>
void TMatrix<T>::SetRow(int row, const TMatrix<T>& vect)
{
	ASSERT(row < h);
	ASSERT(vect.IsVector());
	if (vect.h != 1) TRACE0("TMatrix::SetRow() vector is not row!\n");
	if (vect.Length() != w) TRACE0("TMatrix::SetRow() row lengths don't match!\n");
	int maxI = __min(w, vect.Length());
	int mxIdx = row * w;
	for (int i = 0; i < maxI; i++)
		arVal[mxIdx++] = vect[i];
}
template<class T>
void TMatrix<T>::SetColumn(int col, const TMatrix<T>& vect)
{
	ASSERT(col < w);
	ASSERT(vect.IsVector());
	if (vect.w != 1) TRACE0("TMatrix::SetColumn() vector is not column!\n");
	if (vect.Length() != h) TRACE0("TMatrix::SetColumn() column lengths don't match!\n");
	int maxI = __min(h, vect.Length());
	int mxIdx = col;
	for (int i = 0; i < maxI; i++)
	{
		arVal[mxIdx] = vect[i];
		mxIdx += w;
	}
}

template<class T>
void TMatrix<T>::SetRow(int row, const TVector<T>& vect)
{
	ASSERT(row < h);
	ASSERT(vect.Length() == w);
	int maxI = w;
	int mxIdx = row * w;
	for (int i = 0; i < maxI; i++)
		arVal[mxIdx++] = vect[i];
}
template<class T>
void TMatrix<T>::SetColumn(int col, const TVector<T>& vect)
{
	ASSERT(col < w);
	ASSERT(vect.Length() == h);
	int maxI = h;
	int mxIdx = col;
	for (int i = 0; i < maxI; i++)
	{
		arVal[mxIdx] = vect[i];
		mxIdx += w;
	}
}

template<class T>
void TMatrix<T>::SetRow(int row, const T& val)
{
	ASSERT(row < h);
	int maxI = w;
	int mxIdx = row * w;
	for (int i = 0; i < maxI; i++)
		arVal[mxIdx++] = val;
}
template<class T>
void TMatrix<T>::SetColumn(int col, const T& val)
{
	ASSERT(col < w);
	int maxI = h;
	int mxIdx = col;
	for (int i = 0; i < maxI; i++)
	{
		arVal[mxIdx] = val;
		mxIdx += w;
	}
}

template<class T>
void TMatrix<T>::SetRow(int row, const T* arSrc)
{
	ASSERT(row < h);
	int maxI = w;
	int mxIdx = row * w;
	for (int i = 0; i < maxI; i++)
		arVal[mxIdx++] = arSrc[i];
}
template<class T>
void TMatrix<T>::SetColumn(int col, const T* arSrc)
{
	ASSERT(col < w);
	int maxI = h;
	int mxIdx = col;
	for (int i = 0; i < maxI; i++)
	{
		arVal[mxIdx] = arSrc[i];
		mxIdx += w;
	}
}

template<class T>
void TMatrix<T>::ScaleRowsBy(const TMatrix<T>& vt)
{
	const int rows = vt.Length();
	ASSERT(h == rows && vt.IsVector());
	for (int r = 0; r < rows; r++)
		for (int c = 0; c < w; c++)
			elem(r, c) *= vt[r];
}

template<class T>
void TMatrix<T>::ScaleRowsBy(const TVector2<T>& vt)
{
	const int rows = vt.Length();
	ASSERT(h >= rows);				// allow to scale only part
	for (int r = 0; r < rows; r++)
		for (int c = 0; c < w; c++)
			elem(r, c) *= vt[r];
}
template<class T>
void TMatrix<T>::ScaleRowsBy(const TVector<T>& vt)
{
	const int rows = vt.Length();
	ASSERT(h == rows);
	for (int r = 0; r < rows; r++)
		for (int c = 0; c < w; c++)
			elem(r, c) *= vt[r];
}

template<class T>
void TMatrix<T>::ScaleColumnsBy(const TMatrix<T>& vt)
{
	const int cols = vt.Length();
	ASSERT(w == cols && vt.IsVector());
	for (int c = 0; c < cols; c++)
		for (int r = 0; r < h; r++)
			elem(r, c) *= vt[c];
}

template<class T>
void TMatrix<T>::ScaleColumnsBy(const TVector2<T>& vt)
{
	const int cols = vt.Length();
	ASSERT(w >= cols);				// allow to scale only part
	for (int c = 0; c < cols; c++)
		for (int r = 0; r < h; r++)
			elem(r, c) *= vt[c];
}
template<class T>
void TMatrix<T>::ScaleColumnsBy(const TVector<T>& vt)
{
	const int cols = vt.Length();
	ASSERT(w == cols);
	for (int c = 0; c < cols; c++)
		for (int r = 0; r < h; r++)
			elem(r, c) *= vt[c];
}

template<class T>
void TMatrix<T>::SwapRows(int row1, int row2)
{
	ASSERT(row1 < h && row2 < h);
	row1 *= w;
	row2 *= w;
	for (int c = 0; c < w; c++)
	{
		T tmp = arVal[row1];
		arVal[row1] = arVal[row2];
		arVal[row2] = tmp;
		row1++;
		row2++;
	}
}
template<class T>
void TMatrix<T>::SwapColumns(int col1, int col2)
{
	ASSERT(col1 < w && col2 < w);
	for (int r = 0; r < h; r++)
	{
		T tmp = arVal[col1];
		arVal[col1] = arVal[col2];
		arVal[col2] = tmp;
		col1 += w;
		col2 += w;
	}
}

template<class T>
void TMatrix<T>::CopyRow(int rowDest, int rowSrc, int num)
{
	ASSERT(rowDest <= h-num && rowSrc <= h-num);
	rowDest *= w;
	rowSrc *= w;
	int numElems = num * w;
	for (int c = 0; c < numElems; c++)
		arVal[rowDest++] = arVal[rowSrc++];
}
template<class T>
void TMatrix<T>::CopyColumn(int colDest, int colSrc, int num)
{
	ASSERT(colDest <= w-num && colSrc <= w-num);
	for (int c = 0; c < num; c++)
	{
		int iDest = colDest + c;
		int iSrc = colSrc + c;
		for (int r = 0; r < h; r++)
		{
			arVal[iDest] = arVal[iSrc];
			iDest += w;
			iSrc += w;
		}
	}
}


template<class T>
void TMatrix<T>::BlockScale(const TMatrix<T>& mxBlock, const TMatrix<T>& mxScaler)
{
	SetSize(mxBlock.h*mxScaler.h, mxBlock.w*mxScaler.w);
	int rowsS = mxScaler.h;
	int colsS = mxScaler.w;
	int rowsB = mxBlock.h;
	int colsB = mxBlock.w;
	for (int rS = 0; rS < rowsS; rS++)
		for (int cS = 0; cS < colsS; cS++)
		{
			double valS = mxScaler.elem(rS, cS);
			int idxB = 0;
			for (int rB = 0; rB < rowsB; rB++)
			{
				int idxMx = (rS*rowsB + rB) * w + cS*colsB;
				for (int cB = 0; cB < colsB; cB++)
				//	elem(rS*rowsB + rB, cS*colsB + cB) = mxBlock.elem(rB, cB) * valS;
					arVal[idxMx++] = mxBlock[idxB++] * valS;
			}
		}
}





template<class T>
void TMatrix<T>::Derivative(TMatrix<T>& deriv) const
{
	if (h == 1)
	{
		deriv.SetSize(1, w);
		deriv = 0;
		return;
	}
	deriv.SetSize(h-1, w);
	int derIdx = 0;
	int mxIdx = w;
	for (int row = 1; row < h; row++)
		for (int col = 0; col < w; col++)
			deriv.arVal[derIdx++] = row * arVal[mxIdx++];
}

template<class T>
void TMatrix<T>::Intergral(TMatrix<T>& inter) const
{
	inter.SetSize(h+1, w);
	int intIdx = w;
	int mxIdx = 0;
	for (int row = 1; row <= h; row++)
		for (int col = 0; col < w; col++)
			inter.arVal[intIdx++] = arVal[mxIdx++] / row;
	for (col = 0; col < w; col++)
		inter.arVal[col] = 0;		// set top row (s^0 coef's) to 0
}



//////////////////////////////////////
// Matrix inverse related functions
//////////////////////////////////////







/////////////////////////////////////
// TMatrix associated functions
/////////////////////////////////////

template<class T>
TMatrix<T> abs(const TMatrix<T>& matrix)
{
	TMatrix<T> res(matrix.w, matrix.h);
	int elements = matrix.w * matrix.h;
	for (int i = 0; i < elements; i++)
		res.arVal[i] = fabs(matrix.arVal[i]);
//		res.arVal[i] = matrix.arVal[i] >= 0 ? matrix.arVal[i] : -matrix.arVal[i]);
	return res;
}

template<class T>
TMatrix<T> cross(const TMatrix<T>& a, const TMatrix<T>& b)
{
	ASSERT(a.w * a.h == 3 && b.w * b.h == 3);
	TMatrix<T> res(1, 3);		// (rows, cols) format
	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = a[2]*b[0] - a[0]*b[2];
	res[2] = a[0]*b[1] - a[1]*b[0];
	return res;
}

template<class T>
T dot(const TMatrix<T>& a, const TMatrix<T>& b)
{
	ASSERT(a.w * a.h == b.w * b.h);
	T res = 0;
	int elements = a.w * a.h;
	for (int i = 0; i < elements; i++)
		res += a[i] * b[i];
	return res;
}

template<class T>
TMatrix<T> MultiplyElems(const TMatrix<T>& a, const TMatrix<T>& b)
{
	ASSERT(a.w == b.w && a.h == b.h);
	TMatrix<T> res(a.w, a.h);
	int elements = a.w * a.h;
	for (int i = 0; i < elements; i++)
		res[i] = a[i] * b[i];
	return res;
}






/////////////////////////////////////
// TMatrix associated operators
/////////////////////////////////////

template<class T>
ostream& operator<<(ostream& os, const TMatrix<T>& matrix)
{
//	os << "Rows = " << matrix.h << ", Columns = " << matrix.w << endl;
//	os.precision(3);
	int mxIdx = 0;
	for (int row = 0; row < matrix.h; row++)
	{
		for (int col = 0; col < matrix.w; col++)
		{
			os.width(8);
			os << matrix[mxIdx++];
		}
		if (row != 0 || matrix.h != 1)		// if only one row don't do endl
			os << endl;
	}
	return os;
}

#if defined _AFX			// not too sure about _AFX's definition!
template<class T>
CDumpContext& operator<<(CDumpContext& dc, const TMatrix<T>& matrix)
{
	ostrstream os;
	os << "Rows = " << matrix.h << ", Columns = " << matrix.w << endl;
	os.setf(ios::showpoint | ios::left);
	os.precision(4);
	int mxIdx = 0;
	for (int row = 0; row < matrix.h; row++)
	{
		os << "   ";
		for (int col = 0; col < matrix.w; col++)
		{
			os.width(13);
			os << matrix[mxIdx++];
		}
		os << endl;
	}
	os << ends;		// adds '\0'
	dc << os.str();
	os.rdbuf()->freeze(0);
	return dc;
}
#endif


template<class T>
TMatrix<BOOL_TYPE> operator==(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<BOOL_TYPE> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] == rhs[i];
	return res;
}

template<class T>
TMatrix<BOOL_TYPE> operator==(const TMatrix<T>& lhs, const T& rhs)
{
	TMatrix<BOOL_TYPE> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] == rhs;
	return res;
}

template<class T>
TMatrix<BOOL_TYPE> operator>(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<BOOL_TYPE> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] > rhs[i];
	return res;
}

template<class T>
TMatrix<BOOL_TYPE> operator<(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<BOOL_TYPE> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] < rhs[i];
	return res;
}

template<class T>
TMatrix<T> operator*(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.h);
	TMatrix<T> res(lhs.h, rhs.w);
	int inner = lhs.w;
	T elem;
	for (int row = 0; row < lhs.h; row++)
		for (int col = 0; col < rhs.w; col++)
		{
			elem = 0;
			int mxIdx = row*inner;
			for (int i = 0; i < inner; i++)
				elem += lhs[mxIdx++] * rhs[i*rhs.w + col];
			res[row*res.w + col] = elem;
		}
	return res;
}

template<class T>
TMatrix<T> operator*(const TMatrix<T>& lhs, const TVector<T>& rhs)
{
	ASSERT(lhs.w == 3);
	TMatrix<T> res(lhs.h, 1);
	int inner = lhs.w;
	int mxIdx = 0;
	T elem;
	for (int row = 0; row < lhs.h; row++)
		{
			elem = 0;
			for (int i = 0; i < inner; i++)
				elem += lhs[mxIdx++] * rhs[i];
			res[row] = elem;
		}
	return res;
}

template<class T>
TMatrix<T> operator*(const TMatrix<T>& lhs, double rhs)
{
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = (T)(lhs[i] * rhs);
	return res;
}

template<class T>
TMatrix<T> operator*(const T& lhs, const TMatrix<T>& rhs)
{
	TMatrix<T> res(rhs.h, rhs.w);
	int elements = rhs.h * rhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs * rhs[i];
	return res;
}

template<class T>
TMatrix<T> operator/(const TMatrix<T>& lhs, double rhs)
{
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = (T)(lhs[i] / rhs);
	return res;
}

template<class T>
TMatrix<T> operator+(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] + rhs[i];
	return res;
}

template<class T>
TMatrix<T> operator+(const TMatrix<T>& lhs, const T& rhs)
{
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] + rhs;
	return res;
}

template<class T>
TMatrix<T> operator-(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] - rhs[i];
	return res;
}

template<class T>
TMatrix<T> operator-(const TMatrix<T>& lhs, const T& rhs)
{
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] - rhs;
	return res;
}

template<class T>
TMatrix<T> operator-(const TMatrix<T>& matrix)
{
	TMatrix<T> res(matrix.h, matrix.w);
	int elements = matrix.h * matrix.w;
	for (int i = 0; i < elements; i++)
		res[i] = -matrix[i];
	return res;
}

template<class T>
TMatrix<T> operator||(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] || rhs[i];
	return res;
}

template<class T>
TMatrix<T> operator&&(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
	ASSERT(lhs.w == rhs.w && lhs.h == rhs.h);
	TMatrix<T> res(lhs.h, lhs.w);
	int elements = lhs.h * lhs.w;
	for (int i = 0; i < elements; i++)
		res[i] = lhs[i] && rhs[i];
	return res;
}








#endif // !defined(AFX_MATRIX_H__48FCACC3_FE73_11D4_8C1E_9D2D1A323C2C__INCLUDED_)
