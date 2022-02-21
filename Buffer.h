//////////////////////////
// Buffer.h
//////////////////////////

#if !defined(BUFFER_H__INCLUDED_)
#define BUFFER_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000



///////////////////////////
// Buffers


template<class T> class TFIFO		// in and out can be seperate threads
{
protected:
	T* m_pData;
	int m_sizeBuffer;
	bool m_bBufferNewed;
	int m_idxIn;		// next In
	int m_idxOut;		// next Out - if in == out then empty, must be at least one space free!

public:
	TFIFO();
	virtual ~TFIFO();
	bool SetBufferSize(int size);
	void UseBuffer(T* pBuff, int size);
	int GetBufferSize() const { return m_sizeBuffer; }
	int GetMaxSize() const { return m_sizeBuffer - 1; }	// one extra space
	bool IsEmpty() const { return (m_idxIn == m_idxOut); }
	int GetCount() const;
	int GetFree() const;
	int Add(T& elem);					// returns number of free elements after adding, -1 if already full
	int Add(const T* pElem, int num);	// returns number of free elements after adding, -'ve indicates number no fitted
	T* GetNextToAdd() const { return &m_pData[m_idxIn]; }
	void AddOne() { if (++m_idxIn >= m_sizeBuffer) m_idxIn = 0; }	// adds object in after being set with GetNextToAdd()
	T& Remove();					// ASSERTs if empty
	bool Remove(T& elem);		// returns true if object valid
	void RemoveAll() { m_idxIn = m_idxOut = 0; }
	void Discard() { if (!IsEmpty() && (++m_idxOut >= m_sizeBuffer)) m_idxOut = 0; }	// discard next object
	void Discard(int size);		// discards 'size' objects
	T& GetNext() const { ASSERT(!IsEmpty()); return m_pData[m_idxOut]; }
	bool GetNextBlock(T*& pElem, int& size) const;

	T& GetIter(int iIter) const { ASSERT(!IsEmpty()); return m_pData[iIter]; }
	T* GetNextIter(int& iIter);
	T* GetPrevIter(int& iIter);
	bool CheckIter(int iIter);
	int GetNextOutIter() { ASSERT(!IsEmpty()); return m_idxOut; }
	int GetLastInIter() { ASSERT(!IsEmpty()); return m_idxIn != 0 ? m_idxIn-1 : m_sizeBuffer-1; }

};

template<class T>
TFIFO<T>::TFIFO()
{
	m_pData = NULL;
	m_bBufferNewed = false;
	m_idxIn = m_idxOut = m_sizeBuffer = 0;
}

template<class T>
TFIFO<T>::~TFIFO()
{
	if (m_bBufferNewed && m_pData)
		delete[] m_pData;
}

template<class T>
bool TFIFO<T>::SetBufferSize(int size)
{
	ASSERT(IsEmpty());
	RemoveAll();
	if (m_bBufferNewed && m_pData)
		delete[] m_pData;
	size++;				// allow for one free space
	m_pData = new T[size];
	m_sizeBuffer = size;
	m_bBufferNewed = true;
	return m_pData != NULL;
}

template<class T>
void TFIFO<T>::UseBuffer(T* pBuff, int size)
{
	ASSERT(IsEmpty());
	RemoveAll();
	if (m_bBufferNewed && m_pData)
		delete[] m_pData;
	m_bBufferNewed = false;
	m_pData = pBuff;
	m_sizeBuffer = size;
}

template<class T>
int TFIFO<T>::GetCount() const
{
	int sz = m_idxIn - m_idxOut;
	if (sz < 0)
		sz += m_sizeBuffer;
	return sz;
}

template<class T>
int TFIFO<T>::GetFree() const
{
	int sz = m_idxOut - m_idxIn - 1;		// reserves one space
	if (sz < 0)
		sz += m_sizeBuffer;
	return sz;
}


template<class T>
int TFIFO<T>::Add(T& elem)		// returns number of free elements after adding, -1 if already full
{
	int iFreeAfter = GetFree() - 1;
	if (iFreeAfter >= 0)
	{
		m_pData[m_idxIn] = elem;	// elem assigned before idx increased - safe for Multi Thread
		if (++m_idxIn >= m_sizeBuffer)
			m_idxIn  = 0;
	}
	return iFreeAfter;
}

template<class T>
int TFIFO<T>::Add(const T* pElem, int num)		// returns number of free elements after adding, -1 if already full
{
	int iFreeAfter = GetFree() - num;
	if (iFreeAfter < 0)
		return iFreeAfter;
//		num += iFreeAfter;		// set num to GetFree()
	int idxElem = 0;
	while (num-- > 0)
	{
		m_pData[m_idxIn] = pElem[idxElem++];	// elem assigned before idx increased - safe for Multi Thread
		if (++m_idxIn >= m_sizeBuffer)
			m_idxIn  = 0;
	}
	return iFreeAfter;		// if negative indicates number that didn't get added
}

template<class T>
bool TFIFO<T>::GetNextBlock(T*& pElem, int& size) const
{
	pElem = &m_pData[m_idxOut];
	if (m_idxIn >= m_idxOut)
	{
		size = m_idxIn - m_idxOut;
		return false;			// no more to get
	}
	size = m_sizeBuffer - m_idxOut;
	return m_idxIn != 0;		// there will be more if m_idxIn > 0
}

template<class T>
T& TFIFO<T>::Remove()
{
	ASSERT(!IsEmpty());
	int idxO = m_idxOut;
	if (++m_idxOut >= m_sizeBuffer)
		m_idxOut = 0;
	return m_pData[idxO];
}	// this elem won't be overwritten till next Get(), as one buffer space is always left

template<class T>
bool TFIFO<T>::Remove(T& elem)		// returns true if object valid
{
	if (IsEmpty())
		return false;
	elem = m_pData[m_idxOut];
	if (++m_idxOut >= m_sizeBuffer)
		m_idxOut = 0;
	return true;
}	// this elem won't be overwritten till next Get(), as one buffer space is always left

template<class T>
void TFIFO<T>::Discard(int size)
{
	if (size < GetCount())
	{
		int idxO = m_idxOut + size;	// so m_idxOut isn't temporaly invalid for multi threads
		if (idxO >= m_sizeBuffer)
			idxO -= m_sizeBuffer;
		m_idxOut = idxO;
	}
	else
		m_idxOut = m_idxIn;		// GetCount will be 0
}


/*
template<class T>
T* TFIFO<T>::AddEmpty()		// NOT SAFE for Multi Thread unless used after object copied in with GetNextToAdd()
{
	int idxIn = m_idxIn;
	if (++m_idxIn >= m_sizeBuffer)
		m_idxIn = 0;
	if (m_idxIn == m_idxOut)		// was full!
	{
		m_idxIn = idxIn;
		return NULL;
	}
	return &m_pData[idxIn];
}

template<class T>
void TFIFO<T>::RemoveLastAdded()		// NOT SAFE for Multi Thread
{
	if (--m_idxIn < 0)
		m_idxIn += m_sizeBuffer;
}
*/




template<class T>
T* TFIFO<T>::GetNextIter(int& iIter)
{
	ASSERT(!IsEmpty());
	if (iIter == m_idxIn)		// at head already
		return NULL;
	if (++iIter >= m_sizeBuffer)
		iIter = 0;
	if (iIter == m_idxIn)		// got to head
		return NULL;
	return &m_pData[iIter];
}

template<class T>
T* TFIFO<T>::GetPrevIter(int& iIter)
{
	ASSERT(!IsEmpty());
	if (iIter == m_idxOut)		// at tail already
		return NULL;
	if (--iIter < 0)
		iIter += m_sizeBuffer;
	return m_pData[iIter];
}

template<class T>
bool TFIFO<T>::CheckIter(int iIter)
{
	if (m_idxIn >= m_idxOut)		// not wrapped around
		return iIter <= m_idxIn && iIter >= m_idxOut;
	else
		return iIter <= m_idxIn || iIter >= m_idxOut;
}


#endif	// !defined(BUFFER_H__INCLUDED_)
