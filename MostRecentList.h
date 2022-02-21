//////////////////////////
// MostRecentList.h
//////////////////////////

#if !defined(MOSTRECENTLIST_H__INCLUDED_)
#define MOSTRECENTLIST_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000



/*
	TMostRecentList

	Records the most recently added items for referencing
	Oldest are overwritten once list is full


*/

template<class T> class TMostRecentList
{
protected:
	T* m_pData;
	int m_sizeBuffer;
	bool m_bBufferNewed;
	int m_idxLastIn;		// last In
	int m_iCount;

public:
	TMostRecentList();
	virtual ~TMostRecentList();
	bool SetBufferSize(int size);
	int GetBufferSize() { return m_sizeBuffer; }
	void UseBuffer(T* pBuff, int size);
	void RemoveAll() { m_idxLastIn = m_iCount = 0; }
	bool IsEmpty() { return m_iCount == 0; }
	int GetCount() { return m_iCount; }
	void Add(T& elem);
	T& operator[](int idx);

};


template<class T>
TMostRecentList<T>::TMostRecentList()
{
	m_pData = NULL;
	m_bBufferNewed = false;
	m_sizeBuffer = 0;
	RemoveAll();
}

template<class T>
TMostRecentList<T>::~TMostRecentList()
{
	if (m_bBufferNewed && m_pData)
		delete[] m_pData;
}

template<class T>
bool TMostRecentList<T>::SetBufferSize(int size)
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
void TMostRecentList<T>::UseBuffer(T* pBuff, int size)
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
void TMostRecentList<T>::Add(T& elem)
{
	if (--m_idxLastIn < 0)
		m_idxLastIn += m_sizeBuffer;
	ASSERT(m_idxLastIn >= 0 && m_idxLastIn < m_sizeBuffer);
	m_pData[m_idxLastIn] = elem;
	if (++m_iCount > m_sizeBuffer)
		m_iCount = m_sizeBuffer;
}

template<class T>
T& TMostRecentList<T>::operator[](int idx)	// 0 is last in, 1 prev in etc.
{
	idx += m_idxLastIn;
	if (idx >= m_sizeBuffer)
		idx -= m_sizeBuffer;
	ASSERT(idx >= 0 && idx < m_sizeBuffer);
	return m_pData[idx];
}





#endif	// !defined(MOSTRECENTLIST_H__INCLUDED_)
