
#include "stdafx.h"

#include "StrUtils.h"
//#include <ctype.h>
//#include <string.h>


// Removes whitespace from start and end of string
void RemoveWS(char* str)
{
	int i = 0;
	while (isspace(str[i])) i++;
	int si = i, ei = 0;
	while (str[i])
		if (!isspace(str[i++])) ei = i;
	str[ei] = 0;
	if (si)
		strcpy(str, str+si);
}

// Removes whitespace from start of string
void RemoveWSStart(char* str)
{
	int i = 0;
	while (isspace(str[i])) i++;
	if (i) strcpy(str, str+i);
}

// Removes whitespace from end of string
void RemoveWSEnd(char* str)
{
	int i = 0, ei = 0;
	while (str[i])
		if (!isspace(str[i++])) ei = i;
	str[ei] = 0;
}


// compares string to end of substring, returns length of strSub, or 0 if not equal
int strcmpsub(const char* str, const char* strSub)
{
	for (int i = 0; strSub[i] != 0; i++)
		if (str[i] != strSub[i])
			return 0;
	return i;
}

// case insensitive
int stricmpsub(const char* str, const char* strSub)
{
	for (int i = 0; strSub[i] != 0; i++)
		if (str[i] != strSub[i] && tolower(str[i]) != tolower(strSub[i]))
			return 0;
	return i;
}

int strextract(char*& str, const char* strSub)
{
	for (int i = 0; strSub[i] != 0; i++)
		if (str[i] != strSub[i])
			return 0;
	str += i;
	return i;
}

// case insensitive
int striextract(char*& str, const char* strSub)
{
	for (int i = 0; strSub[i] != 0; i++)
		if (str[i] != strSub[i] && tolower(str[i]) != tolower(strSub[i]))
			return 0;
	str += i;
	return i;
}


////////////////////////////////////////////
// class CFileBuffer

CFileBuffer::CFileBuffer(CFile* pFile)
{
	m_pFile = pFile;
	m_sizeBuffer = sizeof(m_pBuffer) - 1;		// leave spare char for '\0'
	m_pBufferEnd = m_pBuffer + m_sizeBuffer;
	m_pNextLine  = m_pBufferEnd;
	m_iLineCount = 0;
	m_bBufferEndedWithR = false;
}

char* CFileBuffer::GetLine()
{
	char* pLine = m_pNextLine;
	char* pCurr = m_pNextLine;
	char* pBufferEnd = m_pBufferEnd;
	for (;;)
	{
		if (pCurr >= pBufferEnd)
		{
			int offset = pCurr - pLine;
			if (offset == 0 && (m_pBufferEnd != m_pBuffer + m_sizeBuffer))
			{
				m_iLineLength = 0;
				return NULL;				// EOF and no line to return
			}
			// Reload Buffer
			// if a full line is not in buffer move it to start and refill buffer
			if (pLine == m_pBuffer)		// line is longer than buffer
			{
				*pCurr = '\0';			// pCurr is within buffer - buffer is 1 longer than m_sizeBuffer
				m_pNextLine = pCurr;	// keep m_pNextLine == pBufferEnd for next pass
				m_iLineLength = pCurr - pLine;
				break;					// just return buffer full, line count doesn't increment
			}
			int numToMove = m_pBufferEnd - pLine;
			if (numToMove != 0)
				memmove(m_pBuffer, pLine, numToMove);
			pLine = m_pBuffer;
			int numToRead = m_sizeBuffer - numToMove;		// = pCurr - m_pBuffer
			int num = m_pFile->Read(m_pBuffer + numToMove, numToRead);
			ASSERT(num <= numToRead);
			ASSERT(numToRead != 0);
//			m_bGotEOF = (num != numRead);
			m_pBufferEnd = m_pBuffer + numToMove + num;
			pBufferEnd = m_pBufferEnd;
			pCurr = pLine + offset;
		}

		if (*pCurr == '\r')
		{
			*pCurr = '\0';
			m_pNextLine = pCurr + 1;
			m_iLineCount++;
			m_iLineLength = pCurr - pLine;
			pCurr++;
			if (pCurr < pBufferEnd)
			{
				m_bBufferEndedWithR = false;
				if (*pCurr == '\n')		// check if is "\r\n"
					m_pNextLine++;
			}
			else
			{
				ASSERT(pCurr == pBufferEnd);
				m_bBufferEndedWithR = true;
			}
			break;
		}
		else if (*pCurr == '\n')
		{
			if (m_bBufferEndedWithR && pCurr == pLine)		// ignore if "\r\n" pair split across buffer
			{
				m_bBufferEndedWithR = false;
			}
			else
			{
				*pCurr = '\0';
				m_pNextLine = pCurr + 1;
				m_iLineCount++;
				m_iLineLength = pCurr - pLine;
/*				pCurr--
				if (pCurr >= pLine && *pCurr == '\r')		// check if was "\r\n"
				{
					*pCurr = '\0';
					m_iLineLength--;
				}
*/
				m_bBufferEndedWithR = false;
				break;
			}
		}
		pCurr++;
	}
	return pLine;
}


//////////////////////////////////////////////////////
// class CInterpStr
// utility class for interperating stings

CInterpStr::CInterpStr()
{
	m_Str = NULL;
	m_bStringNewed = false;
	m_BuffLen = 0;
	Reset();
	m_FormatDef = FMT_DEC;
}

CInterpStr::CInterpStr(char* str)
{
	m_Str = str;
	m_bStringNewed = false;
	m_BuffLen = strlen(m_Str) + 1;
	Reset();
	m_FormatDef = FMT_DEC;
}

CInterpStr::CInterpStr(int strSize)
{
	m_bStringNewed = false;
	SetStrSize(strSize);
}

CInterpStr::~CInterpStr()
{
	if (m_bStringNewed)
		delete[] m_Str;
}

void CInterpStr::SetToString(char* str)
{
	if (m_bStringNewed)
		delete[] m_Str;
	m_Str = str;
	m_bStringNewed = false;
	m_BuffLen = strlen(m_Str) + 1;
	Reset();
}

void CInterpStr::SetStrSize(int strSize)
{
	if (m_bStringNewed == true)
		delete[] m_Str;
	m_Str = new char[strSize];
	m_bStringNewed = true;
	m_BuffLen = strSize;
	Flush();
	m_FormatDef = FMT_DEC;
}

bool CInterpStr::AtEOnonWS()		// true if at EOL or only WS remaining
{
	int i = m_Idx;
	while (isspace(m_Str[i])) i++;
	return m_Str[i] == 0;
}

int CInterpStr::EatWS()		// returns number of whitespace chars eaten
{
	int si = m_Idx;
	while (isspace(m_Str[m_Idx])) m_Idx++;	// seek to end of whitespace
	return m_Idx - si;
}

bool CInterpStr::Extract(char ch)
{
	if (m_Str[m_Idx] != ch)
		return false;
	m_Idx++;
	return true;
}

int CInterpStr::ExtractString(const char* str)
{
	// returns number of chars in string  if found, or 0
	int i = -1;
	while (str[++i] != 0)
		if (str[i] != m_Str[m_Idx + i])
			{	i = 0; break; }
	m_Idx += i;
	return i;
}

int CInterpStr::ExtractIString(const char* str)
{
	// case insensitive - returns number of chars in string  if found, or 0
	int i = -1;
	while (str[++i] != 0)
		if (tolower(str[i]) != tolower(m_Str[m_Idx + i]))
			{	i = 0; break; }
	m_Idx += i;
	return i;
}

int CInterpStr::IsItemNum()	// the sequence of alphanum are all hex or dec digits
{
	m_ItemStart = m_Idx;
	int i = 0;
	m_FormatSpec = 0;
	m_FormatType = FMT_BIN;
	bool bIsNum = false;
	bool bSigned = false;

	for (;;)
	{
		char ch = m_Str[m_Idx + i];
		if (isxdigit(ch))
		{
			bIsNum = true;
			if (!isdigit(ch))		// Not a decimal digit
				m_FormatType = max(m_FormatType, FMT_HEX);		// has hex digit (>= a)
			else if (ch > '7')
				m_FormatType = max(m_FormatType, FMT_DEC);		// has hex or dec digit (>= 8)
			else if (ch > '1')
				m_FormatType = max(m_FormatType, FMT_OCT);		// has hex, dec or oct digit (>= 2
		}
		else if ((i == 0) && (ch == '-' || ch == '+'))
			bSigned = true;			// Num OK
		else if (ch == 'x' || ch == 'X')
		{					// check previous char is '0' and first digit
			bIsNum = false;			// Will need following digits
			if ( ((i == 1 && !bSigned) || (i == 2 && bSigned))
				&& (m_Str[m_Idx + i - 1] == '0') )
			{
				m_FormatSpec = FMT_HEX;
				m_FormatType = FMT_HEX;
			}
			else
				break;		// Invalid location of 'x' bIsNum already set to false
		}
		else if (!isalnum(ch))		// end of item chars - ok
			break;					// bIsNum is valid flag
		else						// invalid ch - not number
		{
			bIsNum = false;
			break;
		}
		i++;
	}		// for (;;)
	m_ItemLen = bIsNum ? i : 0;
	return m_ItemLen;
}

int CInterpStr::GetNumValue()
{
	if (!m_ItemLen)
		return 0;
	int val = -1;
	char endChar = m_Str[m_Idx + m_ItemLen];
	m_Str[m_Idx + m_ItemLen] = 0;

	int fmt = 0;
	if (m_FormatSpec)
		fmt = m_FormatSpec;
	else if (m_FormatType == FMT_HEX)			// if hex then can only be hex
		fmt = m_FormatType;
	else
		fmt = m_FormatDef;

	if (fmt < m_FormatType)		// check specified format is valid
		fmt = 0;
	if (fmt == FMT_DEC)
		sscanf(&m_Str[m_Idx], "%d", &val);
	else if (fmt == FMT_HEX)
		sscanf(&m_Str[m_Idx], "%x", &val);
	
	m_Str[m_Idx + m_ItemLen] = endChar;
	m_Idx += m_ItemLen;			// advance index
	return val;
}

int CInterpStr::IsItemAlphaNum()	// a sequence of alphanum
{
	m_ItemStart = m_Idx;
	int i = m_Idx;
	while (isalnum(m_Str[i])) i++;
	m_ItemEndChar = m_Str[i];
	m_ItemLen = i - m_ItemStart;
	return m_ItemLen;
}

bool CInterpStr::GetInt(int& num)
{
	int i = m_Idx;
	while (isdigit(m_Str[i])) i++;
	if (i == m_Idx)
		return false;			// no digits
	char chEnd = m_Str[i];
	m_Str[i] = 0;
	sscanf(&m_Str[m_Idx], "%i", &num);
	m_Str[i] = chEnd;
	m_Idx = i;
	return true;
}

int CInterpStr::GetAlphaNumItem(char* dest, int bufflen)	// returns num of chars
{
	int i = 0;
	while (i < bufflen-1 && isalnum(dest[i] = m_Str[i+m_Idx]))
		i++;
	dest[i] = 0;
	m_ItemEndChar = m_Str[i+m_Idx];
	m_Idx += i;				// advance index
	m_ItemLen = i;
	return m_ItemLen;
}


