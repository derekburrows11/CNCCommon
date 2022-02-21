#if !defined STRUTILS_H
#define STRUTILS_H


#define skipWS(str)     while (isspace(*str)) str++
#define extractWS(str)  while (isspace(*str)) str++

#define extractCh(str, ch)  (*str == (ch) ? str++ : 0)


void RemoveWS(char* str);
// Removes whitespace from start and end of string
void RemoveWSStart(char* str);
// Removes whitespace from start of string
void RemoveWSEnd(char* str);
// Removes whitespace from end of string

// compares string to end of substring, returns length of strSub, or 0 if not equal
int strcmpsub(const char* string, const char* strSub);
int strextract(char*& str, const char* strSub);
// case insensitive
int stricmpsub(const char* string, const char* strSub);
int striextract(char*& str, const char* strSub);


////////////////////

class CFileBuffer
{
public:
	CFileBuffer(CFile* pFile);
	char* GetLine();
	int GetLineCount() { return m_iLineCount; }
	int GetLineLength() { return m_iLineLength; }

private:
	CFile* m_pFile;
	char  m_pBuffer[1024];
	char* m_pBufferEnd;
	char* m_pNextLine;
	int m_sizeBuffer;
	int m_iLineCount;
	int m_iLineLength;
	bool m_bBufferEndedWithR;
};

/////////////////////////

class CInterpStr {
public:
	enum {
		FMT_BIN = 1,
		FMT_OCT,
		FMT_DEC,
		FMT_HEX,				// sequential
	};		

 protected:
	char* m_Str;
	bool m_bStringNewed;
	int m_Idx, m_ItemStart, m_ItemLen, m_BuffLen;
	char m_ItemEndChar;
	int m_FormatSpec, m_FormatType, m_FormatDef;

 public:
	CInterpStr();
	CInterpStr(char* string);
	CInterpStr(int strSize);
	~CInterpStr();

	void SetToString(char* str);
	void SetStrSize(int strSize);
	int EatWS();
	bool Extract(char ch);
	int ExtractString(const char* str);
	int ExtractIString(const char* str);
	int IsItemNum();
	int IsItemAlphaNum();
	int IsItemAlpha();
	int GetAlphaNumItem(char* dest, int bufflen);	// returns num of chars
	int GetAlphaItem(char* dest, int bufflen);
	int GetNumItem(char* dest, int bufflen);
	int GetNumValue();							// returns value
	bool GetInt(int& num);

	bool AtEOL() { return m_Str[m_Idx] == 0; }
	bool AtEOnonWS();
	char* Buff() { return m_Str; }
	int BuffLen() { return m_BuffLen; }
	void Reset() { m_Idx = 0; }
	void Flush() { Reset(); m_Str[m_Idx] = 0; }
	void SetIdx(int idx) { m_Idx = idx; }
	int GetIdx() { return m_Idx; }
	char* GetIdxPtr() { return &m_Str[m_Idx]; }

	int operator !() const;
	operator char() const { return m_Str[m_Idx]; }
	operator char*() const { return &m_Str[m_Idx]; }
	void operator +=(int inc) { m_Idx += inc; }

};


#endif	// !defined STRUTILS_H