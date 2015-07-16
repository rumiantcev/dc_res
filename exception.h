// ---------------------------------------------------------------------------

#ifndef exceptionH
#define exceptionH

// ---------------------------------------------------------------------------
class Exception {
public:
	Exception(const AnsiString s);
   ///	Exception(int i);
	Exception(const Exception& e);
	~Exception();
	AnsiString msg() const ;

private:
	AnsiString what;
};

#endif
