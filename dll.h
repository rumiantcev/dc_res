#ifndef _DLL_H_
#define _DLL_H_

#pragma once

// #include <iostream>
#include <string>
#include <exception>
#include <System.hpp>
#include <SysUtils.hpp>

// ===========================================================================
//
// DLL Defs
//
// ===========================================================================
/* #ifdef UNICODE
 #define LoadLibrary  LoadLibraryW
 #else
 #define LoadLibrary  LoadLibraryA
 #endif // !UNICODE
/* */

namespace Dll {

	__declspec(selectany) std::string ErrDllMsg =
		"Ошибка загрузки библиотеки:%s";
	__declspec(selectany) std::string ErrProcMsg =
		"Ошибка получения адреса функции:%s";

	/* class Exception {
	 private:
	 AnsiString message;

	 public:
	 __fastcall Exception(const AnsiString mg) : message(mg) {
	 };
	 __property AnsiString Msg = {read = message, write = message};
	 };
	 */

	class EDllError /* : public Exception */ {
	private:
		std::string m_Name;

	public:
		// __fastcall EDllError(const AnsiString a_Name) {
		// m_Name = a_Name;
		// };
		__fastcall EDllError(const std::string a_Fmt, const std::string a_Name)
			/* : Exception(Format(a_Fmt,ARRAYOFCONST((a_Name)))), m_Name(a_Name) */
		{
			Name = a_Name;
		};

		/* __property */
		std::string Name; // = {read = m_Name, write = m_Name};
	};

	class EDllProcError : public EDllError {
	public:
		// __fastcall EDllProcError(const AnsiString a_Name) : EDllError(a_Name) {
		// };
		__fastcall EDllProcError(const std::string a_Fmt,
			const std::string a_Name) : EDllError(a_Fmt, a_Name) {
		};
	}; /* */

	// ----------------------dll
	class TDll {
		HINSTANCE m_DLLInstance;

	public:
		TDll(const char* a_Name) {
# pragma option push  -w-pia
			if (!(m_DLLInstance = LoadLibraryA(a_Name)))
# pragma option pop
				throw(EDllError(ErrDllMsg, a_Name));
		}

		~TDll() {
			if (m_DLLInstance)
				FreeLibrary(m_DLLInstance);
		}

		operator HINSTANCE() const {
			return m_DLLInstance;
		}
	};

	// ------------ссылка на процедуру в dll
	class TDllProc {
	public:
		TDllProc(const TDll& a_Dll, const char* a_Name) {
# pragma option push -w-pia
			if (!(m_Proc = GetProcAddress(HINSTANCE(a_Dll), a_Name)))
# pragma option pop
				throw(EDllProcError(ErrProcMsg, a_Name));
		}

	public:
		FARPROC m_Proc;
	};

	// -----вызов процедуры без парметров, и результатов
	class TDllProcV0 : public TDllProc {
	public:
		TDllProcV0(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()() {
			typedef void(*TProc)();
			((TProc)m_Proc)();
		}
	};

	// -----вызов процедуры без парметров, с возвратом результата
	template<class R>
	class TDllProc0 : public TDllProc {
	public:
		TDllProc0(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()() {
			typedef R(* TProc)();
			return ((TProc)m_Proc)();
		}
	};

	// -----вызовы с одним параметром
	template<class P1>
	class TDllProcV1 : public TDllProc {
	public:
		TDllProcV1(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1) {
			typedef void(__cdecl*TProc)(P1);
			((TProc)m_Proc)(p1);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1>
	class TDllProc1 : public TDllProc {
	public:
		TDllProc1(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1) {
			typedef R(__cdecl *TProc)(P1);
			return ((TProc)m_Proc)(p1);
		}
	};

	// -----------------------------------------------------------------
	template<class P1>
	class TDllStdProcV1 : public TDllProc {
	public:
		TDllStdProcV1(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1) {
			typedef void(__stdcall*TProc)(P1);
			((TProc)m_Proc)(p1);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1>
	class TDllStdProc1 : public TDllProc {
	public:
		TDllStdProc1(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1) {
			typedef R(__stdcall *TProc)(P1);
			return ((TProc)m_Proc)(p1);
		}
	};

	// вызовы с двумя параметрами

	template<class P1, class P2>
	class TDllProcV2 : public TDllProc {
	public:
		TDllProcV2(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2) {
			typedef void(__cdecl*TProc)(P1, P2);
			((TProc)m_Proc)(p1, p2);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2>
	class TDllProc2 : public TDllProc {
	public:
		TDllProc2(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2) {
			typedef R(__cdecl *TProc)(P1, P2);
			return ((TProc)m_Proc)(p1, p2);
		}
	};

	// -----------------------------------------------------------------
	template<class P1, class P2>
	class TDllStdProcV2 : public TDllProc {
	public:
		TDllStdProcV2(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2) {
			typedef void(__stdcall*TProc)(P1, P2);
			((TProc)m_Proc)(p1, p2);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2>
	class TDllStdProc2 : public TDllProc {
	public:
		TDllStdProc2(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2) {
			typedef R(__stdcall *TProc)(P1, P2);
			return ((TProc)m_Proc)(p1, p2);
		}
	};

	// вызовы с тремя параметрами

	template<class P1, class P2, class P3>
	class TDllProcV3 : public TDllProc {
	public:
		TDllProcV3(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2, P3 p3) {
			typedef void(__cdecl*TProc)(P1, P2, P3);
			((TProc)m_Proc)(p1, p2, p3);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2, class P3>
	class TDllProc3 : public TDllProc {
	public:
		TDllProc3(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2, P3 p3) {
			typedef R(__cdecl *TProc)(P1, P2, P3);
			return ((TProc)m_Proc)(p1, p2, p3);
		}
	};

	// -----------------------------------------------------------------
	template<class P1, class P2, class P3>
	class TDllStdProcV3 : public TDllProc {
	public:
		TDllStdProcV3(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2, P3 p3) {
			typedef void(__stdcall*TProc)(P1, P2, P3);
			((TProc)m_Proc)(p1, p2, p3);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2, class P3>
	class TDllStdProc3 : public TDllProc {
	public:
		TDllStdProc3(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2, P3 p3) {
			typedef R(__stdcall *TProc)(P1, P2, P3);
			return ((TProc)m_Proc)(p1, p2, p3);
		}
	};

	// вызовы с четырьмя параметрами
	template<class P1, class P2, class P3, class P4>
	class TDllProcV4 : public TDllProc {
	public:
		TDllProcV4(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2, P3 p3, P4 p4) {
			typedef void(__cdecl*TProc)(P1, P2, P3, P4);
			((TProc)m_Proc)(p1, p2, p3, p4);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2, class P3, class P4>
	class TDllProc4 : public TDllProc {
	public:
		TDllProc4(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2, P3 p3, P4 p4) {
			typedef R(__cdecl *TProc)(P1, P2, P3, P4);
			return ((TProc)m_Proc)(p1, p2, p3, p4);
		}
	};

	// -----------------------------------------------------------------
	template<class P1, class P2, class P3, class P4>
	class TDllStdProcV4 : public TDllProc {
	public:
		TDllStdProcV4(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		void operator()(P1 p1, P2 p2, P3 p3, P4 p4) {
			typedef void(__stdcall*TProc)(P1, P2, P3, P4);
			((TProc)m_Proc)(p1, p2, p3, p4);
		}
	};
	// -----------------------------------------------------------------

	template<class R, class P1, class P2, class P3, class P4>
	class TDllStdProc4 : public TDllProc {
	public:
		TDllStdProc4(const TDll& a_Dll, const char* a_Name)
			: TDllProc(a_Dll, a_Name) {
		}

		R operator()(P1 p1, P2 p2, P3 p3, P4 p4) {
			typedef R(__stdcall *TProc)(P1, P2, P3, P4);
			return ((TProc)m_Proc)(p1, p2, p3, p4);
		}
	};
	/* */
}
/* */
#endif _DLL_H_
