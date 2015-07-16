// SICx.h

// Copyright © 2000-3000, Andrey A. Meshkov (AL-CHEMIST)
// All rights reserved
//
// http://maalchemist.narod.ru
// maalchemist@yandex.ru
// maalchemist@gmail.com

#ifndef _SICX_H_
#define _SICX_H_

#pragma once

#define SICXAPI __stdcall

// ===========================================================================
//
// SICx Defs
//
// ===========================================================================

// Compiler options flags

#define SIC_OPT_FLAG_OPTIMIZATION       0x01 // 00000001
#define SIC_OPT_FLAG_STACK_FRAME        0x02 // 00000010
#define SIC_OPT_FLAG_STACK_ALIGN        0x04 // 00000100
#define SIC_OPT_FLAG_FPU_FRAME          0x08 // 00001000
#define SIC_OPT_FLAG_CODE_ALIGN         0x10 // 00010000
#define SIC_OPT_FLAG_CODE_REALLOC       0x20 // 00100000

// Compiler options

// --------------------------------
// Code optimization
// --------------------------------
#define SIC_OPT_OPTIMIZATION    SIC_OPT_FLAG_OPTIMIZATION

// Stack frame
// --------------------------------
#define SIC_OPT_STACK_FRAME     SIC_OPT_FLAG_STACK_FRAME

// 16-byte stack alignment
// --------------------------------
#define SIC_OPT_STACK_ALIGN     SIC_OPT_FLAG_STACK_ALIGN

// Mask all FPU exceptions
// Return FPU exception flags in EAX
// Use this option to avoid exception raising on FPU errors
// --------------------------------
#define SIC_OPT_FPU_FRAME       SIC_OPT_FLAG_FPU_FRAME

// Build mode
// x32: 4-byte code alignment
// x64: 8-byte code alignment
// --------------------------------
#define SIC_OPT_CODE_ALIGN      SIC_OPT_FLAG_CODE_ALIGN

// Code reallocation
// --------------------------------
#define SIC_OPT_CODE_REALLOC    SIC_OPT_FLAG_CODE_REALLOC

// Compiler default options x32
// --------------------------------
// Code optimization
// Stack frame
// 4-byte code alignment
// --------------------------------
#define SIC_OPT_DEFAULT_X32   ( SIC_OPT_OPTIMIZATION | \
                                SIC_OPT_STACK_FRAME  | \
                                SIC_OPT_CODE_ALIGN   | \
                                SIC_OPT_CODE_REALLOC )

// Compiler default options x64
// --------------------------------
// Code optimization
// Stack frame
// 16-byte stack alignment
// 8-byte code alignment
// --------------------------------
#define SIC_OPT_DEFAULT_X64   ( SIC_OPT_OPTIMIZATION | \
                                SIC_OPT_STACK_FRAME  | \
                                SIC_OPT_STACK_ALIGN  | \
                                SIC_OPT_CODE_ALIGN   | \
								SIC_OPT_CODE_REALLOC )

// ===========================================================================
//
// SICx Structures
//
// ===========================================================================

//
// use 1 byte packing for the data structures
//
#include <pshpack1.h>

//
// Main SIC structure (x64)
//
typedef struct SIC_Data64 {
	UINT64 fdata; // Function data segment offset
	UINT64 cdata; // Constant data segment offset
	UINT64 vdata; // Variable data segment offset
	UINT64 rdata; // Runtime data segment offset
	UINT64 code; // Code segment offset
	UINT64 data; // Data segment offset
	UINT64 stack; // Stack segment offset
	UINT64 entry; // Entry point
	UINT32 size; // Code size
	UINT32 coops; // Compiler options
	UINT32 tokens; // Scanned tokens count
	UINT32 rpn; // Rpn array item count
	UINT32 fcount; // Functions count
	UINT32 ccount; // Constants count
	UINT32 vcount; // Variables count
	UINT32 rcount; // Runtimes count
	UINT32 ccurs; // Current string cursor
	UINT32 pcurs; // Previous string cursor
	UINT32 ecode; // Error code
	UINT32 rcode; // Return code
	double value; // Return value
} TSIC_Data64, *PSIC_Data64;

//
// Main SIC structure (x86)
//
typedef struct SIC_Data32 {
	UINT32 fdata; // Function data segment offset
	UINT32 cdata; // Constant data segment offset
	UINT32 vdata; // Variable data segment offset
	UINT32 rdata; // Runtime data segment offset
	UINT32 code; // Code segment offset
	UINT32 data; // Data segment offset
	UINT32 stack; // Stack segment offset
	UINT32 entry; // Entry point
	UINT32 size; // Code size
	UINT32 coops; // Compiler options
	UINT32 tokens; // Scanned tokens count
	UINT32 rpn; // Rpn array item count
	UINT32 fcount; // Functions count
	UINT32 ccount; // Constants count
	UINT32 vcount; // Variables count
	UINT32 rcount; // Runtimes count
	UINT32 ccurs; // Current string cursor
	UINT32 pcurs; // Previous string cursor
	UINT32 ecode; // Error code
	UINT32 rcode; // Return code
	double value; // Return value
} TSIC_Data32, *PSIC_Data32;

#ifdef _WIN64
#define TSIC_Data TSIC_Data64
#else
#define TSIC_Data TSIC_Data32
#endif

//
// Common table header
//
typedef struct SIC_TableHeader {
	INT32 icount; // Item count
	INT32 mcount; // Item max count
	INT32 tisize; // Table item size
	INT32 tnsize; // Table item name size
	INT32 titype; // Table item type
	// = 1 - functions
	// = 2 - constants
	// = 3 - variables
	INT32 oooooo; // Padding
} TSIC_TableHeader, *PSIC_TableHeader;

//
// Function table item
//
typedef struct SIC_FunItem {
	CHAR name[52]; // Function name (zero terminated)
	INT16 acount; // Function argument count
	INT16 cosize; // Function code size or flags
	UINT64 offset; // Function offset
} TSIC_FunItem, *PSIC_FunItem;

//
// Function table
//
typedef struct SIC_FunTable {
	TSIC_TableHeader header; // Function table header
	TSIC_FunItem items[1]; // Function item list
} TSIC_FunTable, *PSIC_FunTable;

//
// Constant table item
//
typedef struct SIC_ConItem {
	CHAR name[54]; // Constant name (zero terminated)
	INT16 datype; // Constant data type
	double value; // Constant value
} TSIC_ConItem, *PSIC_ConItem;

//
// Constant table
//
typedef struct SIC_ConTable {
	TSIC_TableHeader header; // Constant table header
	TSIC_ConItem items[1]; // Constant item list
} TSIC_ConTable, *PSIC_ConTable;

//
// Variable table item
//
typedef struct SIC_VarItem {
	CHAR name[54]; // Variable name (zero terminated)
	INT16 datype; // Variable data type
	UINT64 offset; // Variable offset
} TSIC_VarItem, *PSIC_VarItem;

//
// Variable table
//
typedef struct SIC_VarTable {
	TSIC_TableHeader header; // Variable table header
	TSIC_VarItem items[1]; // Variable item list
} TSIC_VarTable, *PSIC_VarTable;

//
// end of  using 1 byte packing for the data structures
//
#include <poppack.h>
// ===========================================================================
//
// SICx Functions
//
// ===========================================================================

VOID SICXAPI sic_cretab(VOID);
//
// Create global tables
//

VOID SICXAPI sic_fretab(VOID);
//
// Destroy global tables
//

DWORD SICXAPI sic_funtac(VOID);
//
// Create global function table
// Assign table header and add built-in functions
//
// <- Result : function table item count or zero on error
//

VOID SICXAPI sic_funtaf(VOID);
//
// Destroy global function table
//

DWORD SICXAPI sic_funloa(VOID);
//
// Load external user defined functions
//
// <- Result : external function item count
//

VOID SICXAPI sic_funulo(VOID);
//
// Unload external user defined functions
//

DWORD SICXAPI sic_contac(VOID);
//
// Create global constant table
// Assign table header and add predefined constants
//
// <- Result : constant table item count or zero on error
//

VOID SICXAPI sic_contaf(VOID);
//
// Destroy global constant table
//

DWORD SICXAPI sic_conloa(VOID);
//
// Load external user defined constants
//
// <- Result : external constant item count
//

VOID SICXAPI sic_conulo(VOID);
//
// Unload external user defined constants
//

DWORD SICXAPI sic_vartac(VOID);
//
// Create global variable table
// Assign table header and add predefined variables
//
// <- Result : variable table item count or zero on error
//

VOID SICXAPI sic_vartaf(VOID);
//
// Destroy global variable table
//

DWORD SICXAPI sic_varloa(VOID);
//
// Load external user defined variables
//
// <- Result : external variable item count
//

VOID SICXAPI sic_varulo(VOID);
//
// Unload external user defined variables
//

DWORD SICXAPI sic_runtac(VOID);
//
// Create global runtime table
//
// <- Result : function table item count
//

VOID SICXAPI sic_runtaf(VOID);
//
// Destroy global runtime table
//

VOID SICXAPI sic_init(LPVOID ASic);
//
// Allocate memory for ASic data segments
// Assign table headers
//
// -> ASic : TSIC_Data structure offset
//

VOID SICXAPI sic_done(LPVOID ASic);
//
// Free memory previously allocated for ASic data and code segments
//
// -> ASic : TSIC_Data structure offset
//

INT SICXAPI sic_afun(LPVOID ASic, LPCSTR AFuname, LPVOID AOffset, INT16 AACount,
	UINT16 AFlags);
//
// Add|set user defined global or local function ( AFuname ) and assign data ( AOffset, AACount, AFlags )
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AFuname : function name
// -> AOffset : function offset
// -> AACount : function argument count
// -> AFlags  : function flags
// x32 : count of 4-byte arguments or
// 0x??FF for stdcall function
// x64 : 4-bit mask for integer and pointer arguments
// ex.: 0010 - second argument is integer or pointer
// 0101 - first and third arguments are integer or pointer
// <- Result : function index or negative value on error
//

INT SICXAPI sic_refun(LPVOID ASic, LPCSTR AFuname, LPCSTR AOrgname,
	BOOL AInvalid);
//
// Rename global or local function
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic     : TSIC_Data structure offset
// -> AFuname  : function new name
// -> AOrgname : function original name
// -> AInvalid : invalidate original function if TRUE
// <- Result   : function index or negative value on error
//

INT SICXAPI sic_aconf(LPVOID ASic, LPCSTR AConame, double AValue);
//
// Add|set user defined global or local float constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_aconi(LPVOID ASic, LPCSTR AConame, INT_PTR AValue);
//
// Add|set user defined global or local integer constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_acons(LPVOID ASic, LPCSTR AConame, LPCSTR AValue);
//
// Add|set user defined global or local string constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_aconp(LPVOID ASic, LPCSTR AConame, LPVOID AValue);
//
// Add|set user defined global or local pointer constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_aconpf(LPVOID ASic, LPCSTR AConame, LPVOID AValue);
//
// Add|set user defined global or local float pointer constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_aconpi(LPVOID ASic, LPCSTR AConame, LPVOID AValue);
//
// Add|set user defined global or local integer pointer constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_aconps(LPVOID ASic, LPCSTR AConame, LPVOID AValue);
//
// Add|set user defined global or local string pointer constant ( AConame ) and assign data ( AValue )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// -> AValue  : constant value
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_recon(LPVOID ASic, LPCSTR AConame, LPCSTR AOrgname,
	BOOL AInvalid);
//
// Rename global or local constant
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic     : TSIC_Data structure offset
// -> AConame  : constant new name
// -> AOrgname : constant original name
// -> AInvalid : invalidate original constant if TRUE
// <- Result   : constant index or negative value on error
//

INT SICXAPI sic_avarf(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local float variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avari(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local integer variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avars(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local string variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avarp(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local pointer variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avarpf(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local float pointer variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avarpi(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local integer pointer variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_avarps(LPVOID ASic, LPCSTR AVaname, LPVOID AOffset);
//
// Add|set user defined global or local string pointer variable ( AVaname ) and assign data ( AOffset )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// -> AOffset : variable offset
// <- Result  : variable index or negative value on error
//

INT SICXAPI sic_revar(LPVOID ASic, LPCSTR AVaname, LPCSTR AOrgname,
	BOOL AInvalid);
//
// Rename global or local variable
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic     : TSIC_Data structure offset
// -> AVaname  : variable new name
// -> AOrgname : variable original name
// -> AInvalid : invalidate original variable if TRUE
// <- Result   : variable index or negative value on error
//

INT SICXAPI sic_invaf(LPVOID ASic, LPCSTR AFuname);
//
// Invalidate global or local function ( AFuname )
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AFuname : function name
// <- Result  : function index or negative value on error
//

INT SICXAPI sic_invac(LPVOID ASic, LPCSTR AConame);
//
// Invalidate global or local constant ( AConame )
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// <- Result  : constant index or negative value on error
//

INT SICXAPI sic_invav(LPVOID ASic, LPCSTR AVaname);
//
// Invalidate global or local variable ( AVaname )
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// <- Result  : variable index or negative value on error
//

VOID SICXAPI sic_patab(LPVOID ASic);
//
// Pack global or local tables ( remove deleted items, decrease and fix table size )
// Deal with global tables if ASic = null or ASic.?data = null
// !!! You can`t add any item to the packed table
//
// -> ASic : TSIC_Data structure offset
//

DWORD SICXAPI sic_pafut(LPVOID ASic);
//
// Pack global or local function table ( remove deleted items, decrease and fix table size )
// Deal with global table if ASic = null or ASic.fdata = null
// !!! You can`t add any item to the packed table
//
// -> ASic   : TSIC_Data structure offset
// <- Result : table item count or zero on error
//

DWORD SICXAPI sic_pacot(LPVOID ASic);
//
// Pack global or local constant table ( remove deleted items, decrease and fix table size )
// Deal with global table if ASic = null or ASic.cdata = null
// !!! You can`t add any item to the packed table
//
// -> ASic   : TSIC_Data structure offset
// <- Result : table item count or zero on error
//

DWORD SICXAPI sic_pavat(LPVOID ASic);
//
// Pack global or local variable table ( remove deleted items, decrease and fix table size )
// Deal with global table if ASic = null or ASic.vdata = null
// !!! You can`t add any item to the packed table
//
// -> ASic   : TSIC_Data structure offset
// <- Result : table item count or zero on error
//

LPVOID SICXAPI sic_gefut(LPVOID ASic);
//
// Get global or local function table offset
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic   : TSIC_Data structure offset
// <- Result : function table offset
//

DWORD SICXAPI sic_gefuc(LPVOID ASic);
//
// Get global or local function item count
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic   : TSIC_Data structure offset
// <- Result : function item count
//

INT SICXAPI sic_gefui(LPVOID ASic, INT AIndex, LPVOID AItem);
//
// Get global or local function item
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic   : TSIC_Data structure offset
// -> AIndex : function index
// <- AItem  : function item
// <- Result : function index or negative value on error
//

INT SICXAPI sic_gefun(LPVOID ASic, LPCSTR AFuname);
//
// Get global or local item ( AFuname ) index in function table
// Deal with global table if ASic = null or ASic.fdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AFuname : function name
// <- Result  : function index or negative value on error
//

LPVOID SICXAPI sic_gecot(LPVOID ASic);
//
// Get global or local constant table offset
// Deal with global table if ASic = null or ASic.cdata = nil
//
// -> ASic   : TSIC_Data structure offset
// <- Result : constant table offset
//

DWORD SICXAPI sic_gecoc(LPVOID ASic);
//
// Get global or local constant item count
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic   : TSIC_Data structure offset
// <- Result : constant item count
//

INT SICXAPI sic_gecoi(LPVOID ASic, INT AIndex, LPVOID AItem);
//
// Get global or local constant item
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic   : TSIC_Data structure offset
// -> AIndex : constant index
// <- AItem  : constant item
// <- Result : constant index or negative value on error
//

INT SICXAPI sic_gecon(LPVOID ASic, LPCSTR AConame);
//
// Get global or local item ( AConame ) index in constant table
// Deal with global table if ASic = null or ASic.cdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AConame : constant name
// <- Result  : constant index or negative value on error
//

LPVOID SICXAPI sic_gevat(LPVOID ASic);
//
// Get global or local variable table offset
// Deal with global table if ASic = null or ASic.vdata = nil
//
// -> ASic   : TSIC_Data structure offset
// <- Result : variable table offset
//

DWORD SICXAPI sic_gevac(LPVOID ASic);
//
// Get global or local variable item count
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic   : TSIC_Data structure offset
// <- Result : variable item count
//

INT SICXAPI sic_gevai(LPVOID ASic, INT AIndex, LPVOID AItem);
//
// Get global or local variable item
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic   : TSIC_Data structure offset
// -> AIndex : variable index
// <- AItem  : variable item
// <- Result : variable index or negative value on error
//

INT SICXAPI sic_gevar(LPVOID ASic, LPCSTR AVaname);
//
// Get global or local item ( AVaname ) index in variable table
// Deal with global table if ASic = null or ASic.vdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> AVaname : variable name
// <- Result  : variable index or negative value on error
//

LPVOID SICXAPI sic_gerut(LPVOID ASic);
//
// Get global or local runtime table offset
// Deal with global table if ASic = null or ASic.rdata = nil
//
// -> ASic   : TSIC_Data structure offset
// <- Result : runtime table offset
//

DWORD SICXAPI sic_geruc(LPVOID ASic);
//
// Get global or local runtime item count
// Deal with global table if ASic = null or ASic.rdata = null
//
// -> ASic   : TSIC_Data structure offset
// <- Result : runtime item count
//

INT SICXAPI sic_gerui(LPVOID ASic, INT AIndex, LPVOID AItem);
//
// Get global or local runtime item
// Deal with global table if ASic = null or ASic.rdata = null
//
// -> ASic   : TSIC_Data structure offset
// -> AIndex : runtime index
// <- AItem  : runtime item
// <- Result : runtime index or negative value on error
//

INT SICXAPI sic_gerun(LPVOID ASic, LPCSTR ARuname);
//
// Get global or local item ( ARuname ) index in runtime table
// Deal with global table if ASic = null or ASic.rdata = null
//
// -> ASic    : TSIC_Data structure offset
// -> ARuname : runtime name
// <- Result  : runtime index or negative value on error
//

DWORD SICXAPI sic_compile(LPVOID ASic, LPCSTR S, DWORD ASop);
//
// Allocate memory for ASic code segment and compile string ( S )
//
// -> ASic : TSIC_Data structure offset
// -> S    : string to compile
// -> ASop : sic compiler options
//
// <- Result      : generated code size or zero on error
// <- ASic.coops  : actual compiler options
// <- ASic.tokens : scanned tokens count
// <- ASic.ccurs  : current string cursor
// <- ASic.pcurs  : previous string cursor
//

DWORD SICXAPI sic_build(LPVOID ASic, LPCSTR S, DWORD ASop);
//
// Allocate memory for ASic code segment and compile string ( S )
// <S> can be complex expression with ';' as delimiter
//
// -> ASic : TSIC_Data structure offset
// -> S    : string to compile
// -> ASop : sic compiler options
//
// <- Result      : generated code size or zero on error
// <- ASic.coops  : actual compiler options
// <- ASic.tokens : scanned tokens count
// <- ASic.ccurs  : current string cursor
// <- ASic.pcurs  : previous string cursor
//

double SICXAPI sic_exec(LPVOID ASic, LPDWORD AError);
//
// Execute code
//
// -> ASic   : TSIC_Data structure offset
//
// <- Result : result or zero on error
// <- AError : error code or zero on success
// AError = 0x00010000 -> null TSIC_Data structure offset
// AError = 0x00020000 -> null code segment offset
// AError = 0x00040000 -> invalid code size
// AError = 0xFF000000 -> EXCEPTION (TODO?)
//
// if compiler option SIC_OPT_FLAG_FPU_EFLAGS is defined
// AError and 0x00000008 = OE flag
// AError and 0x00000004 = ZE flag
// AError and 0x00000001 = IE flag
//

VOID SICXAPI sic_call(LPVOID ASic);
//
// Execute code
//
// -> ASic       : TSIC_Data structure offset
// <- ASic.value : result
//

double SICXAPI sic_cexec(LPVOID ASic, LPCSTR S, LPDWORD ASop, LPDWORD AError);
//
// Compile & execute string
//
// -> ASic   : TSIC_Data structure offset
// -> S      : string to compile
// -> ASop   : sic compiler options offset
//
// <- Result : result or zero on error
// <- ASop   : actual compiler options
// <- AError : error code or zero on success
// AError = 0x00000100 -> compiler error
// AError = 0xFF000000 -> EXCEPTION (TODO?)
//
// if compiler option SIC_OPT_FLAG_FPU_EFLAGS is defined
// AError and 0x00000008 = OE flag
// AError and 0x00000004 = ZE flag
// AError and 0x00000001 = IE flag
//

double SICXAPI sic_bexec(LPVOID ASic, LPCSTR S, LPDWORD ASop, LPDWORD AError);
//
// Compile & execute string
// <S> can be complex expression with ';' as delimiter
//
// -> ASic   : TSIC_Data structure offset
// -> S      : string to compile
// -> ASop   : sic compiler options offset
//
// <- Result : result or zero on error
// <- ASop   : actual compiler options
// <- AError : error code or zero on success
// AError = 0x00000100 -> compiler error
// AError = 0xFF000000 -> EXCEPTION (TODO?)
//
// if compiler option SIC_OPT_FLAG_FPU_EFLAGS is defined
// AError and 0x00000008 = OE flag
// AError and 0x00000004 = ZE flag
// AError and 0x00000001 = IE flag
//

double SICXAPI sic_scexec(LPCSTR S, LPDWORD ASop, LPDWORD AError);
//
// Compile & execute string
//
// -> S      : string to compile
// -> ASop   : sic compiler options offset
//
// <- Result : result or zero on error
// <- ASop   : actual compiler options
// <- AError : error code or zero on success
// AError = 0x00000100 -> compiler error
// AError = 0xFF000000 -> EXCEPTION (TODO?)
//
// if compiler option SIC_OPT_FLAG_FPU_EFLAGS is defined
// AError and 0x00000008 = OE flag
// AError and 0x00000004 = ZE flag
// AError and 0x00000001 = IE flag
//

double SICXAPI sic_sbexec(LPCSTR S, LPDWORD ASop, LPDWORD AError);
//
// Compile & execute string
// <S> can be complex expression with ';' as delimiter
//
// -> S      : string to compile
// -> ASop   : sic compiler options offset
//
// <- Result : result or zero on error
// <- ASop   : actual compiler options
// <- AError : error code or zero on success
// AError = 0x00000100 -> compiler error
// AError = 0xFF000000 -> EXCEPTION (TODO?)
//
// if compiler option SIC_OPT_FLAG_FPU_EFLAGS is defined
// AError and 0x00000008 = OE flag
// AError and 0x00000004 = ZE flag
// AError and 0x00000001 = IE flag
//

INT SICXAPI sic_va_count(VOID);
//
// Variable argument count
//
// <- Result : variable argument count
//

#endif // _SICX_H_
