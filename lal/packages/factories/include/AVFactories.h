/*-----------------------------------------------------------------------
 *
 * File Name: AVFactories.h
 *
 * Author: Finn, L. S.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * AVFactories.h
 *
 * SYNOPSIS
 * #include "AVFactories.h"
 *
 * DESCRIPTION 
 * Provides prototype and status code information for use
 * of CreateVector, CreateArray, DestroyVector and DestroyArray
 *
 * DIAGNOSTICS
 *
 *----------------------------------------------------------------------- 
 */

/* 
 * Header contents go here, in order specified:
 * 
 * 1. Prolog (Comment field with file name, author, revision, etc., as 
 *    specified above)
 * 2. include-loop protection (see below). Note the naming convention!
 */

#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H

/*
 * 3. Includes. This header may include others; if so, they go immediately 
 *    after include-loop protection. Includes should appear in the following 
 *    order: 
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 *    Includes should be double-guarded!
 */

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

/*
 * 4. Header file version string (from CVS; see below). Note the string name. 
 */

NRCSID (AVFACTORIESH, "$Id$");

/*
 * 5. Macros. 
 */

#define CREATEVECTOR_ELENGTH  1
#define CREATEVECTOR_EVPTR    2
#define CREATEVECTOR_EUPTR    4
#define CREATEVECTOR_EMALLOC  8

#define CREATEVECTOR_MSGELENGTH   "Illegal length" /* ELENGTH */
#define CREATEVECTOR_MSGEVPTR    "vector == NULL" /* EVPTR */
#define CREATEVECTOR_MSGEUPTR    "*vector != NULL" /* EUPTR */
#define CREATEVECTOR_MSGEMALLOC   "Malloc failure" /* EMALLOC */

#define DESTROYVECTOR_EVPTR  1
#define DESTROYVECTOR_EUPTR  2
#define DESTROYVECTOR_EDPTR  8

#define DESTROYVECTOR_MSGEVPTR "vector == NULL" /* EVPTR */
#define DESTROYVECTOR_MSGEUPTR "*vector == NULL" /* EUPTR */
#define DESTROYVECTOR_MSGEDPTR "(*vector)->data == NULL" /* EDPTR */

/* 
 * 6. Extern Constant Declarations. These should not be present unless a
 *    specific waiver has been granted. 
 */

/* 
 * 7. Extern Global Variables. These should also not be present unless a 
 *    specific waiver has been granted. 
 */ 

/*
 * 8. Structure, enum, union, etc., typdefs.
 */

/*
 * 9. Functions Declarations (i.e., prototypes).
 */

void CHARCreateVector(Status *, CHARVector **, UINT4);
void I2CreateVector(Status *, INT2Vector **, UINT4);
void I4CreateVector(Status *, INT4Vector **, UINT4);
void I8CreateVector(Status *, INT8Vector **, UINT4);
void U2CreateVector(Status *, UINT2Vector **, UINT4);
void U4CreateVector(Status *, UINT4Vector **, UINT4);
void U8CreateVector(Status *, UINT8Vector **, UINT4);
void CreateVector(Status *, REAL4Vector **, UINT4);
void SCreateVector(Status *, REAL4Vector **, UINT4);
void DCreateVector(Status *, REAL8Vector **, UINT4);
void CCreateVector(Status *, COMPLEX8Vector **, UINT4);
void ZCreateVector(Status *, COMPLEX16Vector **, UINT4);

void I2CreateArray(Status *, INT2Array **, UINT4Vector *);
void I4CreateArray(Status *, INT4Array **, UINT4Vector *);
void I8CreateArray(Status *, INT8Array **, UINT4Vector *);
void U2CreateArray(Status *, UINT2Array **, UINT4Vector *);
void U4CreateArray(Status *, UINT4Array **, UINT4Vector *);
void U8CreateArray(Status *, UINT8Array **, UINT4Vector *);
void CreateArray(Status *, REAL4Array **, UINT4Vector *);
void SCreateArray(Status *, REAL4Array **, UINT4Vector *);
void DCreateArray(Status *, REAL8Array **, UINT4Vector *);
void CCreateArray(Status *, COMPLEX8Array **, UINT4Vector *);
void ZCreateArray(Status *, COMPLEX16Array **, UINT4Vector *);
 
void CHARDestroyVector(Status *, CHARVector **);
void I2DestroyVector(Status *, INT2Vector **);
void I4DestroyVector(Status *, INT4Vector **);
void I8DestroyVector(Status *, INT8Vector **);
void U2DestroyVector(Status *, UINT2Vector **);
void U4DestroyVector(Status *, UINT4Vector **);
void U8DestroyVector(Status *, UINT8Vector **);
void DestroyVector(Status *, REAL4Vector **);
void SDestroyVector(Status *, REAL4Vector **);
void DDestroyVector(Status *, REAL8Vector **);
void CDestroyVector(Status *, COMPLEX8Vector **);
void ZDestroyVector(Status *, COMPLEX16Vector **);

void I2DestroyArray(Status *, INT2Array **);
void I4DestroyArray(Status *, INT4Array **);
void I8DestroyArray(Status *, INT8Array **);
void U2DestroyArray(Status *, UINT2Array **);
void U4DestroyArray(Status *, UINT4Array **);
void U8DestroyArray(Status *, UINT8Array **);
void DestroyArray(Status *, REAL4Array **);
void SDestroyArray(Status *, REAL4Array **);
void DDestroyArray(Status *, REAL8Array **);
void CDestroyArray(Status *, COMPLEX8Array **);
void ZDestroyArray(Status *, COMPLEX16Array **);

#endif

