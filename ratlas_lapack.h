/*
 *  Ratlas_Lapack.h
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#ifndef RATLAS_LAPACK_H
#define RATLAS_LAPACK_H
#include "ratlas.h"

#if defined(__ACCELERATE__) || defined(__VECLIB__)
    #if defined(__LP64__)
        #define LAPACK_INT int
    #else
        #define LAPACK_INT long int
    #endif
    #define LAPACK_DCOMPLEX __CLPK_doublecomplex
#else
    #define LAPACK_INT int
    #define LAPACK_DCOMPLEX Ratlas_Complex
#endif

#ifndef __CLAPACK_H
/* Lapack prototypes */
extern void dgesv_(LAPACK_INT *N, LAPACK_INT *NRHS, double *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, double *B,
        LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void zgesv_(LAPACK_INT *N, LAPACK_INT *NRHS, 
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        LAPACK_DCOMPLEX *B, LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void dposv_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        double *A, LAPACK_INT *LDA, double *B, LAPACK_INT *LDB,
        LAPACK_INT *INFO);
extern void zposv_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_DCOMPLEX *B,
        LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void dsysv_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        double *A, LAPACK_INT *LDA, LAPACK_INT *IPIV, double *B,
        LAPACK_INT *LDB, double *WORK, LAPACK_INT *LWORK,
        LAPACK_INT *INFO);
extern void zsysv_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        LAPACK_DCOMPLEX *B, LAPACK_INT *LDB, LAPACK_DCOMPLEX *WORK,
        LAPACK_INT *LWORK, LAPACK_INT *INFO);
extern void dgetrf_(LAPACK_INT *N, LAPACK_INT *M, double *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, LAPACK_INT
        *INFO);
extern void zgetrf_(LAPACK_INT *N, LAPACK_INT *M, LAPACK_DCOMPLEX
        *A, LAPACK_INT *LDA, LAPACK_INT *IPIV, LAPACK_INT
        *INFO);
extern void dpotrf_(char *UPLO, LAPACK_INT *N, double *A,
        LAPACK_INT *LDA, LAPACK_INT *INFO);
extern void zpotrf_(char *UPLO, LAPACK_INT *N, LAPACK_DCOMPLEX *A,
        LAPACK_INT *LDA, LAPACK_INT *INFO);
extern void dsytrf_(char *UPLO, LAPACK_INT *N, double *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, double *WORK,
        LAPACK_INT *LWORK, LAPACK_INT *INFO);
extern void zsytrf_(char *UPLO, LAPACK_INT *N, LAPACK_DCOMPLEX *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, LAPACK_DCOMPLEX *WORK,
        LAPACK_INT *LWORK, LAPACK_INT *INFO);
extern void dgetrs_( char *TRANS, LAPACK_INT *N, LAPACK_INT
        *NRHS, double *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        double *B, LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void zgetrs_(char *TRANS, LAPACK_INT *N, LAPACK_INT *NRHS,
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        LAPACK_DCOMPLEX *B, LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void dpotrs_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        double *A, LAPACK_INT *LDA, double *B, LAPACK_INT *LDB,
        LAPACK_INT *INFO);
extern void zpotrs_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_DCOMPLEX *B,
        LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void dsytrs_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        double *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        double *B, LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void zsytrs_(char *UPLO, LAPACK_INT *N, LAPACK_INT *NRHS,
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, LAPACK_INT *IPIV,
        LAPACK_DCOMPLEX *B, LAPACK_INT *LDB, LAPACK_INT *INFO);
extern void dgetri_(LAPACK_INT *N, double *A, LAPACK_INT *LDA,
        LAPACK_INT *IPIV, double *WORK, LAPACK_INT *LWORK,
        LAPACK_INT *INFO);
extern void zgetri_(LAPACK_INT *N, LAPACK_DCOMPLEX *A, LAPACK_INT
        *LDA, LAPACK_INT *IPIV, LAPACK_DCOMPLEX *WORK, LAPACK_INT
        *LWORK, LAPACK_INT *INFO);
extern void dpotri_(char *UPLO, LAPACK_INT *N, double *A,
        LAPACK_INT *LDA, LAPACK_INT *INFO);
extern void zpotri_(char *UPLO, LAPACK_INT *N, LAPACK_DCOMPLEX *A,
        LAPACK_INT *LDA, LAPACK_INT *INFO);
extern void dsytri_(char *UPLO, LAPACK_INT *N, double *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, double *WORK,
        LAPACK_INT *INFO);
extern void zsytri_(char *UPLO, LAPACK_INT *N, LAPACK_DCOMPLEX *A,
        LAPACK_INT *LDA, LAPACK_INT *IPIV, LAPACK_DCOMPLEX *WORK,
        LAPACK_INT *INFO);
extern void dgesvd_(char *JOBU, char *JOBVT, LAPACK_INT *M, LAPACK_INT *N, 
        double *A, LAPACK_INT *LDA,  double *S,  double *U,  LAPACK_INT *LDU, 
        double *VT,  LAPACK_INT *LDVT, double *WORK, LAPACK_INT *LWORK, 
        LAPACK_INT *INFO );
extern void zgesvd_(char *JOBU, char *JOBVT, LAPACK_INT *M, LAPACK_INT *N, 
        LAPACK_DCOMPLEX *A, LAPACK_INT *LDA, double *S, LAPACK_DCOMPLEX *U,
        LAPACK_INT *LDU, LAPACK_DCOMPLEX *VT, LAPACK_INT *LDVT,
        LAPACK_DCOMPLEX *WORK, LAPACK_INT *LWORK, double *RWORK,
        LAPACK_INT *INFO);
#endif /* ifdef __CLAPACK_H */

VALUE ratlas_gesv_bang(VALUE self, VALUE astor, VALUE bstor);
VALUE ratlas_posv_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor);
VALUE ratlas_sysv_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor);

VALUE ratlas_getrf_bang(VALUE self, VALUE astor);
VALUE ratlas_potrf_bang(VALUE self, VALUE uplo, VALUE astor);
VALUE ratlas_sytrf_bang(VALUE self, VALUE uplo, VALUE astor);

VALUE ratlas_getrs_bang(VALUE self, VALUE trans, VALUE astor, VALUE rowperm,
        VALUE bstor);
VALUE ratlas_potrs_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor);
VALUE ratlas_sytrs_bang(VALUE self, VALUE uplo, VALUE astor, VALUE rowperm,
        VALUE bstor);

VALUE ratlas_getri_bang(VALUE self, VALUE astor, VALUE rowperm);
VALUE ratlas_potri_bang(VALUE self, VALUE uplo, VALUE astor);
VALUE ratlas_sytri_bang(VALUE self, VALUE uplo, VALUE astor, VALUE rowperm);
VALUE ratlas_gesvd_bang(VALUE self, VALUE astor);

#endif /* ifdef RATLAS_LAPACK_H */
