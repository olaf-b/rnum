/*
 *  Ratlas_Lapack.c
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#include "ratlas_lapack.h"

/* LAPACK functions */

/**
 * computes the solution to a real system of linear equations
 *      A * X = B,
 * where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 * 
 * The LU decomposition with partial pivoting and row interchanges is
 * used to factor A as
 *      A = P * L * U,
 * where P is a permutation matrix, L is unit lower triangular, and U is
 * upper triangular.  The factored form of A is then used to solve the
 * system of equations A * X = B. */
VALUE ratlas_gesv_bang(VALUE self, VALUE astor, VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    int i;
    LAPACK_INT n, nrhs, lda, ldb, info, *ipiv;
    VALUE rowperm;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    if(amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Matrix must be square.");
    if(amat->ncol != bmat->nrow)
        rb_raise(rb_eArgError, "Column-row mismatch for solving Ax=b.");

    n = (LAPACK_INT) amat->nrow;
    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;

    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            dgesv_(&n, &nrhs, amat->data, &lda, ipiv, bmat->data, &ldb, &info);
            break;
        case RATLAS_DCOMPLEX:
            zgesv_(&n, &nrhs, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv, 
                    (LAPACK_DCOMPLEX *) bmat->data, &ldb, &info);
            break;
    }
    rowperm = rb_ary_new2(n);
    for (i = 0; i < n; i++)
    {
        rb_ary_push(rowperm, INT2FIX(ipiv[i]));
    }

    return rb_ary_new3(4, INT2FIX(info), astor, rowperm, bstor);
}



/**
 * POSV computes the solution to a real system of linear equations
 *      A * X = B,
 *   where A is an N-by-N symmetric positive definite matrix and X and B
 *   are N-by-NRHS matrices. */
VALUE ratlas_posv_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    LAPACK_INT n, nrhs, lda, ldb, info;
    char uploflag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    if(amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Matrix must be square.");
    if(amat->ncol != bmat->nrow)
        rb_raise(rb_eArgError, "Column-row mismatch for solving Ax=b.");

    n = (LAPACK_INT) amat->nrow;
    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;
    uploflag = ratlas_blasflag2lapack(FIX2INT(uplo));

    switch (amat->type) {
        case RATLAS_DFLOAT:
            dposv_(&uploflag, &n, &nrhs, amat->data, &lda, 
                    bmat->data, &ldb, &info);
            break;
        case RATLAS_DCOMPLEX:
            zposv_(&uploflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data, &lda, 
                    (LAPACK_DCOMPLEX *) bmat->data, &ldb, &info);
            break;
    }
    return rb_ary_new3(3, INT2FIX(info), astor, bstor);
}


/**
 * SYSV computes the solution to a real system of linear equations
 * *     A * X = B,
 * *  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
 * *  matrices. */
VALUE ratlas_sysv_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    int i;
    LAPACK_INT n, nrhs, lda, ldb, info, *ipiv, lwork=-1;
    char uploflag;
    double qwork, *work; 
    LAPACK_DCOMPLEX cqwork, *cwork; 
    VALUE rowperm;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    if(amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Matrix must be square.");
    if(amat->ncol != bmat->nrow)
        rb_raise(rb_eArgError, "Column-row mismatch for solving Ax=b.");

    n = (LAPACK_INT) amat->nrow;
    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;
    uploflag = ratlas_blasflag2lapack(FIX2INT(uplo));
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);

    switch (amat->type) {
        case RATLAS_DFLOAT:
            dsysv_(&uploflag, &n, &nrhs, amat->data, &lda, ipiv, 
                    bmat->data, &ldb, &qwork, &lwork, 
                    &info);
            lwork = (LAPACK_INT) qwork; 
            work = ALLOCA_N(double, lwork);
            dsysv_(&uploflag, &n, &nrhs, amat->data, &lda, ipiv, 
                    bmat->data, &ldb, work, &lwork, &info);
            break;
        case RATLAS_DCOMPLEX:
            zsysv_(&uploflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data,
                    &lda, ipiv, (LAPACK_DCOMPLEX *) bmat->data, &ldb, 
                    &cqwork, &lwork, &info);
            lwork = (LAPACK_INT) cqwork.r; 
            cwork = ALLOCA_N(LAPACK_DCOMPLEX, lwork);
            zsysv_(&uploflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data,
                    &lda, ipiv, (LAPACK_DCOMPLEX *) bmat->data, &ldb, cwork, 
                    &lwork, &info);
            break;
    }
    
    rowperm = rb_ary_new2(amat->nrow);
    for (i = 0; i < amat->nrow; i++)
    {
        rb_ary_push(rowperm, INT2FIX(ipiv[i]));
    }
    
    return rb_ary_new3(4, INT2FIX(info), astor, rowperm, bstor);
}




/**
 * LU factorization of a general M-by-N matrix A
 *   using partial pivoting with row interchanges.
 * 
 *   The factorization has the form
 *      A = P * L * U
 *   where P is a permutation matrix, L is lower triangular with unit
 *   diagonal elements (lower trapezoidal if m > n), and U is upper
 *   triangular (upper trapezoidal if m < n). */
VALUE ratlas_getrf_bang(VALUE self, VALUE astor)
{
    Ratlas_Matrix *amat;
    int i;
    LAPACK_INT *ipiv, info, lda, n, m;
    VALUE rowperm;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    lda = (LAPACK_INT) amat->nrow;
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    m = (LAPACK_INT) amat->nrow;
    n = (LAPACK_INT) amat->ncol;
    switch (amat->type) {
        case RATLAS_DFLOAT:
            dgetrf_(&m, &n, amat->data, &lda, ipiv, &info);
            break;
        case RATLAS_DCOMPLEX:
            zgetrf_(&m, &n, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv, &info);
            break;
    }

    rowperm = rb_ary_new2(amat->nrow);
    for (i = 0; i < amat->nrow; i++)
    {
        rb_ary_push(rowperm, INT2FIX(ipiv[i]));
    }
    return rb_ary_new3(3, INT2FIX(info), astor, rowperm);
}




/**
 * POTRF computes the Cholesky factorization of a real symmetric
 *   positive definite matrix A.
 * 
 *   The factorization has the form
 *      A = U**T * U,  if UPLO = 'U', or
 *      A = L  * L**T,  if UPLO = 'L',
 *   where U is an upper triangular matrix and L is lower triangular.*/
VALUE ratlas_potrf_bang(VALUE self, VALUE uplo, VALUE astor)
{
    Ratlas_Matrix *amat;
    int uploint;
    LAPACK_INT info, lda, n;
    char uploflag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    lda = (LAPACK_INT) amat->nrow;
    n = (LAPACK_INT) amat->ncol;
    uploint = FIX2INT(uplo);
    uploflag = ratlas_blasflag2lapack(uploint);

    switch (amat->type) {
        case RATLAS_DFLOAT:
            dpotrf_(&uploflag, &n, amat->data, &lda, &info);
            break;
        case RATLAS_DCOMPLEX:
            zpotrf_(&uploflag, &n, (LAPACK_DCOMPLEX *) amat->data, &lda, &info);
            break;
    }
    ratlas_matrix_uplo_nullify(CblasNonUnit, uploint, amat);

    return rb_ary_new3(2, INT2FIX(info), astor);
}




/**
 * SYTRF computes the factorization of a real symmetric matrix A using
 * *  the Bunch-Kaufman diagonal pivoting method.  The form of the
 * *  factorization is
 * *
 * *     A = U*D*U**T  or  A = L*D*L**T
 * *
 * *  where U (or L) is a product of permutation and unit upper (lower)
 * *  triangular matrices, and D is symmetric and block diagonal with
 * *  1-by-1 and 2-by-2 diagonal blocks. 
 * *  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the
 * *  upper triangular part of the matrix A, and the strictly lower triangular
 * *  part of A is not referenced.  If UPLO = 'L', the leading N-by-N lower
 * *  triangular part of A contains the lower triangular part of the matrix A,
 * *  and the strictly upper triangular part of A is not referenced.
 * * 
 * *  On exit, A contains the block diagonal matrix D and the multipliers used
 * *  to obtain the factor U or L
 * * 
 * *  If UPLO = 'U', then A = U*D*U', where
 * *     U = P(n)*U(n)* ... *P(k)U(k)* ...,
 * *  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
 * *  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 * *  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 * *  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
 * *  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 * *
 * *             (   I    v    0   )   k-s
 * *     U(k) =  (   0    I    0   )   s
 * *             (   0    0    I   )   n-k
 * *                k-s   s   n-k
 * *
 * *  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
 * *  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
 * *  and A(k,k), and v overwrites A(1:k-2,k-1:k).
 * *
 * *  If UPLO = 'L', then A = L*D*L', where
 * *     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
 * *  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
 * *  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 * *  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 * *  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
 * *  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 * *
 * *             (   I    0     0   )  k-1
 * *     L(k) =  (   0    I     0   )  s
 * *             (   0    v     I   )  n-k-s+1
 * *                k-1   s  n-k-s+1
 * *
 * *  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
 * *  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
 * *  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
 * */
VALUE ratlas_sytrf_bang(VALUE self, VALUE uplo, VALUE astor)
{
    Ratlas_Matrix *amat;
    int i, uploint;
    LAPACK_INT n, lda, info, *ipiv, lwork=-1;
    char uploflag;
    double qwork, *work; 
    LAPACK_DCOMPLEX cqwork, *cwork; 
    VALUE rowperm;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);

    if(amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Matrix must be square.");

    lda = (LAPACK_INT) amat->nrow;
    n = (LAPACK_INT) amat->ncol;
    uploint = FIX2INT(uplo);
    uploflag = ratlas_blasflag2lapack(uploint);
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            dsytrf_(&uploflag, &n, (double *) amat->data, &lda, ipiv, 
                    &qwork, &lwork, &info);
            lwork = (LAPACK_INT) qwork; 
            work = ALLOCA_N(double, lwork);
            dsytrf_(&uploflag, &n, (double *) amat->data, &lda, ipiv, 
                    work, &lwork, &info);
            break;
        case RATLAS_DCOMPLEX:
            zsytrf_(&uploflag, &n, amat->data, &lda, ipiv, 
                    &cqwork, &lwork, &info);
            lwork = (LAPACK_INT) cqwork.r; 
            cwork = ALLOCA_N(LAPACK_DCOMPLEX, lwork);
            zsytrf_(&uploflag, &n, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv, 
                    cwork, &lwork, &info);
            break;
    }
    
    rowperm = rb_ary_new2(amat->nrow);
    for (i = 0; i < amat->nrow; i++)
    {
        rb_ary_push(rowperm, INT2FIX(ipiv[i]));
    }
    return rb_ary_new3(3, INT2FIX(info), astor, rowperm);
}





/**
 * solves a system of linear equations
 *      A * X = B  or  A' * X = B
 * with a general N-by-N matrix A using the LU factorization computed
 * by getrf.
 * B <- inv( op(A) ) * B, assuming A=LU */
VALUE ratlas_getrs_bang(VALUE self, VALUE trans, VALUE astor, VALUE rowperm,
        VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    struct RArray *pivarr;
    int i, ctrans;
    LAPACK_INT *ipiv, lda, ldb, nrhs, info, n;
    char transflag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;
    n = (LAPACK_INT) amat->ncol;
    ctrans = FIX2INT(trans);
    transflag = ratlas_blasflag2lapack(ctrans);
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    pivarr = RARRAY(rowperm);

    for (i = 0; i < amat->nrow; i++)
    {
        ipiv[i] = (LAPACK_INT) FIX2INT(pivarr->ptr[i]);
    }
    switch (bmat->type) {
        case RATLAS_DFLOAT:
            dgetrs_(&transflag, &n, &nrhs, amat->data,
                    &lda, ipiv, bmat->data, &ldb, &info);
            break;
        case RATLAS_DCOMPLEX:
            zgetrs_(&transflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data,
                    &lda, ipiv, (LAPACK_DCOMPLEX *) bmat->data, &ldb, &info);
            break;
    }
    return rb_ary_new3(2, INT2FIX(info), bstor);
}




/**
 * POTRS solves a system of linear equations A*X = B with a symmetric
 *   positive definite matrix A using the Cholesky factorization
 *   A = U**T*U or A = L*L**T computed by POTRF. */
VALUE ratlas_potrs_bang(VALUE self, VALUE uplo, VALUE astor, VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    LAPACK_INT lda, ldb, nrhs, info, n;
    char uploflag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;
    n = (LAPACK_INT) amat->ncol;
    uploflag = ratlas_blasflag2lapack(FIX2INT(uplo));
    switch (bmat->type) {
        case RATLAS_DFLOAT:
            dpotrs_(&uploflag, &n, &nrhs, amat->data,
                    &lda, bmat->data, &ldb, &info);
            break;
        case RATLAS_DCOMPLEX:
            zpotrs_(&uploflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data,
                    &lda, (LAPACK_DCOMPLEX *) bmat->data, &ldb, &info);
            break;
    }
    return rb_ary_new3(2, INT2FIX(info), bstor);
}




/**
 * SYTRS solves a system of linear equations A*X = B with a real
 * *  symmetric matrix A using the factorization A = U*D*U**T or
 * *  A = L*D*L**T computed by SYTRF */
VALUE ratlas_sytrs_bang(VALUE self, VALUE uplo, VALUE astor, VALUE rowperm,
        VALUE bstor)
{
    Ratlas_Matrix *amat, *bmat;
    struct RArray *pivarr;
    int i;
    LAPACK_INT *ipiv, lda, ldb, nrhs, info, n;
    char uploflag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);

    lda = (LAPACK_INT) amat->nrow;
    ldb = (LAPACK_INT) bmat->nrow;
    nrhs = (LAPACK_INT) bmat->ncol;
    n = (LAPACK_INT) amat->ncol;
    uploflag = ratlas_blasflag2lapack(FIX2INT(uplo));
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    pivarr = RARRAY(rowperm);

    for (i = 0; i < amat->nrow; i++)
    {
        ipiv[i] = (LAPACK_INT) FIX2INT(pivarr->ptr[i]);
    }
    switch (bmat->type) {
        case RATLAS_DFLOAT:
            dsytrs_(&uploflag, &n, &nrhs, amat->data,
                    &lda, ipiv, bmat->data, &ldb, &info);
            break;
        case RATLAS_DCOMPLEX:
            zsytrs_(&uploflag, &n, &nrhs, (LAPACK_DCOMPLEX *) amat->data,
                    &lda, ipiv, (LAPACK_DCOMPLEX *) bmat->data, &ldb, &info);
            break;
    }
    return rb_ary_new3(2, INT2FIX(info), bstor);
}




/**
 * inverse of a matrix using the LU factorization
 * *  computed by getrf.
 * *
 * *  This method inverts U and then computes inv(A) by solving the system
 * *  inv(A)*L = inv(U) for inv(A). */
VALUE ratlas_getri_bang(VALUE self, VALUE astor, VALUE rowperm)
{
    Ratlas_Matrix *amat;
    struct RArray *pivarr;
    int i;
    LAPACK_INT *ipiv, lda, info, n, lwork=-1;
    double qwork, *work;
    LAPACK_DCOMPLEX cqwork, *cwork;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    lda = (LAPACK_INT) amat->nrow;
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    pivarr = RARRAY(rowperm);
    n = (LAPACK_INT) amat->ncol;

    if (amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Can not invert non square matrix.");
    for (i = 0; i < amat->nrow; i++)
    {
        ipiv[i] = (LAPACK_INT) FIX2INT(pivarr->ptr[i]);
    }

    /* Run a query first to find optimal work array by setting lwork=-1 */
    switch (amat->type) {
        case RATLAS_DFLOAT:
            dgetri_(&n, amat->data, &lda, ipiv, &qwork, &lwork,
                    &info);
            lwork = (LAPACK_INT) qwork; 
            work = ALLOCA_N(double, lwork);
            dgetri_(&n, amat->data, &lda, ipiv, work, &lwork,
                    &info);
            break;
        case RATLAS_DCOMPLEX:
            zgetri_(&n, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv, &cqwork, 
                    &lwork, &info);
            lwork = (LAPACK_INT) cqwork.r; 
            cwork = ALLOCA_N(LAPACK_DCOMPLEX, lwork);
            zgetri_(&n, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv, cwork,
                    &lwork, &info);
            break;
    }
    return rb_ary_new3(2, INT2FIX(info), astor);
}




/**
 * POTRI computes the inverse of a real symmetric positive definite
 * *  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
 * *  computed by POTRF. */
VALUE ratlas_potri_bang(VALUE self, VALUE uplo, VALUE astor)
{
    Ratlas_Matrix *amat;
    LAPACK_INT lda, info, n;
    char uplochar;
    int uploint;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    lda = (LAPACK_INT) amat->nrow;
    n = (LAPACK_INT) amat->ncol;
    uploint = FIX2INT(uplo);
    uplochar = ratlas_blasflag2lapack(uploint);

    if (amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Can not invert non square matrix.");

    switch (amat->type) {
        case RATLAS_DFLOAT:
            dpotri_(&uplochar, &n, amat->data, &lda, &info);
            break;
        case RATLAS_DCOMPLEX:
            zpotri_(&uplochar, &n, (LAPACK_DCOMPLEX *) amat->data, &lda, &info);
            break;
    }

#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(uploint, amat);
#endif
    return rb_ary_new3(2, INT2FIX(info), astor);
}




/**
 * SYTRI computes the inverse of a real symmetric indefinite matrix
 * *  A using the factorization A = U*D*U**T or A = L*D*L**T computed by
 * *  SYTRF. */
VALUE ratlas_sytri_bang(VALUE self, VALUE uplo, VALUE astor, VALUE rowperm)
{
    Ratlas_Matrix *amat;
    int i, uploint;
    LAPACK_INT *ipiv, lda, info, n;
    char uplochar;
    struct RArray *pivarr;
    double *work;
    LAPACK_DCOMPLEX *cwork;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    if (amat->nrow != amat->ncol)
        rb_raise(rb_eArgError, "Can not invert non square matrix.");
    
    lda = (LAPACK_INT) amat->nrow;
    ipiv = ALLOCA_N(LAPACK_INT, amat->nrow);
    pivarr = RARRAY(rowperm);
    n = (LAPACK_INT) amat->ncol;
    uploint = FIX2INT(uplo);
    uplochar = ratlas_blasflag2lapack(uploint);

    for (i = 0; i < amat->nrow; i++)
    {
        ipiv[i] = (LAPACK_INT) FIX2INT(pivarr->ptr[i]);
    }

    switch (amat->type) {
        case RATLAS_DFLOAT:
            work = ALLOCA_N(double, n);
            dsytri_(&uplochar, &n, amat->data, &lda, ipiv, work, &info);
            break;
        case RATLAS_DCOMPLEX:
            cwork = ALLOCA_N(LAPACK_DCOMPLEX, n);
            zsytri_(&uplochar, &n, (LAPACK_DCOMPLEX *) amat->data, &lda, ipiv,
                    cwork, &info);
            break;
    }
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(uploint, amat);
#endif
    return rb_ary_new3(2, INT2FIX(info), astor);
}




/**
 * computes the singular value decomposition (SVD) of a real M-by-N
 * matrix A, optionally computing the left and/or right singular
 * vectors. The SVD is written 
 *
 *         A = U * SIGMA * transpose(V) 
 *
 * where SIGMA is an M-by-N matrix which is zero except for its
 * min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and V
 * is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA are
 * the singular values of A; they are real and non-negative, and are
 * returned in descending order.  The first min(m,n) columns of U and
 * V are the left and right singular vectors of A. */
VALUE ratlas_gesvd_bang(VALUE self, VALUE astor)
{
    Ratlas_Matrix *amat, *smat, *umat, *vtmat;
    int sdim;
    LAPACK_INT n, m, lda, ldu, ldvt, info, lwork=-1;
    char jobu, jobvt;
    double qwork, *work, *zeros, *rwork;
    LAPACK_DCOMPLEX cqwork, *cwork, *czeros;
    VALUE sstor=Qnil, ustor=Qnil, vtstor=Qnil;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    
    m = (LAPACK_INT)amat->nrow;
    n = (LAPACK_INT)amat->ncol;
    
    lda = (LAPACK_INT)m;
    ldu = (LAPACK_INT)m;
    ldvt = (LAPACK_INT)n;
    
    jobu = 'A';
    jobvt = 'A';
    
    sdim = ratlas_imin(m, n);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            zeros = xcalloc(sdim, sizeof(double));
            smat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE,
                           sdim, 1, zeros);
            free(zeros);
            zeros = xcalloc(ldu*m, sizeof(double));
            umat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE,
                           ldu, m, zeros);
            free(zeros);
            zeros = xcalloc(n*n, sizeof(double));
            vtmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE,
                            n, n, zeros);
            free(zeros);
            dgesvd_(&jobu, &jobvt, &m, &n, (double *)amat->data, &lda, 
                (double *)smat->data, (double *)umat->data, &ldu,
                vtmat->data, &ldvt, &qwork, &lwork, &info);
            lwork = (LAPACK_INT) qwork;
            work = ALLOCA_N(double, lwork);
            dgesvd_(&jobu, &jobvt, &m, &n, (double *)amat->data, &lda, 
                (double *)smat->data, (double *)umat->data, &ldu, 
                (double *)vtmat->data, &ldvt, work, &lwork, &info);
            sstor = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, smat);
            ustor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, umat);
            vtstor = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, vtmat);
            break;
        case RATLAS_DCOMPLEX:
            czeros = xcalloc(sdim, sizeof(LAPACK_DCOMPLEX));
            smat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE,
                           sdim, 1, czeros);
            free(czeros);
            czeros = xcalloc(ldu*m, sizeof(LAPACK_DCOMPLEX));
            umat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE,
                           ldu, m, czeros);
            free(czeros);
            czeros = xcalloc(n*n, sizeof(LAPACK_DCOMPLEX));
            vtmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE,
                            n, n, czeros);
            free(czeros);
            rwork = xcalloc(5*ratlas_imin(m,n), sizeof(double));
            zgesvd_(&jobu, &jobvt, &m, &n, (LAPACK_DCOMPLEX *)amat->data, &lda, 
                (double *)smat->data, (LAPACK_DCOMPLEX *)umat->data,
                &ldu, vtmat->data, &ldvt, &cqwork, &lwork, rwork, &info);
            lwork = (LAPACK_INT) cqwork.r;
            cwork = ALLOCA_N(LAPACK_DCOMPLEX, lwork);
            zgesvd_(&jobu, &jobvt, &m, &n, (LAPACK_DCOMPLEX *)amat->data, &lda, 
                (double *)smat->data, (LAPACK_DCOMPLEX *)umat->data,
                &ldu, (LAPACK_DCOMPLEX *)vtmat->data, 
                &ldvt, cwork, &lwork, rwork, &info);
            free(rwork);
            sstor = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, smat);
            ustor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, umat);
            vtstor = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, vtmat);
            break;
    }
    
    return rb_ary_new3(5, INT2FIX(info), astor, sstor, ustor, vtstor);
}

