/*
 *  ratlas_cblas.c
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#include "ratlas_cblas.h"




/* ATLAS Level 1 BLAS */


/**
 * scales a vector by a constant
 *     x <- alpha*x */
VALUE ratlas_scal_bang(VALUE self, VALUE alpha, VALUE xstorage)
{
    Ratlas_Matrix *xmat;
    int n;
    double calpha[2];

    Data_Get_Struct(xstorage, Ratlas_Matrix, xmat);
    n = ratlas_imax(xmat->ncol, xmat->nrow);
    if (n != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (xmat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            calpha[0] = NUM2DBL(alpha);
            cblas_dscal(n, calpha[0], xmat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, "Expected complex.");
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zscal(n, calpha, xmat->data, 1);
            break;
    }

    return xstorage;
}




/**
 * constant times a vector plus a vector
 *     y <- alpha*x + y */
VALUE ratlas_axpy_bang(VALUE self, VALUE alpha, VALUE xstorage,  VALUE ystorage)
{
    Ratlas_Matrix *xmat, *ymat;
    int n, ny;
    double calpha[2];


    Data_Get_Struct(xstorage, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystorage, Ratlas_Matrix, ymat);
    n = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    
    if (n != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow)
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != n)
        rb_raise(rb_eArgError, "Vectors must have same size.");
            
    switch (xmat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            calpha[0] = NUM2DBL(alpha);
            cblas_daxpy(n, calpha[0], xmat->data, 1, ymat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            if (ymat->type != RATLAS_DCOMPLEX)
                rb_raise(rb_eArgError, 
                        "Expected complex vector y.");
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, 
                        "Expected complex alpha.");
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zaxpy(n, calpha, xmat->data, 1, ymat->data, 1);
            break;
    }

    return ystorage;    
}




/**
 * inner product <x,y> for real and complex vectors */
VALUE ratlas_dot(VALUE self, VALUE xstor, VALUE ystor)
{
    Ratlas_Matrix *xmat, *ymat;
    int nx, ny;
    double retval[2];
    VALUE retobj=Qnil;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);

    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow)
        rb_raise(rb_eArgError, "Expected vector.");
    if (nx != ny)
        rb_raise(rb_eArgError, "Vectors must have same size.");

    
    switch (xmat->type) {
        case RATLAS_DFLOAT:
            retval[0] = (double) cblas_ddot(nx, xmat->data, 1, 
                    ymat->data, 1);
            retobj =  rb_float_new(retval[0]);
            break;
        case RATLAS_DCOMPLEX:
            if ( (xmat->type != RATLAS_DCOMPLEX) || 
                (ymat->type != RATLAS_DCOMPLEX) )
                rb_raise(rb_eArgError, 
                        "Expected complex vector y.");
            cblas_zdotu_sub(nx, xmat->data, 1, ymat->data, 1, 
                    retval);
            retobj = ratlas_rb_complex_new(retval[0], retval[1]);
            break;
    }
    return retobj;
}




/**
 * inner product <x,x> for complex vector.
 * */
VALUE ratlas_dotc(VALUE self, VALUE xstor, VALUE ystor)
{
    Ratlas_Matrix *xmat, *ymat;
    int nx, ny;
    double retval[2];
    VALUE retobj;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);

    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow)
        rb_raise(rb_eArgError, "Expected vector.");
    if (nx != ny)
        rb_raise(rb_eArgError, "Vectors must have same size.");

    if ( (xmat->type != RATLAS_DCOMPLEX) || 
            (ymat->type != RATLAS_DCOMPLEX))
        rb_raise(rb_eArgError, "Expected complex vectors.");
    cblas_zdotc_sub(nx, xmat->data, 1, ymat->data, 1, retval);
    retobj = ratlas_rb_complex_new(retval[0], retval[1]);

    return retobj;
}




/**
 * euclidean norm of a vector
 *     sqrt( x'*x ) */
VALUE ratlas_nrm2(VALUE self, VALUE xstor)
{
    Ratlas_Matrix *xmat;
    int n;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);

    n = ratlas_imax(xmat->ncol, xmat->nrow);
    if (n != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (xmat->type) {
        case RATLAS_DFLOAT:
            return rb_float_new(cblas_dnrm2(n, xmat->data, 1));
        case RATLAS_DCOMPLEX:
            return rb_float_new(cblas_dznrm2(n, xmat->data, 1));
        default:
            rb_raise(rb_eArgError, "Unknown type in vector.");
    }
}




/**
 * takes the sum of the absolute values */
VALUE ratlas_asum(VALUE self, VALUE xstor)
{
    Ratlas_Matrix *xmat;
    int n;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);

    n = ratlas_imax(xmat->ncol, xmat->nrow);
    if (n != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (xmat->type) {
        case RATLAS_DFLOAT:
            return rb_float_new(cblas_dasum(n, xmat->data, 1));
        case RATLAS_DCOMPLEX:
            return rb_float_new(cblas_dzasum(n, xmat->data, 1));
        default:
            rb_raise(rb_eArgError, "Unknown type in vector.");
    }
}




/**
 * finds the index of element having max. absolute value */
VALUE ratlas_iamax(VALUE self, VALUE xstor)
{
    Ratlas_Matrix *xmat;
    int n;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);

    n = ratlas_imax(xmat->ncol, xmat->nrow);
    if (n != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (xmat->type) {
        case RATLAS_DFLOAT:
            return INT2NUM(cblas_idamax(n, xmat->data, 1));
        case RATLAS_DCOMPLEX:
            return INT2NUM(cblas_izamax(n, xmat->data, 1));
        default:
            rb_raise(rb_eArgError, "Unknown type in vector.");
    }
}




/* ATLAS Level 2 BLAS */

/**
 * matrix-vector operations
 *      y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y */
VALUE ratlas_gemv_bang(VALUE self, VALUE transa, VALUE alpha, VALUE astor,
        VALUE xstor, VALUE beta, VALUE ystor)
{
    Ratlas_Matrix *amat, *xmat, *ymat;
    int nx, ny, ctransa, lda;
    double calpha[2];
    double cbeta[2];
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);
    ctransa = FIX2INT(transa);
    
    /* Checking arguments */
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    if (ny != ymat->ncol * ymat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ctransa == CblasNoTrans)
    {
        if (amat->ncol != xmat->nrow)
            rb_raise(rb_eArgError, "ncol A must be equal nrow x.");
        if (amat->nrow != ymat->nrow) 
            rb_raise(rb_eArgError, "Return vector too small.");
    } else
    {
        if (amat->nrow != xmat->nrow)
            rb_raise(rb_eArgError, 
                "ncol transpose(A) must be equal nrow x.");
        if (amat->ncol != ymat->nrow) 
            rb_raise(rb_eArgError, "Return vector too small.");

    }
    
    lda = amat->nrow;
    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dgemv(CblasColMajor, ctransa, amat->nrow,
                amat->ncol, NUM2DBL(alpha), amat->data, 
                lda, xmat->data, 1, NUM2DBL(beta),
                ymat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zgemv(CblasColMajor, ctransa, amat->nrow,
                amat->ncol, calpha, amat->data, 
                lda, xmat->data, 1, cbeta,
                ymat->data, 1);
            break;
    }
    return ystor;
}




/**
 * matrix-vector  operation
 *      y := alpha*A*x + beta*y,
 * where alpha and beta are scalars, x and y are n element vectors and
 * A is an n by n symmetric matrix. */
VALUE ratlas_symv_bang(VALUE self, VALUE alpha, VALUE astor, VALUE xstor,
        VALUE beta, VALUE ystor)
{
    Ratlas_Matrix *amat, *xmat, *ymat;
    int nx, ny, lda;
    double calpha[2];
    double cbeta[2];
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);
    
    /* Checking arguments */
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    if (ny != ymat->ncol * ymat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (nx != ny) rb_raise(rb_eArgError, "Vectors must have same length.");
    
    lda = amat->nrow;
    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dsymv(CblasColMajor, CblasLower, amat->nrow,
                NUM2DBL(alpha), amat->data, 
                lda, xmat->data, 1, NUM2DBL(beta),
                ymat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zhemv(CblasColMajor, CblasLower, amat->nrow,
                calpha, amat->data, 
                lda, xmat->data, 1, cbeta,
                ymat->data, 1);
            break;
    }
    return ystor;
}




/**
 * matrix-vector operations
 *      x := A*x,   or   x := A'*x,
 * where x is an n element vector and  A is an n by n unit, or non-unit,
 * upper or lower triangular matrix. */
VALUE ratlas_trmv_bang(VALUE self, VALUE astor, VALUE xstor, VALUE uplo,
        VALUE trans, VALUE diag)
{
    Ratlas_Matrix *amat, *xmat;
    int n, nx, lda;
    int cuplo, ctrans, cdiag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    cuplo = FIX2INT(uplo);
    ctrans = FIX2INT(trans);
    cdiag = FIX2INT(diag);
    n = ratlas_imax(amat->ncol, amat->nrow);
    lda = amat->nrow;
    
    /* Checking arguments */
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dtrmv(CblasColMajor, cuplo, ctrans, cdiag,
                    n, amat->data, lda, xmat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            cblas_ztrmv(CblasColMajor, cuplo, ctrans, cdiag,
                    n, amat->data, lda, xmat->data, 1);
            break;
    }
    return xstor;
}




/**
 * solves one of the systems of equations
 *      A*x = b,   or   A'*x = b 
 * where b and x are n element vectors and A is an n by n unit, or
 * non-unit, upper or lower triangular matrix. */
VALUE ratlas_trsv_bang(VALUE self, VALUE astor, VALUE xstor, VALUE uplo,
        VALUE trans, VALUE diag)
{
    Ratlas_Matrix *amat, *xmat;
    int n, nx, lda;
    int cuplo, ctrans, cdiag;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    cuplo = FIX2INT(uplo);
    ctrans = FIX2INT(trans);
    cdiag = FIX2INT(diag);
    n = ratlas_imax(amat->ncol, amat->nrow);
    lda = amat->nrow;
    
    /* Checking arguments */
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");

    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dtrsv(CblasColMajor, cuplo, ctrans, cdiag,
                    n, amat->data, lda, xmat->data, 1);
            break;
        case RATLAS_DCOMPLEX:
            cblas_ztrsv(CblasColMajor, cuplo, ctrans, cdiag,
                    n, amat->data, lda, xmat->data, 1);
            break;
    }
    return xstor;
}




/**
 *  rank 1 operation
 *       A := alpha*x*y' + A */
VALUE ratlas_ger_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor)
{
    Ratlas_Matrix *xmat, *ymat, *amat;
    int nx, ny, lda;
    double calpha[2];
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    lda = amat->nrow;

    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if ( (amat->nrow != nx) || (amat->ncol != ny) )
        rb_raise(rb_eArgError, "Matrix a is too small to hold x*y'.");
    if ( (amat->type != xmat->type) || (amat->type != ymat->type) )
        rb_raise(rb_eArgError, 
                "Matrix and vectors must have same type.");
    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dger(CblasColMajor, nx, ny, NUM2DBL(alpha),
                    xmat->data, 1, ymat->data, 1, 
                    amat->data, lda);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zgeru(CblasColMajor, nx, ny, calpha,
                    xmat->data, 1, ymat->data, 1, 
                    amat->data, lda);
            break;

    }
    return astor;
}




/**
 * rank 1 operation
 *      A := alpha*x*conjg( y' ) + A */
VALUE ratlas_gerc_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor)
{
    Ratlas_Matrix *xmat, *ymat, *amat;
    int nx, ny, lda;
    double calpha[2];
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    lda = amat->nrow;

    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if ( (amat->nrow != nx) || (amat->ncol != ny) )
        rb_raise(rb_eArgError, "Matrix a is too small to hold x*y'.");
    if ( (amat->type != xmat->type) || (amat->type != ymat->type) )
        rb_raise(rb_eArgError, 
                "Matrix and vectors must have same type.");
    if (amat->type != RATLAS_DCOMPLEX)
        rb_raise(rb_eArgError, "Expected complex inputs.");
    
    calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
    calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
    cblas_zgerc(CblasColMajor, nx, ny, calpha,
            xmat->data, 1, ymat->data, 1, 
            amat->data, lda);
    return astor;
}




/**
 * symmetric rank 2 operation
 *      A := alpha*x*y' + alpha*y*x' + A,
 *   where alpha is a scalar, x and y are n element vectors and A is an n
 *   by n symmetric matrix */
VALUE ratlas_syr_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE astor)
{
    Ratlas_Matrix *xmat, *amat;
    int nx, lda;
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    lda = amat->nrow;

    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if ( (amat->nrow != nx) )
        rb_raise(rb_eArgError, "Matrix a is too small to hold x*x'.");
    if ( (amat->type != xmat->type) )
        rb_raise(rb_eArgError, 
                "Matrix and vector must have same type.");
    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dsyr(CblasColMajor, CblasLower, nx,
                    NUM2DBL(alpha), xmat->data, 1, 
                    amat->data, lda);
            break;
        case RATLAS_DCOMPLEX:
            cblas_zher(CblasColMajor, CblasLower, nx,
                    NUM2DBL(alpha), xmat->data, 1, 
                    amat->data, lda);
            break;

    }
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, amat);
#endif
    return astor;
}




VALUE ratlas_syr2_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor)
{
    Ratlas_Matrix *xmat, *ymat, *amat;
    int nx, ny, lda;
    double calpha[2];
    
    Data_Get_Struct(xstor, Ratlas_Matrix, xmat);
    Data_Get_Struct(ystor, Ratlas_Matrix, ymat);
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    nx = ratlas_imax(xmat->ncol, xmat->nrow);
    ny = ratlas_imax(ymat->ncol, ymat->nrow);
    lda = amat->nrow;

    if (nx != xmat->ncol * xmat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != ymat->ncol * ymat->nrow) 
        rb_raise(rb_eArgError, "Expected vector.");
    if (ny != nx)
        rb_raise(rb_eArgError, "Vectors must have same length.");
    if ( (amat->nrow != nx) || (amat->ncol != ny) )
        rb_raise(rb_eArgError, "Matrix a is too small to hold x*y'.");
    if ( (amat->type != xmat->type) || (amat->type != ymat->type) )
        rb_raise(rb_eArgError, 
                "Matrix and vectors must have same type.");
    switch (amat->type) {
        case RATLAS_DFLOAT:
            cblas_dsyr2(CblasColMajor, CblasLower, nx, 
                    NUM2DBL(alpha),
                    xmat->data, 1, ymat->data, 1, 
                    amat->data, lda);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zher2(CblasColMajor, CblasLower, nx, 
                    calpha,
                    xmat->data, 1, ymat->data, 1, 
                    amat->data, lda);
            break;

    }

#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, amat);
#endif
    return astor;
}




/* ATLAS Level 3 BLAS */

/**
 * matrix-matrix operations
 * 
 *      C := alpha*op( A )*op( B ) + beta*C */
VALUE ratlas_gemm_bang(VALUE self, VALUE transa, VALUE transb, VALUE alpha,
        VALUE astor, VALUE bstor, VALUE beta, VALUE cstor)
{
    int lda, ldb, ldc, ctransa, ctransb, n, k, m;
    Ratlas_Matrix *amat, *bmat, *cmat;
    double calpha[2], cbeta[2];

    /* alpha and beta are scalars, and A, B and C are matrices, with 
     * op( A ) *  an m by k matrix,  op( B )  a  k by n matrix and  
     * C an m by n matrix.
     */

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

    lda = amat->nrow;
    ldb = bmat->nrow;
    ldc = cmat->nrow;
    
    ctransa = FIX2INT(transa);
    ctransb = FIX2INT(transb);

    if (ctransa > CblasNoTrans)
    {
        m = amat->ncol;
        k = amat->nrow;
    } else
    {
        m = amat->nrow;
        k = amat->ncol;
    }

    if (ctransb > CblasNoTrans)
    {
        if (bmat->ncol != k)
            rb_raise(rb_eArgError, 
                    "Dimension mismatch for op(A)*B'.");
        n = bmat->nrow;

    } else
    {
        if (bmat->nrow != k)
            rb_raise(rb_eArgError, 
                    "Dimension mismatch for op(A)*B.");
        n = bmat->ncol;
    }
    if ((cmat->nrow != m) || (cmat->ncol != n))
            rb_raise(rb_eArgError, "Dimension mismatch for C.");

    switch (cmat->type) {
        case RATLAS_DFLOAT:
            cblas_dgemm(CblasColMajor, ctransa, ctransb,
                    m, n, k, NUM2DBL(alpha),
                    amat->data, lda, 
                    bmat->data, ldb,
                    NUM2DBL(beta), cmat->data, ldc);
            break;    
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zgemm(CblasColMajor, ctransa, ctransb,
                    m, n, k, calpha, amat->data,
                    lda, bmat->data, ldb,
                    cbeta, cmat->data, ldc);
            break;
    }
    
    return cstor;
}




/**
 * matrix-matrix operations
 *      C := alpha*A*B + beta*C,
 *   or
 *      C := alpha*B*A + beta*C,
 *   A is a symmetric matrix and  B and C are  m by n matrices */
VALUE ratlas_symm_bang(VALUE self, VALUE side, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *bmat, *cmat;
    double calpha[2], cbeta[2];
    int cside, n, m, lda, ldb, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

    if ( (amat->ncol != bmat->nrow)    || (amat->nrow != bmat->ncol) )
        rb_raise(rb_eArgError, "Dimension mismatch.");
    
    lda = amat->nrow;
    ldb = bmat->nrow;
    ldc = cmat->nrow;
    m = cmat->nrow;
    n = cmat->ncol;
    
    /* Left means C := alpha*A*B + beta*C
     * Right means C := alpha*B*A + beta*C */
    cside = FIX2INT(side);

    switch (cmat->type) {
        case RATLAS_DFLOAT:
            cblas_dsymm(CblasColMajor, cside, CblasLower, m,
                    n, NUM2DBL(alpha), amat->data,
                    lda, bmat->data, ldb, NUM2DBL(beta),
                    cmat->data, ldc);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zsymm(CblasColMajor, cside, CblasLower, m,
                    n, calpha, amat->data,
                    lda, bmat->data, ldb, cbeta,
                    cmat->data, ldc);
            break;
    }
    return cstor;
}




/**
 * complex matrix-matrix operations
 *      C := alpha*A*B + beta*C,
 *   or
 *      C := alpha*B*A + beta*C
 * A is an hermitian matrix and  B and C are m by n matrices */
VALUE ratlas_hemm_bang(VALUE self, VALUE side, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *bmat, *cmat;
    double calpha[2], cbeta[2];
    int cside, n, m, lda, ldb, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

    if ( (amat->ncol != bmat->nrow)    || (amat->nrow != bmat->ncol) )
        rb_raise(rb_eArgError, "Dimension mismatch.");
    
    lda = amat->nrow;
    ldb = bmat->nrow;
    ldc = cmat->nrow;
    m = cmat->nrow;
    n = cmat->ncol;
    
    /* Left means C := alpha*A*B + beta*C
     * Right means C := alpha*B*A + beta*C */
    cside = FIX2INT(side);

    calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
    calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
    cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
    cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
    cblas_zhemm(CblasColMajor, cside, CblasLower, m, n, calpha, amat->data,
            lda, bmat->data, ldb, cbeta, cmat->data, ldc);
    return cstor;
}




/**
 * symmetric rank k operations
 *      C := alpha*A*A' + beta*C,
 *   or
 *      C := alpha*A'*A + beta*C
 * only lower triangular part of C is referenced. */
VALUE ratlas_syrk_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat;
    double calpha[2], cbeta[2];
    int ctrans, n, k, lda, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

     /* C := alpha*A*A' + beta*C,
     * or
     * C := alpha*A'*A + beta*C,
     * where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     * and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     * in the second case. */
        
    ctrans = FIX2INT(trans);
    if (ctrans == CblasNoTrans) 
    {
        n = amat->nrow;
        k = amat->ncol;
    } else
    {
        n = amat->ncol;
        k = amat->nrow;
    }
    if ( (cmat->nrow != n) || (cmat->ncol != n) )
        rb_raise(rb_eArgError, "C is not nxn");

    lda = amat->nrow;
    ldc = cmat->nrow;

    switch (cmat->type) {
        case RATLAS_DFLOAT:
            cblas_dsyrk(CblasColMajor, CblasLower, ctrans,
                    n, k, NUM2DBL(alpha), amat->data,
                    lda, NUM2DBL(beta),
                    cmat->data, ldc);

            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zsyrk(CblasColMajor, CblasLower, ctrans,
                    n, k, calpha, amat->data,
                    lda, cbeta,
                    cmat->data, ldc);
            break;
    }
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, cmat);
#endif
    return cstor;
}




/**
 * hermitian rank k operations
 *  
 *      C := alpha*A*conjg( A' ) + beta*C,
 *  
 *   or
 *       C := alpha*conjg( A' )*A + beta*C
 * only lower triangular part of C is referenced. */
VALUE ratlas_herk_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat;
    int ctrans, n, k, lda, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

     /* C := alpha*A*AH + beta*C,
     * or
     * C := alpha*AH*A + beta*C,
     * where  alpha and beta  are scalars, C is an  n by n hermitian matrix
     * and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     * in the second case. */
        
    ctrans = FIX2INT(trans);
    if (ctrans == CblasNoTrans) 
    {
        n = amat->nrow;
        k = amat->ncol;
    } else
    {
        n = amat->ncol;
        k = amat->nrow;
    }
        
    if ( (cmat->nrow != n) || (cmat->ncol != n) )
        rb_raise(rb_eArgError, "C is not nxn");

    lda = amat->nrow;
    ldc = cmat->nrow;

    cblas_zherk(CblasColMajor, CblasLower, ctrans,
                n, k, NUM2DBL(alpha), amat->data,
                lda, NUM2DBL(beta),
                cmat->data, ldc);
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, cmat);
#endif
    return cstor;
}




/**
 * symmetric rank 2k operations
 *      C := alpha*A*B' + alpha*B*A' + beta*C,
 *   or
 *      C := alpha*A'*B + alpha*B'*A + beta*C 
 * only lower triangular part of C is referenced. */
VALUE ratlas_syr2k_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *bmat, *cmat;
    double calpha[2], cbeta[2];
    int ctrans, n, k, lda, ldb, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

    ctrans = FIX2INT(trans);
    lda = amat->nrow;
    ldb = bmat->nrow;
    ldc = cmat->nrow;

    if (ctrans == CblasNoTrans) 
    {
        n = amat->nrow;
        k = amat->ncol;
    } else
    {
        n = amat->ncol;
        k = amat->nrow;
    }
    if ( (amat->nrow != bmat->nrow) || (amat->ncol != bmat->ncol) )
        rb_raise(rb_eArgError, "Matrices must have equal dimension.");
    if ( (cmat->nrow != n) || (cmat->ncol != n) )
        rb_raise(rb_eArgError, "C is not nxn");    

    switch (cmat->type) {
        case RATLAS_DFLOAT:
            cblas_dsyr2k(CblasColMajor, CblasLower, ctrans,
                    n, k, NUM2DBL(alpha), amat->data,
                    lda, bmat->data, ldb, NUM2DBL(beta),
                    cmat->data, ldc);
            break;
        case RATLAS_DCOMPLEX:
            calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cbeta[0] = NUM2DBL(rb_ivar_get(beta, id_atre));
            cbeta[1] = NUM2DBL(rb_ivar_get(beta, id_atim));
            cblas_zsyr2k(CblasColMajor, CblasLower, ctrans,
                    n, k, calpha, amat->data,
                    lda, bmat->data, ldb, cbeta,
                    cmat->data, ldc);
            break;
    }
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, cmat);
#endif
    return cstor;
}



/**
 * hermitian rank 2k operations
 *      C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
 *   or
 *      C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C
 * only lower triangular part of C is referenced. */
VALUE ratlas_her2k_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor)
{
    Ratlas_Matrix *amat, *bmat, *cmat;
    double calpha[2];
    int ctrans, n, k, lda, ldb, ldc;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(bstor, Ratlas_Matrix, bmat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);

    lda = amat->nrow;
    ldb = bmat->nrow;
    ldc = cmat->nrow;
    ctrans = FIX2INT(trans);
    if (ctrans == CblasNoTrans) 
    {
        n = amat->nrow;
        k = amat->ncol;
    } else
    {
        n = amat->ncol;
        k = amat->nrow;
    }
        
    if ( (cmat->nrow != n) || (cmat->ncol != n) )
        rb_raise(rb_eArgError, "C is not nxn");    

    calpha[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
    calpha[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
    cblas_zher2k(CblasColMajor, CblasLower, ctrans,
            n, k, calpha, amat->data,
            lda, bmat->data, ldb, NUM2DBL(beta),
            cmat->data, ldc);
#ifdef PROPAGATESYM
    ratlas_matrix_sym2dense(CblasLower, cmat);
#endif
    return cstor;
}
