/*
 *  ratlas_func.c
 *  Ruby Numeric module 
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#include "ratlas_func.h"
#define DEBUG 1



const size_t ratlas_sizeof[RATLAS_NTYPES] = {
    sizeof(double), // FLOAT
    sizeof(Ratlas_Complex)  // COMPLEX
    };




char ratlas_blasflag2lapack(int cblasflag)
{
    switch (cblasflag) {
        case CblasNoTrans:
            return 'N';
        case CblasTrans:
            return 'T';
        case CblasConjTrans:
            return 'C';
        case CblasUpper:
            return 'U';
        case CblasLower:
            return 'L';
        case CblasNonUnit:
            return 'N';
        case CblasUnit:
            return 'U';
        case CblasLeft:
            return 'L';
        case CblasRight:
            return 'R';
        default:
            return 'N';
    }
}




int ratlas_imax(int a, int b)
{
    return (a>=b) ? a : b;
}




int ratlas_imin(int a, int b)
{
    return (a<=b) ? a : b;
}




int ratlas_swap_major(int dim1, int dim2, long int idx)
{
    return (idx%dim2)*dim1 + idx/dim2; 
}



/*
 * Equivalent to Complex::new(re, im) 
 */
VALUE ratlas_rb_complex_new(double re, double im)
{
    VALUE argv[2];

    argv[0] = rb_float_new(re);
    argv[1] = rb_float_new(im);
    return rb_class_new_instance(2, argv, complex_class);
}




VALUE ratlas_get_complex(Ratlas_Complex *c, int offset)
{
    return ratlas_rb_complex_new(c[offset].r, c[offset].i);
}




static VALUE ratlas_get_float(double *dataptr, int offset)
{
    return rb_float_new(dataptr[offset]);
}




void ratlas_set_complex(Ratlas_Complex *c, int offset, VALUE newval)
{
    c[offset].r = NUM2DBL(rb_ivar_get(newval, id_atre));
    c[offset].i = NUM2DBL(rb_ivar_get(newval, id_atim));
}




static void ratlas_set_float(double *dataptr, int offset, VALUE newval)
{
    dataptr[offset] = NUM2DBL(newval);
}


/**
 * Store a ruby row array into a column major matrix */
void ratlas_store_complex_row(int argc, VALUE *argv, Ratlas_Matrix *rmat, int i)
{
    int j, m;
    double re, im;
    Ratlas_Complex *c;

    m = rmat->ncol;
    c = (Ratlas_Complex *) rmat->data;

    for (j = 0; j < argc; j++)
    {
        re = rb_ivar_get(argv[j], id_atre);
        im = rb_ivar_get(argv[j], id_atim);
        c[i+j*m].r = NUM2DBL(re);
        c[i+j*m].i = NUM2DBL(im);
    }
}




void ratlas_store_complex_column(int argc, VALUE *argv, Ratlas_Matrix *rmat,
        int j)
{
    int i, offset;
    double re, im;
    Ratlas_Complex *c;

    offset = j * rmat->ncol;
    c = (Ratlas_Complex *) rmat->data;

    for (i = 0; i < argc; i++)
    {
        re = rb_ivar_get(argv[i], id_atre);
        im = rb_ivar_get(argv[i], id_atim);
        c[i+offset].r = NUM2DBL(re);
        c[i+offset].i = NUM2DBL(im);
    }
}




/**
 * Store a ruby row array into a column major matrix */
void ratlas_store_float_row(int argc, VALUE *argv, Ratlas_Matrix *rmat, int i)
{
    int j, m;
    double *r;
    
    m = rmat->nrow;
    r = (double *) rmat->data;

    for (j = 0; j < argc; j++)
    {
        ratlas_assertnum(*argv);
        r[i+j*m] = NUM2DBL(*(argv++));
    }
}



void ratlas_store_float_column(int argc, VALUE *argv, Ratlas_Matrix *rmat,
                     int j)
{
    int i, offset;
    double *r;

    offset = j * rmat->ncol;
    r = (double *) rmat->data;
    for (i = 0; i < argc; i++)
    {
        ratlas_assertnum(*argv);
        r[offset+i] = NUM2DBL(*(argv++));
    }
}




Ratlas_Matrix* ratlas_matrix_alloc (int type, int matrixtype, int nrow, 
        int ncol, void *dataptr)
{
    unsigned long int nelem;
    size_t elemsize;
    Ratlas_Matrix *rmat;
    
    switch (matrixtype) {
        case RATLAS_DENSE:
            nelem = nrow*ncol;
            break;
        case RATLAS_SYMPACK:
            if (nrow != ncol)
                rb_raise(rb_eArgError, "dim1 != dim2 for symmetric matrix.");
            nelem = (nrow*(nrow+1))/2;
            break;
        default:
            nelem = nrow*ncol;
            break;
    }
    
    rmat = ALLOC(Ratlas_Matrix);
    rmat->nrow = nrow;
    rmat->ncol = ncol;
    rmat->type = type;
    rmat->matrixtype = matrixtype;

    elemsize = ratlas_sizeof[type];
    if (nrow<1 || ncol<1) 
    {
           rmat->data = NULL;
        rmat->memsize = 0;
        return rmat;
    }
    

    rmat->memsize = nelem*elemsize;
/*    rmat->data = xcalloc(nelem, elemsize); */
    rmat->data = xmalloc(rmat->memsize); 

    if (dataptr != NULL)
        memcpy((void *)rmat->data, (void *)dataptr, rmat->memsize);

    return rmat;
}




Ratlas_Matrix* ratlas_tripack_alloc (int type, int matrixtype, int nrow, 
        int ncol, int upper, void *dataptr)
{
    unsigned long int nelem;
    size_t elemsize;
    Ratlas_Matrix *rmat;
    
    if (matrixtype != RATLAS_TRIPACK)
        rb_raise(rb_eArgError, "Matrix is not triangular packed.");
    if (upper)
    {
        nelem = (ncol*(ncol+1))/2;
        if ( ncol > nrow ) nelem -= ((ncol-nrow) * (ncol-nrow +1))/2;
    } else
    {
        nelem = (nrow*(nrow+1))/2;
        if ( nrow > ncol ) nelem -= ((nrow-ncol) * (nrow-ncol +1))/2;
    }
    
    rmat = ALLOC(Ratlas_Matrix);
    rmat->nrow = nrow;
    rmat->ncol = ncol;
    rmat->type = type;
    rmat->matrixtype = matrixtype;

    elemsize = ratlas_sizeof[type];
    if (nrow<1 || ncol<1) 
    {
        rmat->data = NULL;
        rmat->memsize = 0;
        return rmat;
    }
    

    rmat->memsize = nelem*elemsize;
    rmat->data = xcalloc(nelem, elemsize);

       if (dataptr != NULL)
    {
        memcpy((void *)rmat->data, (void *)dataptr, rmat->memsize);
    }

    return rmat;
}




void ratlas_matrix_free(Ratlas_Matrix *rmat)
{
    if (!rmat)
    {
#ifdef DEBUG
    fprintf(stderr, "DEBUG: free of NULL\n");
#endif
        abort();
    }
    xfree(rmat->data);
    xfree(rmat);
}




void ratlas_matrix_ptr_free(Ratlas_Matrix *rmat)
{
    xfree(rmat);
}




/**
  * Check if Fixnum or Float.  Return 0 if not.  */
int ratlas_isnumeric(VALUE obj)
{
    int retval;

    switch (TYPE(obj)) {
        case T_FIXNUM:
            retval = 1;
            break;
        case T_FLOAT:
            retval = 1;
            break;
        default:
            retval = 0;
            break;
    }
    return retval;
}




void ratlas_assertnum(VALUE obj)
{
    if(!ratlas_isnumeric(obj))
        rb_raise(rb_eTypeError, "Expect Numeric type.");
}




/*
 * Check for Float or Complex.  Return -1 if unknown type.
 */
int ratlas_get_type(VALUE obj)
{
    if (TYPE(obj) == T_FLOAT)
        return RATLAS_DFLOAT;
    else if (rb_respond_to(obj, id_imag))
        return RATLAS_DCOMPLEX;
    else
        return -1;
}




int ratlas_check_type(VALUE obj, int ratlastype)
{
    if (ratlas_get_type(obj) == ratlastype)
        return 1;
    return 0;
}




/*
 * Create new matrix from ruby array. Ruby array must be rank 2.
 */
VALUE ratlas_matrix_new(VALUE self, VALUE arg)
{
/*     struct RArray *array1, *array2; */
    VALUE storage, array2;
    Ratlas_Matrix *rmat;
    int size1, size2, i;

    Check_Type(arg, T_ARRAY);
/*     array1 = RARRAY(arg); */
/*     size1 = array1->len; */
    size1 = RARRAY_LEN(arg);
    if (size1<1) 
        rb_raise(rb_eArgError, "Non-empty array expected.");
    
/*     Check_Type(*array1->ptr, T_ARRAY); */
    array2 = *RARRAY_PTR(arg);
    Check_Type(array2, T_ARRAY);
/*     array2 = RARRAY(*array1->ptr); */
/*     size2 = array2->len; */
    size2 = RARRAY_LEN(array2);
    if (size2<1) 
        rb_raise(rb_eArgError, "Matrix must have at least 1 column.");

    rmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size1, size2,
            NULL);
/*     ratlas_store_float_row(size2, array2->ptr, rmat, 0); */
    ratlas_store_float_row(size2, RARRAY_PTR(array2), rmat, 0);
    for (i = 1; i < size1; i++)
    {
        array2 = RARRAY_PTR(arg)[i];
        ratlas_store_float_row(size2, RARRAY_PTR(array2), rmat, i);
    }
    
    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    return storage;
}




VALUE ratlas_memadr(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    return INT2NUM((unsigned long int) mat->data);
}




VALUE ratlas_memsize(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    return INT2NUM(mat->memsize);
}




VALUE ratlas_matrix_from_memadr(VALUE self, VALUE nrow, VALUE ncol,
               VALUE memadr)
{
    Ratlas_Matrix *mat;

    Check_Type(nrow, T_FIXNUM);    
    Check_Type(ncol, T_FIXNUM);    
    Check_Type(memadr, T_FIXNUM);

    mat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE,     FIX2INT(nrow),
                   FIX2INT(ncol), NULL);
    mat->data = (void *) FIX2INT(memadr);
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_ptr_free,
            mat);
}




VALUE ratlas_cmatrix_new(VALUE self, VALUE arg)
{
/*     struct RArray *array1, *array2; */
    VALUE storage, array2;
    Ratlas_Matrix *rmat;
    int size1=0, size2=0, i=0;

    Check_Type(arg, T_ARRAY);
    size1 = RARRAY_LEN(arg);
    if (size1<1) 
        rb_raise(rb_eArgError, "Non-empty array expected.");
    
    array2 = *RARRAY_PTR(arg);
    Check_Type(array2, T_ARRAY);
    size2 = RARRAY_LEN(array2);
    if (size2<1) 
        rb_raise(rb_eArgError, "Matrix must have at least 1 column.");

    rmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
            RATLAS_DENSE, size1, size2, NULL);
    ratlas_store_complex_row(size2, RARRAY_PTR(array2), rmat, 0);
    for (i = 1; i < size1; i++)
    {
        array2 = RARRAY_PTR(arg)[i];
        ratlas_store_complex_row(size2, RARRAY_PTR(array2), rmat, i);
    }
    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);
    return storage;
}




VALUE ratlas_vector_new(VALUE self, VALUE arg)
{
/*     struct RArray *array; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size;

    Check_Type(arg, T_ARRAY);
/*     array = RARRAY(arg); */
    size = RARRAY_LEN(arg);
    if (size<1) 
        rb_raise(rb_eArgError, "Non-empty array expected.");
    
    rmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size, 1, NULL);

    ratlas_store_float_column(size, RARRAY_PTR(arg), rmat, 0);
     storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);
    return storage;
}




VALUE ratlas_cvector_new(VALUE self, VALUE arg)
{
/*     struct RArray *array; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size;

    Check_Type(arg, T_ARRAY);
/*     array = RARRAY(arg); */
    size = RARRAY_LEN(arg);
    if (size<1) rb_raise(rb_eArgError, "Non-empty array expected.");
    rmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE, size, 1, NULL);
    ratlas_store_complex_column(size, RARRAY_PTR(arg), rmat, 0);
    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);
    return storage;
}




VALUE ratlas_zeros_new(VALUE self, VALUE arg)
{
/*     struct RArray *carg; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size1, size2, arglen;
    unsigned long int nelem;
    double *zeros;

    Check_Type(arg, T_ARRAY);
/*     carg = RARRAY(arg); */
    arglen = RARRAY_LEN(arg);
    if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
    Check_Type(*RARRAY_PTR(arg), T_FIXNUM);
    size1 = NUM2INT(*RARRAY_PTR(arg));
    if (arglen == 2) {
        Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
        size2 = NUM2INT(RARRAY_PTR(arg)[1]);
    } else {
        size2 = 1;
    }

    nelem = size1*size2;
    if (nelem < 1) rb_raise(rb_eArgError, "Cannot make empty matrix.");
    
    zeros = xcalloc(size1*size2, sizeof(double));
    rmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size1, size2,
            (void *) zeros);
    free(zeros);

    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    return storage;
}
                    


VALUE ratlas_czeros_new(VALUE self, VALUE arg)
{
/*     struct RArray *carg; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size1, size2, arglen;
    unsigned long int nelem;
    double *zeros;

    Check_Type(arg, T_ARRAY);
/*     carg = RARRAY(arg); */
    arglen = RARRAY_LEN(arg);
    if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
    Check_Type(*RARRAY_PTR(arg), T_FIXNUM);
    size1 = NUM2INT(*RARRAY_PTR(arg));
    if (arglen == 2) {
        Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
        size2 = NUM2INT(RARRAY_PTR(arg)[1]);
    } else {
        size2 = 1;
    }

    nelem = size1*size2;
    if (nelem < 1) rb_raise(rb_eArgError, "Cannot make empty matrix.");
    
    zeros = xcalloc(2*size1*size2, sizeof(double));
    rmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE, size1, size2,
            (void *) zeros);
    free(zeros);

    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    return storage;
}




VALUE ratlas_ones_new(VALUE self, VALUE arg)
{
/*     struct RArray *carg; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size1, size2, arglen;
    unsigned long int nelem;
    double *zeros;

    Check_Type(arg, T_ARRAY);
    arglen = RARRAY_LEN(arg);
    if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
    Check_Type(*RARRAY_PTR(arg), T_FIXNUM);
    size1 = NUM2INT(*RARRAY_PTR(arg));
    if (arglen == 2) {
        Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
        size2 = NUM2INT(RARRAY_PTR(arg)[1]);
    } else {
        size2 = 1;
    }

    nelem = size1*size2;
    if (nelem < 1) rb_raise(rb_eArgError, "Cannot make empty matrix.");
    
    zeros = xcalloc(size1*size2, sizeof(double));
    rmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size1, size2,
            (void *) zeros);
    free(zeros);

    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    ratlas_add_bang(self, storage, rb_float_new(1.0));

    return storage;
}




VALUE ratlas_eye_new(VALUE self, VALUE arg)
{
/*     struct RArray *carg; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size1, size2, ndiag, arglen;
    unsigned long int nelem, i;
    double *zeros, *r;

    Check_Type(arg, T_ARRAY);
/*     carg = RARRAY(arg); */
    arglen = RARRAY_LEN(arg);
    if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
    Check_Type(*RARRAY_PTR(arg), T_FIXNUM);
    size1 = NUM2INT(*RARRAY_PTR(arg));
    if (arglen == 2) {
        Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
        size2 = NUM2INT(RARRAY_PTR(arg)[1]);
    } else {
        size2 = 1;
    }

    nelem = size1*size2;
    if (nelem < 1) rb_raise(rb_eArgError, "Cannot make empty matrix.");
    
    zeros = xcalloc(size1*size2, sizeof(double));
    rmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size1, size2,
            (void *) zeros);
    free(zeros);

    ndiag = ratlas_imin(size1, size2);

    r = (double *) rmat->data;
    for (i = 0; i < ndiag; i++) {
        r[rmat->nrow*i + i] = 1.0;
    }
    
    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    return storage;
}



VALUE ratlas_indgen_bang(VALUE self, VALUE stor, VALUE start)
{
    Ratlas_Matrix *mat;
    int i, j;
    unsigned long int pos;
    double val, cstart, *r;
    Ratlas_Complex *c;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    ratlas_assertnum(start);
    cstart = NUM2DBL(start);

    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < mat->nrow; i++)
            {
                for (j = 0; j < mat->ncol; j++)
                {
                    pos = j*mat->nrow + i;
                    val = i*mat->ncol + j + cstart;
                    r[pos] = val;
                }
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < mat->nrow; i++)
            {
                for (j = 0; j < mat->ncol; j++)
                {
                    pos = j*mat->nrow + i;
                    val = i*mat->ncol + j + cstart;
                    c[pos].r = val;
                    c[pos].i = 0.0;
                }
            }
            break;
    }
    return stor;
}



VALUE ratlas_cones_new(VALUE self, VALUE arg)
{
/*     struct RArray *carg; */
    VALUE storage;
    Ratlas_Matrix *rmat;
    int size1, size2, arglen;
    unsigned long int nelem;
    double *zeros;

    Check_Type(arg, T_ARRAY);
    arglen = RARRAY_LEN(arg);
    if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
    Check_Type(*RARRAY_PTR(arg), T_FIXNUM);
    size1 = NUM2INT(*RARRAY_PTR(arg));
    if (arglen == 2) {
        Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
        size2 = NUM2INT(RARRAY_PTR(arg)[1]);
    } else {
        size2 = 1;
    }

    nelem = size1*size2;
    if (nelem < 1) rb_raise(rb_eArgError, "Cannot make empty matrix.");
    
    zeros = xcalloc(size1*size2, sizeof(Ratlas_Complex));
    rmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE, size1, size2,
            (void *) zeros);
    free(zeros);

    storage = Data_Wrap_Struct(ratlas_storage_class, 0,
            ratlas_matrix_free, rmat);

    ratlas_add_bang(self, storage, ratlas_rb_complex_new(1.0, 0.0));

    return storage;
}




/*
 * Create a nxn matrix from the n-length vector with vector elements
 * on the diagonal. */
VALUE ratlas_vec2diag(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat, *vec;
    int i, pos;
    VALUE newstor=Qnil;
    double *zeros, *rm, *rv;
    Ratlas_Complex *cm, *cv;

    Data_Get_Struct(stor, Ratlas_Matrix, vec);
    switch (vec->type) {
        case RATLAS_DFLOAT:
            zeros = xcalloc(vec->nrow*vec->nrow, sizeof(double));
            mat = ratlas_matrix_alloc(RATLAS_DFLOAT,
                RATLAS_DENSE, vec->nrow, vec->nrow, zeros);
            free(zeros);
            rm = (double *) mat->data;
            rv = (double *) vec->data;
            for (i = 0; i < vec->nrow; i++)
            {
                rm[mat->nrow*i + i] = rv[i];
            }
            newstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, mat);
            break;
        case RATLAS_DCOMPLEX:
            zeros = xcalloc(vec->nrow*vec->ncol*2, sizeof(double));
            mat = ratlas_matrix_alloc(RATLAS_DCOMPLEX,
                    RATLAS_DENSE, vec->nrow, vec->nrow, 
                    zeros);
            free(zeros);
            cm = (Ratlas_Complex *) mat->data;
            cv = (Ratlas_Complex *) vec->data;
            for (i = 0; i < vec->nrow; i++)
            {
                pos = (mat->nrow*i + i);
                cm[pos].r = cv[i].r;
                cm[pos].i = cv[i].i;
            }
            newstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, mat);
            break;
    }
    return newstor;
}



/*
 * Read the diagonal of a full storage matrix and return as a vector. */ 
VALUE ratlas_diag2vec(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat, *vec;
    int ndiag, i;
    VALUE newstor=Qnil;
    double *rv, *rm;
    Ratlas_Complex *cv, *cm;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    ndiag = ratlas_imin(mat->nrow, mat->ncol);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            vec = ratlas_matrix_alloc(RATLAS_DFLOAT,
                    RATLAS_DENSE, ndiag, 1, NULL);
            rv = (double *) vec->data;
            rm = (double *) mat->data;
            for (i = 0; i < ndiag; i++)
            {
                rv[i] = rm[mat->nrow*i + i];
            }
            newstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, vec);
            break;
        case RATLAS_DCOMPLEX:
            vec = ratlas_matrix_alloc(RATLAS_DCOMPLEX,
                    RATLAS_DENSE, ndiag, 1, NULL);
            cv = (Ratlas_Complex *) vec->data;
            cm = (Ratlas_Complex *) mat->data;
            for (i = 0; i < ndiag; i++)
            {
                cv[i].r = cm[mat->nrow*i + i].r;
                cv[i].i = cm[mat->nrow*i + i].i;
            }
            newstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, vec);
            break;
    }
    return newstor;
}




VALUE ratlas_diag_new(VALUE self, VALUE stor)
{
    unsigned long int i;
    int ndiag;
    Ratlas_Matrix *mat, *dmat;
    VALUE dstor=Qnil;
    double *rd, *rm;
    Ratlas_Complex *cd, *cm;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);

    ndiag = ratlas_imin(mat->nrow, mat->ncol);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            dmat = ratlas_matrix_alloc(RATLAS_DFLOAT,
                    RATLAS_DENSE, ndiag, 1, NULL);
            rd = (double *) dmat->data;
            rm = (double *) mat->data;
            for (i = 0; i < ndiag; i++)
            {
                rd[i] = rm[mat->nrow*i + i];
            }
            dstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, dmat);
            break;
        case RATLAS_DCOMPLEX:
            dmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX,
                    RATLAS_DENSE, ndiag, 1, NULL);
            cd = (Ratlas_Complex *) dmat->data;
            cm = (Ratlas_Complex *) mat->data;
            for (i = 0; i < ndiag; i++)
            {
                cd[i].r = cm[mat->nrow*i + i].r;
                cd[i].i = cm[mat->nrow*i + i].i;
            }
            dstor = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, dmat);
            break;
    }
    
    return dstor;
}




VALUE ratlas_vector_to_a (VALUE self, VALUE storage)
{
    Ratlas_Matrix *rmat;
    int i;
    double *r, re, im;
    Ratlas_Complex *c;
    VALUE array;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    if (rmat == NULL)
        rb_raise(rb_eArgError, "Expected RatlasStorage-object.");
    if (rmat->ncol != 1) 
        rb_raise(rb_eArgError, "Storage object should have 1 column.");
    array = rb_ary_new2(rmat->nrow);

    switch (rmat->type) {
        case RATLAS_DFLOAT:
            r = (double *) rmat->data;
            for (i = 0; i < rmat->nrow; i++)
                rb_ary_push(array, rb_float_new(r[i])); 
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) rmat->data;
            for (i = 0; i < rmat->nrow; i++)
            {
                re = c[i].r;
                im = c[i].i;
                rb_ary_push(array, ratlas_rb_complex_new(re, im));
            }
            break;
    }
    
    return array;
}




VALUE ratlas_matrix_to_a (VALUE self, VALUE storage)
{
    Ratlas_Matrix *rmat;
    VALUE array1, array2;
    double *r, re, im;
    int i, j, m;
    Ratlas_Complex *c;
    
    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    if (rmat == NULL)
        rb_raise(rb_eArgError, "Expected RatlasStorage-object.");
    if (rmat->ncol < 1)
        rb_raise(rb_eArgError, "Storage object has no columns.");
        
    array1 = rb_ary_new2(rmat->nrow);
    m = rmat->nrow;
    switch (rmat->type) {
        case RATLAS_DFLOAT:
            r = (double *) rmat->data;
            for (i = 0; i < rmat->nrow; i++)
            {
                array2 = rb_ary_new2(rmat->ncol);
                for (j = 0; j < rmat->ncol; j++)
                {
                    rb_ary_push(array2, rb_float_new(r[i+m*j]));
                }
                rb_ary_push(array1, array2);
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) rmat->data;
            for (i = 0; i < rmat->nrow; i++)
            {
                array2 = rb_ary_new2(rmat->ncol);
                for (j = 0; j < rmat->ncol; j++)
                {
                    re = c[i+m*j].r;
                    im = c[i+m*j].i;
                    rb_ary_push(array2, ratlas_rb_complex_new(re, im));
                }
                rb_ary_push(array1, array2);
            }
            break;
        default:
            rb_raise(rb_eArgError, "Unsupported storage type.");
            break;
    }
    return array1;
}




VALUE ratlas_get_column(VALUE self, VALUE storage, VALUE colidx)
{
    Ratlas_Matrix *mat, *retmat;
    int j;
    unsigned long int nbytes;
    double *rsrc;
    Ratlas_Complex *csrc;
    
    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    Check_Type(colidx, T_FIXNUM);
    j = FIX2INT(colidx);
    if (j < 0) j += mat->ncol;
    if (j > mat->ncol - 1)
        rb_raise(rb_eIndexError, "Invalid column index.");

    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, 
            mat->nrow, 1, NULL);
    nbytes = mat->nrow * ratlas_sizeof[mat->type];
    switch(mat->type) {
        case RATLAS_DFLOAT:
            rsrc = (double *)mat->data + j*mat->nrow;
            memcpy((void *)retmat->data, (void *)rsrc, nbytes);
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *)mat->data + j*mat->nrow;
            memcpy((void *)retmat->data, (void *)csrc, nbytes);
            break;
        default:
            rb_raise(rb_eArgError, "Unsupported storage type.");
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, retmat);
}




VALUE ratlas_get_columns(VALUE self, VALUE storage, VALUE colidx)
{
/*     struct RArray *arr; */
    Ratlas_Matrix *mat, *retmat;
    int i, j, arrlen;
    VALUE col;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;
    unsigned long int nbytes, colskip=0;
            
    Check_Type(colidx, T_ARRAY);
    arrlen = RARRAY_LEN(colidx);
    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, mat->nrow,
            arrlen, NULL);
    nbytes = mat->nrow * ratlas_sizeof[mat->type];
    colskip = mat->nrow;

    switch (mat->type) {
        case RATLAS_DFLOAT:
            rdest = (double *)retmat->data;
            for (i = 0; i < arrlen; i++)
            {
                col = RARRAY_PTR(colidx)[i];
                if (TYPE(col) != T_FIXNUM) goto error;
                j = FIX2INT(col);
                if (j < 0) j += mat->ncol;
                if (j > mat->ncol -1) goto error;    
                rsrc = (double *)mat->data + j*colskip;
                rdest += i * colskip;
                memcpy((void *)rdest, (void *)rsrc, nbytes);
            }
            break;
        case RATLAS_DCOMPLEX:
               cdest = (Ratlas_Complex *)retmat->data;
            for (i = 0; i < arrlen; i++)
            {
                col = RARRAY_PTR(colidx)[i];
                if (TYPE(col) != T_FIXNUM) goto error;
                j = FIX2INT(col);
                if (j < 0) j += mat->ncol;
                if (j > mat->ncol -1) goto error;    
                csrc = (Ratlas_Complex *)mat->data + j*colskip;
                cdest += i * colskip;
                memcpy((void *)cdest, (void *)csrc, nbytes);
            }
            break;
 }

    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retmat);
    
error:
    ratlas_matrix_free(retmat);
    rb_raise(rb_eArgError, "Could not get columns.  Illegal indices?");
}




VALUE ratlas_get_row(VALUE self, VALUE storage, VALUE rowidx)
{
    Ratlas_Matrix *mat, *retmat;
    int i, j;
    unsigned long int offset;
    double *rdest, *rsrc;
    Ratlas_Complex *cdest, *csrc;

    Data_Get_Struct(storage, Ratlas_Matrix, mat);
        Check_Type(rowidx, T_FIXNUM);
    i = FIX2INT(rowidx);
    if (i < 0) i += mat->nrow;
    if (i > mat->nrow -1)
        rb_raise(rb_eIndexError, "Invalid row index.");
    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, 1, mat->ncol, 
            NULL);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            rsrc = (double *) mat->data;
            rdest = (double *) retmat->data;
            for (j = 0; j < mat->ncol; j++)
            {
                offset = i + j*mat->nrow;
                rdest[j] = rsrc[offset];
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *) mat->data;
            cdest = (Ratlas_Complex *) retmat->data;
            for (j = 0; j < mat->ncol; j++)
            {
                offset = i + j*mat->nrow;
                cdest[j].r = csrc[offset].r;
                cdest[j].i = csrc[offset].i;
            }
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retmat);
}





VALUE ratlas_get_rows(VALUE self, VALUE storage, VALUE rowidx)
{
/*     struct RArray *arr; */
    int i, j, k, arrlen;
    unsigned long int isrc, idest;
    VALUE row;
    Ratlas_Matrix *mat, *retmat;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;
    
    Check_Type(rowidx, T_ARRAY);
/*     arr = RARRAY(rowidx); */
    arrlen = RARRAY_LEN(rowidx);
    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, arrlen, 
            mat->ncol, NULL);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            rsrc = (double *)mat->data;
            rdest = (double *)retmat->data;
            for (k = 0; k < arrlen; k++)
            {
                row = RARRAY_PTR(rowidx)[k];
                if (TYPE(row) != T_FIXNUM) goto error;
                i = FIX2INT(row);
                if (i < 0) i += mat->nrow;
                if (i > mat->nrow -1) goto error;    
                for (j = 0; j < mat->ncol; j++)
                {
                    isrc = i + j*mat->nrow;
                    idest = k + j*retmat->nrow;
                    rdest[idest] = rsrc[isrc];
                }
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *)mat->data;
            cdest = (Ratlas_Complex *)retmat->data;
            for (k = 0; k < arrlen; k++)
            {
                row = RARRAY_PTR(rowidx)[k];
                if (TYPE(row) != T_FIXNUM) goto error;
                i = FIX2INT(row);
                if (i < 0) i += mat->nrow;
                if (i > mat->nrow -1) goto error;    
                for (j = 0; j < mat->ncol; j++)
                {
                    isrc = i + j*mat->nrow;
                    idest = k + j*retmat->nrow;
                    cdest[idest].r = csrc[isrc].r;
                    cdest[idest].i = csrc[isrc].i;
                }
            }
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retmat);

error:
    ratlas_matrix_free(retmat);
    rb_raise(rb_eArgError, "Could not get rows.  Illegal index?");
}




VALUE ratlas_set_column_bang(VALUE self, VALUE storage, VALUE colidx,
        VALUE newval)
{
    Ratlas_Matrix *rmat;
/*     struct RArray *cnewval; */
    int i, j;
    unsigned long int offset;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    
    Check_Type(colidx, T_FIXNUM);
    j = FIX2INT(colidx);
    if (j < 0) j += rmat->ncol;
    if (j > rmat->ncol - 1)
        rb_raise(rb_eIndexError, "Invalid column index.");
/*     cnewval = RARRAY(newval); */
    if (RARRAY_LEN(newval) != rmat->nrow)
        rb_raise(rb_eArgError, "Wrong column size.");
    switch (rmat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < rmat->nrow; i++) {
                ratlas_assertnum(RARRAY_PTR(newval)[i]);
                offset = i + j*rmat->nrow;
                ratlas_set_float(rmat->data, offset, RARRAY_PTR(newval)[i]);
            }
            break;
        case RATLAS_DCOMPLEX:
            for (i = 0; i < rmat->nrow; i++) {
                offset = i + j*rmat->nrow;
                ratlas_set_complex(rmat->data, offset, RARRAY_PTR(newval)[i]);
            }
            break;
    }
    return storage;    
}




VALUE ratlas_set_columns_bang(VALUE self, VALUE storage, VALUE colidx,
        VALUE newval)
{
/*     struct RArray *ccolidx, *cnewval; */
    int j;

    Check_Type(colidx, T_ARRAY);
/*     ccolidx = RARRAY(colidx); */
/*     cnewval = RARRAY(newval); */
    if (RARRAY_LEN(newval) != RARRAY_LEN(colidx))
        rb_raise(rb_eArgError, "Index and argument size must match.");
    
    for (j = 0;  j < RARRAY_LEN(colidx); j++)
    {
        ratlas_set_column_bang(self, storage, RARRAY_PTR(colidx)[j],
                RARRAY_PTR(newval)[j]);
    }    
    return storage;
        
}




VALUE ratlas_set_row_bang(VALUE self, VALUE storage, VALUE rowidx, VALUE newval)
{
    Ratlas_Matrix *rmat;
/*     struct RArray *cnewval; */
    int i, j;
    unsigned long int offset;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    
    Check_Type(rowidx, T_FIXNUM);
    i = FIX2INT(rowidx);
    if (i < 0) i += rmat->nrow;
    if (i > rmat->nrow - 1)
        rb_raise(rb_eIndexError, "Invalid row index.");
/*     cnewval = RARRAY(newval); */
    if (RARRAY_LEN(newval) != rmat->ncol)
        rb_raise(rb_eArgError, "Wrong row size.");
    switch (rmat->type) {
        case RATLAS_DFLOAT:
            for (j = 0; j < rmat->ncol; j++) {
                ratlas_assertnum(RARRAY_PTR(newval)[j]);
                offset = i + j*rmat->nrow;
                ratlas_set_float(rmat->data, offset, RARRAY_PTR(newval)[j]);
            }
            break;
        case RATLAS_DCOMPLEX:
            for (j = 0; j < rmat->ncol; j++) {
                offset = i + j*rmat->nrow;
                ratlas_set_complex(rmat->data, offset, RARRAY_PTR(newval)[j]);
            }
            break;
    }
    return storage;

}




VALUE ratlas_set_rows_bang(VALUE self, VALUE storage, VALUE rowidx, VALUE newval)
{
/*     struct RArray *crowidx, *cnewval; */
    int i;

    Check_Type(rowidx, T_ARRAY);
/*     crowidx = RARRAY(rowidx); */
/*     cnewval = RARRAY(newval); */
    if (RARRAY_LEN(newval) != RARRAY_LEN(rowidx))
        rb_raise(rb_eArgError, "Index and argument size must match.");
    
    for (i = 0;  i < RARRAY_LEN(rowidx); i++) {
        ratlas_set_row_bang(self, storage, RARRAY_PTR(rowidx)[i],
                RARRAY_PTR(newval)[i]);
    }    
    return storage;
        
}




VALUE ratlas_size(VALUE self, VALUE storage)
{
    VALUE i, j;
    Ratlas_Matrix *rmat;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    i = INT2FIX(rmat->nrow);
    j = INT2FIX(rmat->ncol);
    return rb_ary_new3(2, i, j);
}



/*
 * Get one value from Ratlas_Matrix. */
VALUE ratlas_get_one(VALUE self, VALUE storage, VALUE index)
{
    long int nelem, idx;
    Ratlas_Matrix *rmat;
    
    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    nelem = rmat->nrow * rmat->ncol;
    idx = FIX2INT(index);
    if (idx < 0)
        idx += nelem;
    if ((idx > nelem - 1) || (idx < 0))
        rb_raise(rb_eArgError, "Index out of range.");
    idx = ratlas_swap_major(rmat->nrow, rmat->ncol, idx);
    switch (rmat->type) {
        case RATLAS_DCOMPLEX:
            return ratlas_get_complex(rmat->data, idx);
        case RATLAS_DFLOAT:
            return ratlas_get_float(rmat->data, idx);
        default:
            rb_raise(rb_eArgError, "Unsupported storage type.");
    }
}



/*
 * Using a Array of indices, get values from Ratlas_Matrix.  Return as
 * vector.  The index run continously in Row Major order.  */
VALUE ratlas_get_many(VALUE self, VALUE storage, VALUE ary)
{
/*     struct RArray *ary; */
    VALUE idx;
    int i, j, arylen;
    long int nelem;
    Ratlas_Matrix *retvec, *rmat;
    double *rv, *rm;
    Ratlas_Complex *cv, *cm;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    nelem = rmat->nrow * rmat->ncol;
    
/*     ary = RARRAY(idx); */
    arylen = RARRAY_LEN(ary);
    switch (rmat->type) { 
        case  RATLAS_DFLOAT:
            retvec = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                    RATLAS_DENSE, arylen, 1, NULL);
            rv = (double *) retvec->data;
            rm = (double *) rmat->data;
            for (j = 0; j < arylen; j++)
            {
                idx = RARRAY_PTR(ary)[j];
                Check_Type(idx, T_FIXNUM);
                i = FIX2INT(idx);
                if (i < 0) i += nelem;
                if ((i > nelem - 1) || (i < 0))
                    rb_raise(rb_eArgError, "Index out of range.");
                i = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                rv[j] = rm[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            retvec = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
                    RATLAS_DENSE, arylen, 1, NULL);
            cv = (Ratlas_Complex *) retvec->data;
            cm = (Ratlas_Complex *) rmat->data;
            for (j = 0; j < arylen; j++) {
                idx = RARRAY_PTR(ary)[j];
                Check_Type(idx, T_FIXNUM);
                i = FIX2INT(idx);
                if (i < 0) i += nelem;
                if ((i > nelem - 1) || (i < 0))
		    rb_raise(rb_eArgError, "Index out of range.");
                i = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                cv[j].r = cm[i].r;
                cv[j].i = cm[i].i;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retvec);
}




VALUE ratlas_get_by_range(VALUE self, VALUE storage, VALUE range)
{
    int first, last, i, j=0, range_len, idx;
    unsigned long int nelem;
    Ratlas_Matrix *retvec, *rmat;
    double *rv, *rm;
    Ratlas_Complex *cv, *cm;

    if (rb_obj_is_kind_of(range, range_class) != Qtrue)
        rb_raise(rb_eArgError, "Expect Range.");
    first = FIX2INT(rb_funcall(range, id_first, 0));
    last = FIX2INT(rb_funcall(range, id_last, 0));

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    nelem = rmat->nrow * rmat->ncol;
    if (first < 0) first += nelem;
    if (last < 0) last += nelem;
    if (first < 0 || last < 0) rb_raise(rb_eArgError, "Illegal range.");
    if ( rb_funcall(range, id_exclude_end, 0) == Qtrue)
        last -= 1;
    if (first > last) rb_raise(rb_eArgError, "Expect an increasing range.");
    if (last > nelem -1 || first > nelem -1) 
        rb_raise(rb_eArgError, "Illegal range.");
    range_len = last - first + 1;

    switch (rmat->type) { 
        case  RATLAS_DFLOAT:
            retvec = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                    RATLAS_DENSE, range_len, 1, NULL);
            rv = (double *) retvec->data;
            rm = (double *) rmat->data;
            for (i = first; i <= last; i++) {
                idx = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                rv[j++] = rm[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            retvec = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
                    RATLAS_DENSE, range_len, 1, NULL);
            cv = (Ratlas_Complex *) retvec->data;
            cm = (Ratlas_Complex *) rmat->data;
            for (i = first; i <= last ; i++) {
                idx = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                cv[j].r = cm[idx].r;
                cv[j].i = cm[idx].i;
                j++;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retvec);
}




/*
 * Get one value from 2d array based on row and column index. */
VALUE ratlas_get2d_one(VALUE self, VALUE storage, VALUE row, VALUE col)
{
    Ratlas_Matrix *rmat;
    unsigned long int idx;
    int m, n;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    Check_Type(row, T_FIXNUM);
    Check_Type(col, T_FIXNUM);    
    m = FIX2INT(row);
    n = FIX2INT(col);
    if (m < 0) m += rmat->nrow;
    if (m > rmat->nrow - 1) rb_raise(rb_eArgError, "Index out of range.");
    if (n < 0) n += rmat->ncol;
    if (n > rmat->ncol - 1) rb_raise(rb_eArgError, "Index out of range.");

    idx = n + m * rmat->ncol;
    return ratlas_get_one(self, storage, INT2FIX(idx));
}





VALUE ratlas_get2d_many(VALUE self, VALUE storage, VALUE rows, VALUE cols)
{
    Ratlas_Matrix *mat, *retmat;
/*     struct RArray *crows, *ccols; */
    int i, j, irow, icol, rowlen, collen;
    unsigned long int isrc, idest;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;

    Data_Get_Struct(storage, Ratlas_Matrix, mat);
/*     crows = RARRAY(rows); */
/*     ccols = RARRAY(cols); */
    rowlen = RARRAY_LEN(rows);
    collen = RARRAY_LEN(cols);
    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, rowlen, collen,
		    NULL);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            rsrc = (double *) mat->data;
            rdest = (double *) retmat->data;
            for (i = 0; i < rowlen; i++) {
                irow = FIX2INT(RARRAY_PTR(rows)[i]);
                if (irow > (mat->nrow-1)) goto error;
                for (j = 0; j < collen; j++) {
                    icol = FIX2INT(RARRAY_PTR(cols)[j]);
                    if (icol > (mat->ncol-1)) goto error;
                    idest = i + j*retmat->nrow;
                    isrc = irow + icol*mat->nrow;
                    rdest[idest] = rsrc[isrc];
                }
            }
            break;
        case RATLAS_DCOMPLEX:
            for (i = 0; i < rowlen; i++) {
                csrc = (Ratlas_Complex *) mat->data;
                cdest = (Ratlas_Complex *) retmat->data;
                irow = FIX2INT(RARRAY_PTR(rows)[i]);
                if (irow > (mat->nrow-1)) goto error;
                for (j = 0; j < collen; j++) {
                    icol = FIX2INT(RARRAY_PTR(cols)[j]);
                    if (icol > (mat->ncol-1)) goto error;
                    idest = i + j*retmat->nrow;
                    isrc = irow + icol*mat->nrow;
                    cdest[idest].r = csrc[isrc].r;
                    cdest[idest].i = csrc[isrc].i;
                }
            }
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retmat);

error:
    ratlas_matrix_free(retmat);
    rb_raise(rb_eArgError, "Unable to get elements.  Illegal indices?");
}




VALUE ratlas_get2d_by_range(VALUE self, VALUE storage, VALUE rowrange, 
        VALUE colrange)
{
    int colfirst, collast, rowfirst, rowlast, i, j, ir=0, ic=0, 
        colrange_len, rowrange_len;
    unsigned long int isrc, idest;
    Ratlas_Matrix *retmat, *rmat;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;

    if ((rb_obj_is_kind_of(rowrange, range_class) != Qtrue) ||
        (rb_obj_is_kind_of(colrange, range_class) != Qtrue) )
        rb_raise(rb_eArgError, "Expect Range.");
    colfirst = FIX2INT(rb_funcall(colrange, id_first, 0));
    collast = FIX2INT(rb_funcall(colrange, id_last, 0));
    rowfirst = FIX2INT(rb_funcall(rowrange, id_first, 0));
    rowlast = FIX2INT(rb_funcall(rowrange, id_last, 0));

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    if (colfirst < 0) colfirst += rmat->ncol;
    if (rowfirst < 0) rowfirst += rmat->nrow;
    if (collast < 0) collast += rmat->ncol;
    if (rowlast < 0) rowlast += rmat->nrow;
    if ( rb_funcall(colrange, id_exclude_end, 0) == Qtrue) collast -= 1;
    if ( rb_funcall(rowrange, id_exclude_end, 0) == Qtrue) rowlast -= 1;

    if (colfirst > collast)
               rb_raise(rb_eArgError, "Expect an increasing range.");
    if (rowfirst > rowlast)
               rb_raise(rb_eArgError, "Expect an increasing range.");
    if (collast > (rmat->ncol -1))
               rb_raise(rb_eArgError, "Range exceeds elements.");
    if (rowlast > (rmat->nrow -1))
               rb_raise(rb_eArgError, "Range exceeds elements.");

    colrange_len = collast - colfirst + 1;
    rowrange_len = rowlast - rowfirst + 1;
    switch (rmat->type) { 
        case  RATLAS_DFLOAT:
            retmat = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                RATLAS_DENSE, rowrange_len, colrange_len, NULL);
            rsrc = (double *) rmat->data;
            rdest = (double *) retmat->data;
            for (i = rowfirst; i <= rowlast; i++) {
                ic = 0;
                for (j = colfirst; j <= collast; j++) {
                    idest = ir + ic * retmat->nrow;
                    isrc = i + j * rmat->nrow;
                    rdest[idest] = rsrc[isrc];
                    ic++;
                }
                ir++;
            }
            break;
        case RATLAS_DCOMPLEX:
            retmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
                RATLAS_DENSE, rowrange_len, colrange_len, NULL);
            csrc = (Ratlas_Complex *) rmat->data;
            cdest = (Ratlas_Complex *) retmat->data;
            for (i = rowfirst; i <= rowlast ; i++) {
                ic = 0;
                for (j = colfirst; j <= collast; j++) {
                    idest = ir + ic * retmat->nrow;
                    isrc = i + j * rmat->nrow;
                    cdest[idest].r = csrc[isrc].r;
                    cdest[idest].i = csrc[isrc].i;
                    ic++;
                }
                ir++;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            retmat);
}





VALUE ratlas_set_one_bang(VALUE self, VALUE storage, VALUE index, VALUE newval)
{
    unsigned long int nelem, idx;
    Ratlas_Matrix *rmat;
    
    idx = FIX2INT(index);
    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    nelem = rmat->nrow * rmat->ncol;
    if (idx < 0) idx += nelem;
    if ((idx > nelem - 1) || (idx < 0))
        rb_raise(rb_eArgError, "Index out of range.");
    idx = ratlas_swap_major(rmat->nrow, rmat->ncol, idx);
    switch (rmat->type) {
        case RATLAS_DCOMPLEX:
            ratlas_set_complex(rmat->data, idx, newval);
            return storage;
        case RATLAS_DFLOAT:
            ratlas_set_float(rmat->data, idx, newval);
            return storage;
        default:
            rb_raise(rb_eArgError, "Unsupported storage type.");
    }
}



VALUE ratlas_set_many_bang(VALUE self, VALUE storage, VALUE ary, VALUE newvals)
{
/*     struct RArray *ary; */
    VALUE idx;
    int i, j, arylen;
    long int nelem;
    Ratlas_Matrix *rvec, *rmat;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;
    
    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    nelem = rmat->nrow * rmat->ncol;
/*     ary = RARRAY(idx); */
    Data_Get_Struct(newvals, Ratlas_Matrix, rvec);
    if (rvec->ncol != 1) rb_raise(rb_eArgError, "Input must be vector.");
    
    arylen = RARRAY_LEN(ary);
    switch (rmat->type) { 
        case  RATLAS_DFLOAT:
            rsrc = (double *)rvec->data;
            rdest = (double *)rmat->data;
            for (j = 0; j < arylen; j++) {
                idx = RARRAY_PTR(ary)[j];
                Check_Type(idx, T_FIXNUM);
                i = FIX2INT(idx);
                if (i < 0) i += nelem;
                if ((i > nelem - 1) || (i < 0))
                    rb_raise(rb_eArgError, "Index out of range.");
                i = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                rdest[i] = rsrc[j];
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *) rvec->data;
            cdest = (Ratlas_Complex *) rmat->data;
            for (j = 0; j < arylen; j++) {
                idx = RARRAY_PTR(ary)[j];
                Check_Type(idx, T_FIXNUM);
                i = FIX2INT(idx);
                if (i < 0) i += nelem;
                if ((i > nelem - 1) || (i < 0))
                    rb_raise(rb_eArgError, "Index out of range.");
                i = ratlas_swap_major(rmat->nrow, rmat->ncol, i);
                cdest[i].r = csrc[j].r;
                cdest[i].i = csrc[j].i;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return storage;
}




VALUE ratlas_set_by_range_bang(VALUE self, VALUE storage, VALUE range, 
        VALUE newvals)
{
    int first, last, i, j=0, range_len, idx;
    long int nmat, nnewmat;
    Ratlas_Matrix *newmat, *mat;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;

    if (rb_obj_is_kind_of(range, range_class) != Qtrue)
        rb_raise(rb_eArgError, "Expect Range.");
    first = FIX2INT(rb_funcall(range, id_first, 0));
    last = FIX2INT(rb_funcall(range, id_last, 0));

    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    Data_Get_Struct(newvals, Ratlas_Matrix, newmat);

    if ( (mat->type != newmat->type) || 
            (mat->matrixtype != newmat->matrixtype) )
        rb_raise(rb_eArgError, "Type must be equal.");
    nmat = mat->nrow * mat->ncol;
    nnewmat = newmat->nrow * newmat->ncol;
    if (first < 0) first += nmat;
    if (last < 0) last += nmat;
    if (first < 0 || last < 0) rb_raise(rb_eArgError, "Illegal range.");
    if ( rb_funcall(range, id_exclude_end, 0) == Qtrue)
        last -= 1;
    if (first > last) rb_raise(rb_eArgError, "Expect an increasing range.");
    if (last > nmat -1) rb_raise(rb_eArgError, "Range exceeds elements.");
    range_len = last - first + 1;
    if (range_len != nnewmat) rb_raise(rb_eArgError, 
            "Incorrect size of new values array.");

    switch (mat->type) { 
        case  RATLAS_DFLOAT:
            rsrc = (double *) newmat->data;
            rdest = (double *) mat->data;
            for (i = first; i <= last; i++) {
                idx = ratlas_swap_major(mat->nrow, mat->ncol, i);
                rdest[idx] = rsrc[j];
                j++;
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *) newmat->data;
            cdest = (Ratlas_Complex *) mat->data;
            for (i = first; i <= last ; i++) {
                idx = ratlas_swap_major(mat->nrow, mat->ncol, i);
                cdest[idx].r = csrc[j].r;
                cdest[idx].i = csrc[j].i;
                j++;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return storage;
}




VALUE ratlas_set2d_one_bang(VALUE self, VALUE storage, VALUE row, VALUE col,
        VALUE newval)
{
    Ratlas_Matrix *rmat;
    long int idx;
    int m, n;

    Data_Get_Struct(storage, Ratlas_Matrix, rmat);
    m = FIX2INT(row);
    n = FIX2INT(col);
    if (m < 0) m += rmat->nrow;
    if (m > rmat->nrow - 1) rb_raise(rb_eArgError, "Index out of range.");
    if (n < 0) n += rmat->ncol;
    if (n > rmat->ncol - 1) rb_raise(rb_eArgError, "Index out of range.");

    idx = n + m * rmat->ncol;
    return ratlas_set_one_bang(self, storage, INT2FIX(idx), newval);
}




VALUE ratlas_set2d_many_bang(VALUE self, VALUE storage, VALUE rows, VALUE cols,
        VALUE newval)
{
    Ratlas_Matrix *mat, *newmat;
/*     struct RArray *rowarr, *colarr; */
    int ir, ic, row, col, rowlen, collen;
    unsigned long int isrc=0, idest;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;

    Check_Type(rows, T_ARRAY);
    Check_Type(cols, T_ARRAY);
    rowlen = RARRAY_LEN(rows);
    collen = RARRAY_LEN(cols);
/*     rowarr = RARRAY(rows); */
/*     colarr = RARRAY(cols); */
    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    Data_Get_Struct(newval, Ratlas_Matrix, newmat);
    if (mat->type != newmat->type || mat->matrixtype != newmat->matrixtype)
               rb_raise(rb_eArgError, "Conflicting types for matrices.");
    if (rowlen != newmat->nrow || collen != newmat->ncol)
        rb_raise(rb_eArgError, "Incorrect size of argument.");
    switch (mat->type) {
        case RATLAS_DFLOAT:    
            rsrc = (double *) newmat->data;
            rdest = (double *) mat->data;
            for (ic=0; ic < collen; ic++) {
                Check_Type(RARRAY_PTR(cols)[ic], T_FIXNUM);
                col = FIX2INT(RARRAY_PTR(cols)[ic]);
                if (col < 0) col += mat->ncol;
                if (col > mat->ncol -1)
		       	rb_raise(rb_eArgError, "Index out of range.");
                for (ir=0; ir < rowlen; ir++) {
                    Check_Type(RARRAY_PTR(rows)[ir], T_FIXNUM);
                    row = FIX2INT(RARRAY_PTR(rows)[ir]);
                    if (row < 0) row += mat->nrow;
                    if (row > mat->nrow -1)
                        rb_raise(rb_eArgError, "Index out of range.");    
                    idest = row + col * mat->nrow;
                    rdest[idest] = rsrc[isrc];
                    isrc++;
                }
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *)newmat->data;
            cdest = (Ratlas_Complex *)mat->data;
            for (ic=0; ic < collen; ic++) {
                Check_Type(RARRAY_PTR(cols)[ic], T_FIXNUM);
                col = FIX2INT(RARRAY_PTR(cols)[ic]);
                if (col < 0) col += mat->ncol;
                if (col > mat->ncol -1)
                    rb_raise(rb_eArgError, "Index out of range.");
                for (ir=0; ir < rowlen; ir++) {
                    Check_Type(RARRAY_PTR(rows)[ir], T_FIXNUM);
                    row = FIX2INT(RARRAY_PTR(rows)[ir]);
                    if (row < 0) row += mat->nrow;
                    if (row > mat->nrow -1)
                        rb_raise(rb_eArgError, "Index out of range.");    
                    idest = row + col * mat->nrow;
                    cdest[idest].r = csrc[isrc].r;
                    cdest[idest].i = csrc[isrc].i;
                    isrc++;
                }
            }
            break;
    }
    return storage;
}




VALUE ratlas_set2d_by_range_bang(VALUE self, VALUE storage, VALUE rowrange, 
        VALUE colrange, VALUE newvals)
{
    int colfirst, collast, rowfirst, rowlast, i, j, ir=0, ic=0, 
        colrange_len, rowrange_len;
    unsigned long int isrc, idest;
    Ratlas_Matrix *mat, *newmat;
    double *rsrc, *rdest;
    Ratlas_Complex *csrc, *cdest;

    if (    (rb_obj_is_kind_of(rowrange, range_class) != Qtrue) ||
        (rb_obj_is_kind_of(colrange, range_class) != Qtrue) )
        rb_raise(rb_eArgError, "Expect Range.");
    colfirst = FIX2INT(rb_funcall(colrange, id_first, 0));
    collast = FIX2INT(rb_funcall(colrange, id_last, 0));
    rowfirst = FIX2INT(rb_funcall(rowrange, id_first, 0));
    rowlast = FIX2INT(rb_funcall(rowrange, id_last, 0));

    Data_Get_Struct(storage, Ratlas_Matrix, mat);
    Data_Get_Struct(newvals, Ratlas_Matrix, newmat);

    if (mat->type != newmat->type || mat->matrixtype != newmat->matrixtype)
        rb_raise(rb_eArgError, "Conflicting types.");
    if (colfirst < 0) colfirst += mat->ncol;
    if (rowfirst < 0) rowfirst += mat->nrow;
    if (collast < 0) collast += mat->ncol;
    if (rowlast < 0) rowlast += mat->nrow;
    if ( rb_funcall(colrange, id_exclude_end, 0) == Qtrue)
        collast -= 1;
    if ( rb_funcall(rowrange, id_exclude_end, 0) == Qtrue)
        rowlast -= 1;
    if (colfirst > collast)
               rb_raise(rb_eArgError, "Expect an increasing range.");
    if (rowfirst > rowlast)
               rb_raise(rb_eArgError, "Expect an increasing range.");
    if (collast > (mat->ncol -1))
               rb_raise(rb_eArgError, "Range exceeds elements.");
    if (rowlast > (mat->nrow -1))
               rb_raise(rb_eArgError, "Range exceeds elements.");

    colrange_len = collast - colfirst + 1;
    rowrange_len = rowlast - rowfirst + 1;
    if (rowrange_len * colrange_len != newmat->nrow *newmat->ncol)
        rb_raise(rb_eArgError, "Argument does not match ranges.");
    if (rowrange_len != newmat->nrow || colrange_len != newmat->ncol)
        rb_raise(rb_eArgError, "Incorrect size of argument.");

    switch (mat->type) { 
        case  RATLAS_DFLOAT:
            rsrc = (double *) newmat->data;
            rdest = (double *) mat->data;
            for (i = rowfirst; i <= rowlast; i++) {
                ic = 0;
                for (j = colfirst; j <= collast; j++) {
                    isrc = ir + ic * newmat->nrow;
                    idest = i + j * mat->nrow;
                    rdest[idest] = rsrc[isrc];
                    ic++;
                }
                ir++;
            }
            break;
        case RATLAS_DCOMPLEX:
            csrc = (Ratlas_Complex *) newmat->data;
            cdest = (Ratlas_Complex *) mat->data;
            for (i = rowfirst; i <= rowlast ; i++) {
                ic = 0;
                for (j = colfirst; j <= collast; j++) {
                    isrc = ir + ic * newmat->nrow;
                    idest = i + j * mat->nrow;
                    cdest[idest].r = csrc[isrc].r;
                    cdest[idest].i = csrc[isrc].i;
                    ic++;
                }
                ir++;
            }
            break;
        default:
            rb_raise(rb_eArgError, 
                "Unsupported storage type.");
    }
    return storage;
}




/*
 * Flattenes  and concatenates matrices and vectors. */
VALUE ratlas_concat(VALUE self, VALUE stor1, VALUE stor2)
{
    Ratlas_Matrix *mat1, *mat2, *newmat;
    long int nelem;

    Data_Get_Struct(stor1, Ratlas_Matrix, mat1);
    Data_Get_Struct(stor2, Ratlas_Matrix, mat2);

    if (mat1->type != mat2->type)
        rb_raise(rb_eArgError, "Can not concatenate unequal types.");
    nelem = mat1->nrow*mat1->ncol + mat2->nrow*mat2->ncol;    
    newmat = ratlas_matrix_alloc (mat1->type, mat1->matrixtype, 
            nelem, 1, NULL);
    memcpy((void *)newmat->data, (void *)mat1->data, mat1->memsize);
    memcpy((void *)newmat->data + mat1->memsize, (void *)mat2->data, 
            mat2->memsize);
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   newmat);
}




VALUE ratlas_hcat(VALUE self, VALUE stor1, VALUE stor2)
{
    Ratlas_Matrix *mat1, *mat2, *newmat;

    Data_Get_Struct(stor1, Ratlas_Matrix, mat1);
    Data_Get_Struct(stor2, Ratlas_Matrix, mat2);

    if (mat1->type != mat2->type)
        rb_raise(rb_eArgError, "Can not concatenate unequal types.");
    if (mat1->nrow != mat2->nrow)
        rb_raise(rb_eArgError, "Need equal number of rows for hcat.");
    
    newmat = ratlas_matrix_alloc (mat1->type, mat1->matrixtype, 
            mat1->nrow, mat1->ncol + mat2->ncol, NULL);
    memcpy((void *)newmat->data, (void *)mat1->data, mat1->memsize);
    memcpy((void *)newmat->data + mat1->memsize, (void *)mat2->data, 
            mat2->memsize);
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   newmat);
}





VALUE ratlas_vcat(VALUE self, VALUE stor1, VALUE stor2)
{
    Ratlas_Matrix *mat1, *mat2, *newmat;
    long int offset, offset1, offset2;
    int j, mat1_colsize, mat2_colsize, newmat_colsize, elemsize;

    Data_Get_Struct(stor1, Ratlas_Matrix, mat1);
    Data_Get_Struct(stor2, Ratlas_Matrix, mat2);

    if (mat1->type != mat2->type)
        rb_raise(rb_eArgError, "Can not concatenate unequal types.");
    if (mat1->ncol != mat2->ncol)
        rb_raise(rb_eArgError, 
                "Need equal number of columns for vcat.");
    
    newmat = ratlas_matrix_alloc (mat1->type, mat1->matrixtype, 
            mat1->nrow + mat2->nrow, mat1->ncol, NULL);
    switch (mat1->type) {
        case RATLAS_DFLOAT:
            elemsize = ratlas_sizeof[RATLAS_DFLOAT];
            break;
        case RATLAS_DCOMPLEX:
            elemsize = ratlas_sizeof[RATLAS_DCOMPLEX];
            break;
        default:
            elemsize = ratlas_sizeof[RATLAS_DFLOAT];
    }
    mat1_colsize = mat1->nrow*elemsize;
    mat2_colsize = mat2->nrow*elemsize;
    newmat_colsize = mat1_colsize + mat2_colsize;
    for (j = 0; j < newmat->ncol; j++) {
        offset = j * newmat_colsize;
        offset1 = j * mat1_colsize;
        offset2 = j * mat2_colsize;
        memcpy( (void *)newmat->data + offset, 
                (void *)mat1->data + offset1, mat1_colsize );
        memcpy( (void *)newmat->data + offset + mat1_colsize,
                       (void *)mat2->data + offset2, mat2_colsize );
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   newmat);
}




VALUE ratlas_map_bang(VALUE self, VALUE stor, VALUE block)
{
    VALUE val;
    Ratlas_Matrix *mat;
    unsigned long int i, nelem;
    double *r;
    Ratlas_Complex *c;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat);

    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < nelem; i++) {
                val = rb_float_new(r[i]);
                val = rb_funcall2(block, id_call, 1, &val);
                r[i] = NUM2DBL(val);
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                val = ratlas_rb_complex_new(c[i].r, c[i].i);
                val = rb_funcall2(block, id_call, 1, &val);
                c[i].r = NUM2DBL(rb_ivar_get(val, id_atre));
                c[i].i = NUM2DBL(rb_ivar_get(val, id_atim));
            }
            break;
    }
    return stor;
}




VALUE ratlas_each(VALUE self, VALUE stor, VALUE block)
{
    VALUE val;
    Ratlas_Matrix *mat;
    unsigned long int i, nelem;
    int m, n;
    double *r;
    Ratlas_Complex *c, ci;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat);

    m = mat->nrow;
    n = mat->ncol;
    nelem = m*n;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < nelem; i++) {
                val = rb_float_new(r[ratlas_swap_major(m,n,i)]);
                rb_funcall2(block, id_call, 1, &val);
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                ci = c[ratlas_swap_major(m,n,i)];
                val = ratlas_rb_complex_new(ci.r, ci.i);
                rb_funcall2(block, id_call, 1, &val);
            }
            break;
    }
    return stor;
}



VALUE ratlas_each_with_index(VALUE self, VALUE stor, VALUE block)
{
    VALUE val, idx, ary;
    Ratlas_Matrix *mat;
    unsigned long int i, nelem;
    int m, n;
    double *r;
    Ratlas_Complex *c, ci;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat);

    m = mat->nrow;
    n = mat->ncol;
    nelem = m*n;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < nelem; i++) {
                val = rb_float_new(r[ratlas_swap_major(m,n,i)]);
                idx = INT2FIX(i);
                ary = rb_ary_new3(2, val, idx);
                rb_funcall2(block, id_call, 1, &ary);
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                ci = c[ratlas_swap_major(m,n,i)];
                val = ratlas_rb_complex_new(ci.r, ci.i);
                idx = INT2FIX(i);
                ary = rb_ary_new3(2, val, idx);
                rb_funcall2(block, id_call, 1, &ary);
            }
            break;
    }
    return stor;
}




VALUE ratlas_each_with_ijindex(VALUE self, VALUE stor, VALUE block)
{
    VALUE ary;
    Ratlas_Matrix *mat;
    unsigned long int nelem;
    int i, j, m, n;
    double *r;
    Ratlas_Complex *c, ci;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat);

    m = mat->nrow;
    n = mat->ncol;
    nelem = m*n;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < m; i++)
              for (j = 0; j < n; j++) {
                ary = rb_ary_new3(3, rb_float_new(r[j*m+i]), 
                    INT2FIX(i), INT2FIX(j));
                rb_funcall2(block, id_call, 1, &ary);
              }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < m; i++)
              for (j = 0; j < n; j++) {
                ci = c[j*m+i];
                
                ary = rb_ary_new3(3, ratlas_rb_complex_new(ci.r, ci.i),
                   INT2FIX(i), INT2FIX(j));
                rb_funcall2(block, id_call, 1, &ary);
              }
            break;
    }
    return stor;
}




VALUE ratlas_zip_bang(VALUE self, VALUE stor, VALUE stors, VALUE block)
{
    VALUE *tmparr, val;
    Ratlas_Matrix *mat0, **mats;
    unsigned long int i, nelem;
/*     struct RArray *argarr; */
    int iarg, storslen;
    Ratlas_Complex *ctmp, *cmat0;
    
    Data_Get_Struct(stor, Ratlas_Matrix, mat0);
    nelem = mat0->nrow * mat0->ncol;
/*     argarr = RARRAY(stors); */
    storslen = RARRAY_LEN(stors);
    mats = ALLOCA_N(Ratlas_Matrix*, storslen);

    for (i = 0; i < storslen; i++)
    {
        Data_Get_Struct(RARRAY_PTR(stors)[i], Ratlas_Matrix, mats[i]);
        if ( (mats[i]->nrow != mat0->nrow) || 
                (mats[i]->ncol != mat0->ncol) || 
                (mats[i]->type != mat0->type) || 
                (mats[i]->matrixtype != mat0->matrixtype) )
            rb_raise(rb_eArgError, 
                    "Matrix not of same size or type.");
    }
    
    tmparr = ALLOCA_N(VALUE, storslen + 1);

    switch (mat0->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
            {
                tmparr[0] = rb_float_new( ((double *)mat0->data)[i] );
                for (iarg = 0; iarg < storslen; iarg++)
                {
                    tmparr[iarg+1] = rb_float_new(
                            ((double *)mats[iarg]->data)[i] );    
                }
            val = rb_funcall2(block, id_call, storslen + 1, tmparr);
            ( (double *)mat0->data )[i] = NUM2DBL(val);
            }
            break;
        case RATLAS_DCOMPLEX:
            cmat0 = (Ratlas_Complex *) mat0->data;
            for (i = 0; i < nelem; i++) {
                tmparr[0] = ratlas_rb_complex_new(cmat0[i].r, cmat0[i].i);
                for (iarg = 0; iarg < storslen; iarg++) {
                    ctmp = (Ratlas_Complex *) mats[iarg]->data;
                    tmparr[iarg+1] = ratlas_rb_complex_new(
                        ctmp[i].r, ctmp[i].i);    
                }
                val = rb_funcall2(block, id_call, storslen + 1, tmparr);
                cmat0[i].r = NUM2DBL(rb_ivar_get(val, id_atre));
                cmat0[i].i = NUM2DBL(rb_ivar_get(val, id_atim));
            }
            break;
    }
    return self;
}




/**
 * Clone a matrix and all its contents */
Ratlas_Matrix* ratlas_matrix_clone(Ratlas_Matrix *orig)
{
    Ratlas_Matrix *copy;

    copy = ratlas_matrix_alloc(orig->type, orig->matrixtype, orig->nrow,
            orig->ncol, orig->data);
    return copy;
}



/**
 * Duplicate the matrix container.  Don't care about copying contents. */
Ratlas_Matrix* ratlas_matrix_dup(Ratlas_Matrix *orig)
{
    Ratlas_Matrix *copy;

    copy = ratlas_matrix_alloc(orig->type, orig->matrixtype, orig->nrow,
            orig->ncol, NULL);
    return copy;
}



/**
 * Clone a storage object and its contents. */
VALUE ratlas_storage_clone(VALUE self, VALUE oldstore)
{
    VALUE newstore;
    Ratlas_Matrix *orig, *copy;

    Data_Get_Struct(oldstore, Ratlas_Matrix, orig);
    copy = ratlas_matrix_clone(orig);

    newstore = Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            copy);
    return newstore;
}




/**
 * Duplicate the storage object and the matrix container, but don't care
 * about the matrix contents */
VALUE ratlas_storage_dup(VALUE self, VALUE oldstore)
{
    VALUE newstore;
    Ratlas_Matrix *orig, *copy;

    Data_Get_Struct(oldstore, Ratlas_Matrix, orig);
    copy = ratlas_matrix_dup(orig);

    newstore = Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            copy);
    return newstore;
}




VALUE ratlas_storage_alloc(VALUE self, VALUE arg)
{
    Ratlas_Matrix *mat;
/*     struct RArray *carg; */
    int size1, size2, arglen;

    if (TYPE(arg) == T_FIXNUM) {
        size1 = FIX2INT(arg);
        size2 = 1;
    } else {
        Check_Type(arg, T_ARRAY);
	arglen = RARRAY_LEN(arg);
/*         carg = RARRAY(arg); */
        if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
        Check_Type(RARRAY_PTR(arg)[0], T_FIXNUM);
        size1 = FIX2INT(RARRAY_PTR(arg)[0]);
        if (arglen == 2) {
            Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
            size2 = FIX2INT(RARRAY_PTR(arg)[1]);
        } else {
            size2 = 1;
        }
    }    
    mat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, size1, size2, NULL);

    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free, mat);
}




VALUE ratlas_complex_storage_alloc(VALUE self, VALUE arg)
{
    Ratlas_Matrix *mat;
/*     struct RArray *carg; */
    int size1, size2, arglen;

    if (TYPE(arg) == T_FIXNUM) {
        size1 = FIX2INT(arg);
        size2 = 1;
    } else {
        Check_Type(arg, T_ARRAY);
	arglen = RARRAY_LEN(arg);
/*         carg = RARRAY(arg); */
        if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
        Check_Type(RARRAY_PTR(arg)[0], T_FIXNUM);
        size1 = FIX2INT(RARRAY_PTR(arg)[0]);
        if (arglen == 2) {
            Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
            size2 = FIX2INT(RARRAY_PTR(arg)[1]);
        } else {
            size2 = 1;
        }
    }    
    mat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE, 
            size1, size2, NULL);

    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            mat);
}




VALUE ratlas_re(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat, *newmat;
    unsigned long int nelem, i;
    VALUE retval=Qnil;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;

    switch (mat->type) {
        case RATLAS_DFLOAT:
            retval = ratlas_storage_clone(self, stor);
            break;
        case RATLAS_DCOMPLEX:
            newmat = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                    RATLAS_DENSE, mat->nrow, mat->ncol,
                    NULL);
            for (i = 0; i < nelem; i++) {
                ((double *)newmat->data)[i] = 
                    ((Ratlas_Complex *)mat->data)[i].r;
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0,
                ratlas_matrix_free, newmat);
            break;
    }
    return retval;
}




VALUE ratlas_im(VALUE self, VALUE stor)
{
    Ratlas_Matrix *mat, *newmat;
    unsigned long int nelem, i;
    VALUE newstor;
    double *zeros;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    zeros = xcalloc(nelem, ratlas_sizeof[RATLAS_DFLOAT]);
    newmat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE, mat->nrow, 
            mat->ncol, (void *) zeros);
    free(zeros);
    newstor = Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            newmat);

    switch (mat->type) {
        case RATLAS_DFLOAT:
            break;
        case RATLAS_DCOMPLEX:
            for (i = 0; i < nelem; i++) {
                ((double *)newmat->data)[i] = 
                    ((Ratlas_Complex *)mat->data)[i].i;
            }
            break;
    }
    return newstor;
}




VALUE ratlas_reshape_bang(VALUE self, VALUE stor, VALUE arg)
{
    Ratlas_Matrix *mat;
    unsigned long int nelem;
    int size1, size2, arglen;
/*     struct RArray *carg; */

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;

    if (TYPE(arg) == T_FIXNUM) {
        size1 = FIX2INT(arg);
        if (size1 != nelem)
            rb_raise(rb_eArgError, "Illegal reshape.");
        size2 = 1;
    } else
    {
        Check_Type(arg, T_ARRAY);
/*         carg = RARRAY(arg); */
	arglen = RARRAY_LEN(arg);
        if (arglen > 2) rb_raise(rb_eArgError, "Unexpected dimensions.");
        Check_Type(RARRAY_PTR(arg)[0], T_FIXNUM);
        size1 = NUM2INT(RARRAY_PTR(arg)[0]);
        if (arglen == 2) {
            Check_Type(RARRAY_PTR(arg)[1], T_FIXNUM);
            size2 = NUM2INT(RARRAY_PTR(arg)[1]);
        } else {
            size2 = 1;
        }
    }    
    
    if (nelem != size1*size2) rb_raise(rb_eArgError, "Illegal reshape.");
    mat->nrow = size1;
    mat->ncol = size2;
    return stor;
}




/* Ratlas operations */

VALUE ratlas_mul_bang(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat;
    double a, ca[2];
    unsigned long int nelem;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            a = NUM2DBL(alpha);    
            cblas_dscal(nelem, a, (double *)mat->data, 1); 
            break;
        case RATLAS_DCOMPLEX:
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, "Expected complex.");
            ca[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            ca[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zscal(nelem, ca, (Ratlas_Complex *)mat->data, 1); 
            break;
    }
    return stor;
}



VALUE ratlas_mul(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat, *retmat;
    double a, ca[2];
    unsigned long int nelem;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_clone (mat); 
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            a = NUM2DBL(alpha);    
            cblas_dscal(nelem, a, (double *)retmat->data, 1); 
            break;
        case RATLAS_DCOMPLEX:
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, "Expected complex.");
            ca[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            ca[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cblas_zscal(nelem, ca, (Ratlas_Complex *)retmat->data, 1); 
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




VALUE ratlas_div_bang(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat;
    double a, ca[2];
    unsigned long int nelem;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            a = NUM2DBL(alpha);
            if (a == 0.0)
                rb_raise(rb_eZeroDivError, "Division by zero.");
            a = 1/a;
            cblas_dscal(nelem, a, (double *)mat->data, 1); 
            break;
        case RATLAS_DCOMPLEX:
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, "Expected complex.");
            ca[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            ca[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            if (ca[0]*ca[0] + ca[1]*ca[1] == 0.0)
                rb_raise(rb_eZeroDivError, "Division by zero.");
            ca[0] = 1.0/ca[0];
            ca[1] = 1.0/ca[1];
            cblas_zscal(nelem, ca, (Ratlas_Complex *)mat->data, 1); 
            break;
    }
    return stor;
}




VALUE ratlas_div(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat, *retmat;
    double a, ca[2];
    unsigned long int nelem;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_clone(mat); 
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            a = NUM2DBL(alpha);
            if (a == 0.0)
                rb_raise(rb_eZeroDivError, "Division by zero.");
            a = 1/a;
            cblas_dscal(nelem, a, (double *)retmat->data, 1); 
            break;
        case RATLAS_DCOMPLEX:
            if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
                rb_raise(rb_eArgError, "Expected complex.");
            ca[0] = NUM2DBL(rb_ivar_get(alpha, id_atre));
            ca[1] = NUM2DBL(rb_ivar_get(alpha, id_atim));
            if (ca[0]*ca[0] + ca[1]*ca[1] == 0.0)
                rb_raise(rb_eZeroDivError, "Division by zero.");
            ca[0] = 1.0/ca[0];
            ca[1] = 1.0/ca[1];
            cblas_zscal(nelem, ca, (Ratlas_Complex *)retmat->data, 1); 
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




VALUE ratlas_add_bang(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat;
    double c, d;
    unsigned long int nelem, i;
    Ratlas_Complex *cdest;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            c = NUM2DBL(alpha);
            for (i = 0; i < nelem; i++) {
                ((double *)mat->data)[i] += c;
            }
            break;
        case RATLAS_DCOMPLEX:
            //if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
            //    rb_raise(rb_eArgError, "Expected complex.");
            c = NUM2DBL(rb_ivar_get(alpha, id_atre));
            d = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cdest = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                cdest[i].r += c;
                cdest[i].i += d;
            }
            break;
    }
    return stor;
}




VALUE ratlas_add(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat, *retmat;
    double c, d, *rsrc, *rdest;
    unsigned long int nelem, i;
    Ratlas_Complex *csrc, *cdest;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_alloc (mat->type, mat->matrixtype, mat->nrow, 
            mat->ncol, NULL); 
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            c = NUM2DBL(alpha);
            rsrc = (double *) mat->data;
            rdest = (double *) retmat->data;
            for (i = 0; i < nelem; i++) {
                rdest[i] = rsrc[i] + c;
            }
            break;
        case RATLAS_DCOMPLEX:
            //if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
            //    rb_raise(rb_eArgError, "Expected complex.");
            c = NUM2DBL(rb_ivar_get(alpha, id_atre));
            d = NUM2DBL(rb_ivar_get(alpha, id_atim));
            csrc = (Ratlas_Complex *) mat->data;
            cdest = (Ratlas_Complex *) retmat->data;
            for (i = 0; i < nelem; i++) {
                cdest[i].r = csrc[i].r + c;
                cdest[i].i = csrc[i].i + d;
            }
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}



VALUE ratlas_sub_bang(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat;
    double c, d;
    unsigned long int nelem, i;
    Ratlas_Complex *cdest;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            c = NUM2DBL(alpha);
            for (i = 0; i < nelem; i++)
                ((double *)mat->data)[i] -= c;
            break;
        case RATLAS_DCOMPLEX:
            //if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
            //    rb_raise(rb_eArgError, "Expected complex.");
            c = NUM2DBL(rb_ivar_get(alpha, id_atre));
            d = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cdest = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                cdest[i].r -= c;
                cdest[i].i -= d;
            }
            break;
    }
    return stor;
}




VALUE ratlas_sub(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat, *retmat;
    double c, d;
    unsigned long int nelem, i;
    Ratlas_Complex *cdest, *csrc;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    retmat = ratlas_matrix_alloc(mat->type, mat->matrixtype, mat->nrow, 
            mat->ncol, NULL); 
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            c = NUM2DBL(alpha);
            for (i = 0; i < nelem; i++)
                ((double *)retmat->data)[i] = ((double *)mat->data)[i] - c;
            break;
        case RATLAS_DCOMPLEX:
            //if (!ratlas_check_type(alpha, RATLAS_DCOMPLEX))
            //    rb_raise(rb_eArgError, "Expected complex.");
            c = NUM2DBL(rb_ivar_get(alpha, id_atre));
            d = NUM2DBL(rb_ivar_get(alpha, id_atim));
            cdest = (Ratlas_Complex *) retmat->data;
            csrc = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++)
            {
                cdest[i].r = csrc[i].r - c;
                cdest[i].i = csrc[i].i - d;
            }
            break;
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




VALUE ratlas_rem_bang(VALUE self, VALUE stor, VALUE denom)
{
    Ratlas_Matrix *mat;
    unsigned long int nelem, i;
    int cdenom, cnum;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    if (mat->type != RATLAS_DFLOAT)
        rb_raise(rb_eArgError, "Reminder only for real matrices.");
    Check_Type(denom, T_FIXNUM);
    
    nelem = mat->nrow * mat->ncol;
    cdenom = FIX2INT(denom);
    
    for (i = 0; i < nelem; i++)
    {
        cnum = ((double *)mat->data)[i];
        ((double *)mat->data)[i] = (double) (cnum % cdenom);
    }
    return stor;
}




static void apac(Ratlas_Matrix *amat, double *alpha, Ratlas_Matrix *cmat)
{
    int j;
    double *ra, *rc;
    Ratlas_Complex *ca, *cc;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (j = 0; j < amat->ncol; j++) {
                cblas_daxpy(amat->nrow, alpha[0], rc + j*cmat->nrow,
                        1, ra + j*amat->nrow, 1);
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (j = 0; j < amat->ncol; j++) {
                cblas_zaxpy(amat->nrow, alpha, cc + j*cmat->nrow,
                        1, ca + j*amat->nrow, 1);
            }
            break;
    }
}

/**
 * A <- A + alpha*C' */
static void apact(Ratlas_Matrix *amat, void *alpha, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    double a, b, c, d, *ra, *rc;
    Ratlas_Complex *ca, *cc, calpha;    

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            a = *((double *)alpha);
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ra[i] += a*rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            calpha = *((Ratlas_Complex *) alpha);
            a = calpha.r;
            b = calpha.i;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                c = cc[idx].r;
                d = cc[idx].i;
                ca[i].r += (a*c - b*d);
                ca[i].i += (a*d + b*c);
            }
            break;
    }
}





/**
 * A' <- A' + alpha*C */
static void atpac(Ratlas_Matrix *amat, void *alpha, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    double a, b, c, d, *ra, *rc;
    Ratlas_Complex *ca, *cc, calpha;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            a = *((double *)alpha);
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ra[idx] += a*rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            calpha = *((Ratlas_Complex *)alpha);
            a = calpha.r;
            b = calpha.i;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                c = cc[i].r;
                d = cc[i].i;
                ca[idx].r += (a*c - b*d);
                ca[idx].i += (a*d + b*c);
            }
            break;
    }
}




/**
 * A <- A + C */
static void apc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, i=0;
    Ratlas_Complex *ca, *cc;

    nelem = amat->nrow * amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)amat->data)[i] += ((double *)cmat->data)[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                ca[i].r += cc[i].r;
                ca[i].i += cc[i].i;
            }
            break;
    }
}




/**
 * A' <- A' + C */
static void atpc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i, ire;
    int m, n;
    Ratlas_Complex *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[idx] += ((double *)cmat->data)[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                ire = i;
                ca[idx].r += cc[ire].r;
                ca[idx].i += cc[ire].i;
            }
            break;
    }
}




/**
 * A <- A + C' */
static void apct_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i, ire;
    int m, n;
    Ratlas_Complex *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[i] += ((double *)cmat->data)[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                ire = i;
                ca[ire].r += cc[idx].r;
                ca[ire].i += cc[idx].i;
            }
            break;
    }
}




/*
 * op(A) <- op(A)+op(C) */
VALUE ratlas_madd_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{

    Ratlas_Matrix *amat, *cmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    /* Check dimentions.  Identify if both or none have been transposed.*/
    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        apc_bang(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa > CblasNoTrans)
            atpc_bang(amat, cmat);
        else
            apct_bang(amat, cmat);
    }
    return astor;
}




/**
 * returns A + C */
static Ratlas_Matrix* apc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, i=0;
    Ratlas_Matrix *retmat;
    double *ra, *rc, *rr;
    Ratlas_Complex *ca, *cc, *cr;

    nelem = amat->nrow * amat->ncol;
    retmat = ratlas_matrix_alloc (amat->type, amat->matrixtype, amat->nrow, 
            amat->ncol, NULL);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            rr = (double *) retmat->data;
            for (i = 0; i < nelem; i++)
                rr[i] = ra[i] + rc[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            cr = (Ratlas_Complex *) retmat->data;
            for (i = 0; i < nelem; i++) {
                cr[i].r = ca[i].r + cc[i].r;
                cr[i].i = ca[i].i + cc[i].i;
            }
            break;
    }
    return retmat;
}




/**
 * returns A' + C */
static Ratlas_Matrix* atpc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i, ire;
    int m, n;
    Ratlas_Matrix *retmat;
    double *ra, *rr, *rc;
    Ratlas_Complex *ca, *cr, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    retmat = ratlas_matrix_alloc (cmat->type, cmat->matrixtype, cmat->nrow, 
            cmat->ncol, NULL); 
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            rr = (double *) retmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[idx] + rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            cr = (Ratlas_Complex *) retmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                ire = i;
                cr[i].r = ca[idx].r + cc[i].r;
                cr[i].i = ca[idx].i + cc[i].i;
            }
            break;
    }
    return retmat;
}




/**
 * returns A + C' */
static Ratlas_Matrix* apct(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    Ratlas_Matrix *retmat;
    double *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;

    retmat = ratlas_matrix_alloc (amat->type, amat->matrixtype, amat->nrow, 
            amat->ncol, NULL); 
        
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[i] + rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                cr[i].r = ca[i].r + cc[idx].r;
                cr[i].i = ca[i].i + cc[idx].i;
            }
            break;
    }
    return retmat;
}




/*
 * returns op(A)+op(C) */
VALUE ratlas_madd(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{

    Ratlas_Matrix *amat, *cmat, *retmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    /* Check dimentions.  Identify if both or none have been transposed.*/
    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");

        retmat = apc(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa > CblasNoTrans)
            retmat = atpc(amat, cmat);
        else
            retmat = apct(amat, cmat);
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




static void amc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, i=0;
    Ratlas_Complex *ca, *cc;

    nelem = amat->nrow * amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)amat->data)[i] -= ((double *)cmat->data)[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                ca[i].r -= cc[i].r;
                ca[i].i -= cc[i].i;
            }
            break;
    }
}




static void atmc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    Ratlas_Complex *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[idx] -= ((double *)cmat->data)[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                ca[idx].r -= cc[i].r;
                ca[idx].i -= cc[i].i;
            }
            break;
    }
}




static void amct_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    Ratlas_Complex *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[i] -= ((double *)cmat->data)[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n);
                ca[i].r -= cc[idx].r;
                ca[i].i -= cc[idx].i;
            }
            break;
    }
}





/*
 * op(A) <- op(A) - op(C) */
VALUE ratlas_msub_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    /* Check dimentions.  Identify if both or none have been transposed.*/
    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        amc_bang(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa == CblasTrans)
            atmc_bang(amat, cmat);
        else
            amct_bang(amat, cmat);
    }
    return astor;
}



/**
 * returns A - C */
static Ratlas_Matrix* amc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, i=0;
    Ratlas_Matrix *retmat;
    double *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;

    nelem = amat->nrow * amat->ncol;
    retmat = ratlas_matrix_alloc (amat->type, amat->matrixtype, amat->nrow, 
            amat->ncol, NULL);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++)
                rr[i] = ra[i] - rc[i];
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                cr[i].r = ca[i].r - cc[i].r;
                cr[i].i = ca[i].i - cc[i].i;
            }
            break;
    }
    return retmat;
}




/**
 * returns A' - C */
static Ratlas_Matrix* atmc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    Ratlas_Matrix *retmat;
    double *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;
    
    retmat = ratlas_matrix_alloc (cmat->type, cmat->matrixtype, cmat->nrow, 
            cmat->ncol, NULL); 
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[idx] - rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                cr[i].r = ca[idx].r - cc[i].r;
                cr[i].i = ca[idx].i - cc[i].i;
            }
            break;
    }
    return retmat;
}




/**
 * returns A - C' */
static Ratlas_Matrix* amct(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int nelem, idx, i;
    int m, n;
    Ratlas_Matrix *retmat;
    double *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;

    m = amat->nrow;
    n = amat->ncol;
    nelem = m*n;

    retmat = ratlas_matrix_alloc (amat->type, amat->matrixtype, amat->nrow, 
            amat->ncol, NULL); 
        
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[i] - rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                cr[i].r = ca[i].r - cc[idx].r;
                cr[i].i = ca[i].i - cc[idx].i;
            }
            break;
    }
    return retmat;
}




/*
 * returns op(A)-op(C) */
VALUE ratlas_msub(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{

    Ratlas_Matrix *amat, *cmat, *retmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    /* Check dimentions.  Identify if both or none have been transposed.*/
    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");

        retmat = amc(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa > CblasNoTrans)
            retmat = atmc(amat, cmat);
        else
            retmat = amct(amat, cmat);
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




/**
 * A <- A.*C */
static void mul_ac_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem;
    double a, b, c, d;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)amat->data)[i] *= ((double *)cmat->data)[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                a = ca[i].r;
                b = ca[i].i;
                c = cc[i].r;
                d = cc[i].i;
                ca[i].r = a*c - b*d;
                ca[i].i = a*d + b*c;
            }
            break;
    }
}



/**
 * A' <- A'.*C */
static void mul_atc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[idx] *= ((double *)cmat->data)[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n); 
                a = ca[idx].r;
                b = ca[idx].i;
                c = cc[i].r;
                d = cc[i].i;
                ca[idx].r = a*c - b*d;
                ca[idx].i = a*d + b*c;
            }
            break;
    }    
}




/**
 * A <- A.*C' */
static void mul_act_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ((double *)amat->data)[i] *= ((double *)cmat->data)[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = (i%n * m + i/n); 
                a = ca[i].r;
                b = ca[i].i;
                c = cc[idx].r;
                d = cc[idx].i;
                ca[i].i = a*c - b*d;
                ca[i].r = a*d + b*c;
            }
            break;
    }    
}




/**
 * op(A) <- op(A).*op(C) */
VALUE ratlas_mmul_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        mul_ac_bang(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa == CblasTrans)
            mul_atc_bang(amat, cmat);
        else
            mul_act_bang(amat, cmat);
    }
    return astor;
}



/**
 * returns A.*C */
static Ratlas_Matrix* mul_ac(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    retmat = ratlas_matrix_dup(amat);    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++)
                rr[i] = ra[i] * rc[i];
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                a = ca[i].r;
                b = ca[i].i;
                c = cc[i].r;
                d = cc[i].i;
                cr[i].r = a*c - b*d;
                cr[i].i = a*d + b*c;
            }
            break;
    }
    return retmat;
}



/**
 * returns A'.*C */
static Ratlas_Matrix* mul_atc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;
    
    retmat = ratlas_matrix_dup(cmat);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[idx] * rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[idx].r;
                b = ca[idx].i;
                c = cc[i].r;
                d = cc[i].i;
                cr[idx].r = a*c - b*d;
                cr[idx].i = a*d + b*c;
            }
            break;
    }
    return retmat;
}




/**
 * returns A.*C' */
static Ratlas_Matrix* mul_act(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;

    retmat = ratlas_matrix_dup(amat);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[i] * rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[i].r;
                b = ca[i].i;
                c = cc[idx].r;
                d = cc[idx].i;
                cr[i].r = a*c - b*d;
                cr[i].i = a*d + b*c;
            }
            break;
    }
    return retmat;    
}




/**
 * returns op(A).*op(C) */
VALUE ratlas_mmul(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat, *retmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        retmat = mul_ac(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa == CblasTrans)
            retmat = mul_atc(amat, cmat);
        else
            retmat = mul_act(amat, cmat);
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




/**
 * A <- A./C */
static void div_ac_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem;
    double a, b, c, d, *ra, *rc;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++)
                ra[i] = ra[i]/rc[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                a = ca[i].r;
                b = ca[i].i;
                c = cc[i].r;
                d = cc[i].i;
                ca[i].r = (a*c + b*d)/(c*c + d*d);
                ca[i].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }
}



/**
 * A' <- A'./C */
static void div_atc_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *ra, *rc;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ra[idx] = ra[idx]/rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[idx].r;
                b = ca[idx].i;
                c = cc[i].r;
                d = cc[i].i;
                ca[idx].r = (a*c + b*d)/(c*c + d*d);
                ca[idx].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }    
}




/**
 * A <- A./C' */
static void div_act_bang(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *ra, *rc;
    Ratlas_Complex *ca, *cc;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;

    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                ra[i] = ra[i]/rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[i].r;
                b = ca[i].i;
                c = cc[idx].r;
                d = cc[idx].i;
                ca[i].r = (a*c + b*d)/(c*c + d*d);
                ca[i].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }    
}




/**
 * op(A) <- op(A)./op(C) */
VALUE ratlas_mdiv_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        div_ac_bang(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa > CblasNoTrans)
            div_atc_bang(amat, cmat);
        else
            div_act_bang(amat, cmat);
    }
    return astor;
}




/**
 *  A./C */
static Ratlas_Matrix* div_ac(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    retmat = ratlas_matrix_dup(amat);    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++)
                rr[i] = ra[i]/rc[i];
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                a = ca[i].r;
                b = ca[i].i;
                c = cc[i].r;
                d = cc[i].i;
                cr[i].r = (a*c + b*d)/(c*c + d*d);
                cr[i].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }
    return retmat;
}



/**
 *  A'./C */
static Ratlas_Matrix* div_atc(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;
    retmat = ratlas_matrix_dup(cmat);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[idx] = ra[idx]/rc[i];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[idx].r;
                b = ca[idx].i;
                c = cc[i].r;
                d = cc[i].i;
                cr[idx].r = (a*c + b*d)/(c*c + d*d);
                cr[idx].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }
    return retmat;    
}




/**
 *  A./C' */
static Ratlas_Matrix* div_act(Ratlas_Matrix *amat, Ratlas_Matrix *cmat)
{
    unsigned long int i, nelem, idx;
    int m, n;
    double a, b, c, d, *rr, *ra, *rc;
    Ratlas_Complex *cr, *ca, *cc;
    Ratlas_Matrix *retmat;

    nelem = amat->ncol * amat->nrow;
    m = amat->nrow;
    n = amat->ncol;
    retmat = ratlas_matrix_dup(amat);
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rr = (double *) retmat->data;
            ra = (double *) amat->data;
            rc = (double *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                rr[i] = ra[i]/rc[idx];
            }
            break;
        case RATLAS_DCOMPLEX:
            cr = (Ratlas_Complex *) retmat->data;
            ca = (Ratlas_Complex *) amat->data;
            cc = (Ratlas_Complex *) cmat->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n; 
                a = ca[i].r;
                b = ca[i].i;
                c = cc[idx].r;
                d = cc[idx].i;
                cr[i].r = (a*c + b*d)/(c*c + d*d);
                cr[i].i = (b*c - a*d)/(c*c + d*d);
            }
            break;
    }
    return retmat;    
}




/**
 * op(A)./op(C) */
VALUE ratlas_mdiv(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor)
{
    Ratlas_Matrix *amat, *cmat, *retmat;
    int ctransa, ctransc;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(cstor, Ratlas_Matrix, cmat);
    ctransa = FIX2INT(transa);
    ctransc = FIX2INT(transc);

    if (amat->type != cmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    if ( ctransa+ctransc == 2*ctransa) {
        if ((amat->nrow != cmat->nrow) || (amat->ncol != cmat->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        retmat = div_ac(amat, cmat);
    } else {
        if ((amat->nrow != cmat->ncol) || (amat->ncol != cmat->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        if (ctransa > CblasNoTrans)
            retmat = div_atc(amat, cmat);
        else
            retmat = div_act(amat, cmat);
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
                   retmat);
}




VALUE ratlas_pow_bang(VALUE self, VALUE stor, VALUE alpha)
{
    Ratlas_Matrix *mat;
    double c;
    unsigned long int nelem, i;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            ratlas_assertnum(alpha);
            c = NUM2DBL(alpha);
            /* if (modf(c, ip) > 0.0)
            {
                rb_raise(rb_eArgError, 
                "pow() not implemented for complex results.");
            }
            */
            for (i = 0; i < nelem; i++)
                ((double*)mat->data)[i] = pow(((double *)mat->data)[i], c);
            break;
        case RATLAS_DCOMPLEX:
            rb_raise(rb_eArgError, 
                    "pow() not implemented for complex.");
            break;
    }
    return stor;
}




VALUE ratlas_sqrt_bang(VALUE self, VALUE stor)
{
    unsigned long int nelem, i;
    Ratlas_Matrix *mat;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)mat->data)[i] = sqrt(((double *)mat->data)[i]); 
            break;
        case RATLAS_DCOMPLEX:
            rb_raise(rb_eArgError, "sqrt() not implemented for complex.");
            break;
    }
    return stor;
}




VALUE ratlas_log_bang(VALUE self, VALUE stor)
{
    unsigned long int nelem, i;
    Ratlas_Matrix *mat;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)mat->data)[i] = log(((double *)mat->data)[i]); 
            break;
        case RATLAS_DCOMPLEX:
            rb_raise(rb_eArgError, 
                    "log() not implemented for complex.");
            break;
    }
    return stor;
}




VALUE ratlas_exp_bang(VALUE self, VALUE stor)
{
    unsigned long int nelem, i;
    Ratlas_Matrix *mat;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            for (i = 0; i < nelem; i++)
                ((double *)mat->data)[i] = exp(((double *)mat->data)[i]); 
            break;
        case RATLAS_DCOMPLEX:
            rb_raise(rb_eArgError, 
                    "exp() not implemented for complex.");
            break;
    }
    return stor;
}



/*
 * Sum all elements. */  
VALUE ratlas_sum(VALUE self, VALUE stor)
{
    unsigned long int nelem, i;
    Ratlas_Matrix *mat;
    double resum=0.0, imsum=0.0, *r;
    Ratlas_Complex *c;
    VALUE retval=Qnil;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    nelem = mat->nrow * mat->ncol;
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            for (i = 0; i < nelem; i++)
                resum += r[i];
            retval = rb_float_new(resum);
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            for (i = 0; i < nelem; i++) {
                resum += c[i].r;
                imsum += c[i].i;
            }
            retval = ratlas_rb_complex_new(resum, imsum);
            break;
    }
    return retval;
}




/*
 * Sum over rows: sum_i A(i,j). */
VALUE ratlas_rowsum(VALUE self, VALUE stor)
{
    unsigned long int idx;
    int i, j;
    Ratlas_Matrix *mat, *summat;
    double *zeros, *rs, *rm;
    Ratlas_Complex *cs, *cm;
    VALUE retval=Qnil;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            zeros = xcalloc(mat->ncol, sizeof(double));
            summat = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                    RATLAS_DENSE, 1, mat->ncol, zeros);
            free(zeros);
            rs = (double *) summat->data;
            rm = (double *) mat->data;
            for (j = 0; j < mat->ncol; j++) {
                for (i = 0; i < mat->nrow; i++) {
                    idx = j*mat->nrow + i;
                    rs[j] += rm[idx];
                }
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, summat);
            break;
        case RATLAS_DCOMPLEX:
            zeros = xcalloc(mat->ncol, ratlas_sizeof[RATLAS_DCOMPLEX]);
            summat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
                    RATLAS_DENSE, 1, mat->ncol, zeros);
            free(zeros);
            cs = (Ratlas_Complex *) summat->data;
            cm = (Ratlas_Complex *) mat->data;
            for (j = 0; j < mat->ncol; j++) {
                for (i = 0; i < mat->nrow; i++) {
                    idx = j*mat->nrow + i;
                    cs[j].r += cm[idx].r;
                    cs[j].i += cm[idx].i;
                }
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, summat);
            break;
    }
    return retval;
}




/*
 * Sum over columns: sum_j A(i,j). */
VALUE ratlas_colsum(VALUE self, VALUE stor)
{
    unsigned long int idx;
    int i, j;
    Ratlas_Matrix *mat, *summat;
    double *zeros, *rs, *rm;
    Ratlas_Complex *cs, *cm;
    VALUE retval=Qnil;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    switch (mat->type) {
        case RATLAS_DFLOAT:
            zeros = xcalloc(mat->nrow, sizeof(double));
            summat = ratlas_matrix_alloc(RATLAS_DFLOAT, 
                    RATLAS_DENSE, mat->nrow, 1, zeros);
            free(zeros);
            rs = (double *) summat->data;
            rm = (double *) mat->data;
            for (i = 0; i < mat->nrow; i++) {
                for (j = 0; j < mat->ncol; j++) {
                    idx = j*mat->nrow + i;
                    rs[i] += rm[idx];
                }
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, summat);
            break;
        case RATLAS_DCOMPLEX:
            zeros = xcalloc(mat->nrow, 2*sizeof(double));
            summat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, 
                    RATLAS_DENSE, mat->nrow, 1, zeros);
            free(zeros);
            cs = (Ratlas_Complex *) summat->data;
            cm = (Ratlas_Complex *) mat->data;
            for (i = 0; i < mat->nrow; i++) {
                for (j = 0; j < mat->ncol; j++) {
                    idx = j*mat->nrow + i;
                    cs[i].r += cm[idx].r;
                    cs[i].i += cm[idx].i;
                }
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, summat);
            break;
    }
    return retval;
} 



static VALUE colon_AtB(Ratlas_Matrix *mat1, Ratlas_Matrix *mat2)
{
    unsigned long int nelem, idx, i;
    int m, n;
    double a, b, c, d;
    double resum=0.0, imsum=0.0, *r1, *r2;
    Ratlas_Complex *c1, *c2;

    m =  mat1->ncol;
    n = mat1->nrow;
    nelem = mat1->ncol * mat1->nrow;

    switch (mat1->type) {
        case RATLAS_DFLOAT:
            r1 = (double *) mat1->data;
            r2 = (double *) mat2->data;
            for (i = 0; i < nelem; i++) {
                idx =  i%n * m + i/n;
                resum += r1[idx] * r2[i];
            }
            return rb_float_new(resum);
            break;
        case RATLAS_DCOMPLEX:
            c1 = (Ratlas_Complex *) mat1->data;
            c2 = (Ratlas_Complex *) mat2->data;
            for (i = 0; i < nelem; i++) {
                idx = i%n * m + i/n;
                a = c1[idx].r;
                b = c1[idx].i;
                c = c2[i].r;
                d = c2[i].i;
                resum += a*c - b*d;
                imsum += a*d + b*c;
            }
            return ratlas_rb_complex_new(resum, imsum);
            break;
        default:
            rb_raise(rb_eArgError, "Unsupported type in matrix.");
            break;
    }
}




static VALUE colon_AB(Ratlas_Matrix *mat1, Ratlas_Matrix *mat2)
{
    unsigned long int nelem, i;
    double a, b, c, d;
    double resum=0.0, imsum=0.0, *r1, *r2;
    Ratlas_Complex *c1, *c2;
    
    nelem = mat1->nrow * mat1->ncol;
    switch (mat1->type) {
        case RATLAS_DFLOAT:
            r1 = (double *) mat1->data;
            r2 = (double *) mat2->data;
            for (i = 0; i < nelem; i++)
                resum += r1[i] * r2[i];
            return rb_float_new(resum);
            break;
        case RATLAS_DCOMPLEX:
            c1 = (Ratlas_Complex *) mat1->data;
            c2 = (Ratlas_Complex *) mat2->data;
            for (i = 0; i < nelem; i++) {
                a = c1[i].r;
                b = c1[i].i;
                c = c2[i].r;
                d = c2[i].i;
                resum += a*c - b*d;
                imsum += a*d + b*c;
            }
            return ratlas_rb_complex_new(resum, imsum);
            break;
        default:
            rb_raise(rb_eArgError, "Unsupported type in matrix.");
            break;
    }
}




/**
 * s <-- sum_i sum_j op(A)ij * Bij */
VALUE ratlas_colon(VALUE self, VALUE trans1, VALUE stor1, VALUE trans2, 
        VALUE stor2)
{
    Ratlas_Matrix *mat1, *mat2;
    int ctrans1, ctrans2;
    
    Data_Get_Struct(stor1, Ratlas_Matrix, mat1);
    Data_Get_Struct(stor2, Ratlas_Matrix, mat2);
    ctrans1 = FIX2INT(trans1);
    ctrans2 = FIX2INT(trans2);
    
    if (mat1->type != mat2->type)
        rb_raise(rb_eArgError, "Matrices must be of same type.");

    /* Check dimentions.  Identify if both or none have been transposed.*/
    if ( ctrans1 != ctrans2) {
        if ((mat1->nrow != mat2->ncol) || (mat1->ncol != mat2->nrow))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        return colon_AtB(mat1, mat2);
    } else {
        if ((mat1->nrow != mat2->nrow) || (mat1->ncol != mat2->ncol))
            rb_raise(rb_eArgError, 
                    "Matrices must have same dimensions.");
        return colon_AB(mat1, mat2);
    }
}




/*
 * A <- D*A  where D is a diagonal matrix, stored as a column vector. */
VALUE ratlas_damul_bang(VALUE self, VALUE astor, VALUE dstor)
{
    Ratlas_Matrix *amat, *dmat;
    int i, j;
    unsigned long int cursor;
    double a, b, c, d, *rd, *ra;
    Ratlas_Complex *cd, *ca;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(dstor, Ratlas_Matrix, dmat);

    if (dmat->nrow != amat->nrow)
        rb_raise(rb_eArgError, "Dimension mismatch.");
    if (amat->type != dmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");
    
    switch (amat->type) {
        case RATLAS_DFLOAT:
            rd = (double *) dmat->data;
            ra = (double *) amat->data;
            for (i = 0; i < dmat->nrow; i++) {
                d = rd[i];
                for (j = 0; j < amat->ncol; j++) {
                    ra[i+j*amat->nrow] *= d;
                }
            }
            break;
        case RATLAS_DCOMPLEX:
            cd = (Ratlas_Complex *) dmat->data;
            ca = (Ratlas_Complex *) amat->data;
            for (i = 0; i < dmat->nrow; i++) {
                c = cd[i].r;
                d = cd[i].i;
                for (j = 0; j < amat->ncol; j++) {
                    cursor = i+j*amat->nrow;
                    a = ca[cursor].r;
                    b = ca[cursor].i;
                    ca[cursor].r = a*c - b*d;
                    ca[cursor].i = a*d + b*c;
                }
            }
        default:
            rb_raise(rb_eArgError, "Unknown type.");
    }
    return astor;
}




/*
 * A <- A*D  where D is a diagonal matrix, stored as a column vector. */
VALUE ratlas_admul_bang(VALUE self, VALUE astor, VALUE dstor)
{
    Ratlas_Matrix *amat, *dmat;
    int j, i, pos;
    double *ra, *rd, a, b, c, d;
    Ratlas_Complex *ca, *cd;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(dstor, Ratlas_Matrix, dmat);

    if (amat->ncol != dmat->nrow)
        rb_raise(rb_eArgError, "Dimension mismatch.");
    if (amat->type != dmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");

    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rd = (double *) dmat->data;
            for (j = 0; j < amat->ncol; j++)
                cblas_dscal(amat->nrow, rd[j], ra+j*amat->nrow, 1);
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cd = (Ratlas_Complex *) dmat->data;
            for (j = 0; j < amat->ncol; j++) {
                c = cd[j].r;
                d = cd[j].i;
                for (i = 0; i < amat->nrow; i++) {
                    pos = j*amat->nrow + i;
                    a = ca[pos].r;
                    b = ca[pos].i;
                    ca[pos].r = a*c - b*d;
                    ca[pos].i = a*d + b*c;
                }
            }
            break;
        default:
            rb_raise(rb_eArgError, "Unknown type.");
    }
    return astor;

}




/*
 * A <-- A+D, where D is a diagonal matrix with the diagonal stored as a vector
 * (1-column matrix) */
VALUE ratlas_dadd_bang(VALUE self, VALUE astor, VALUE dstor)
{
    Ratlas_Matrix *dmat, *amat;
    int i, nrow, ncol, cursor;
    double *ra, *rd;
    Ratlas_Complex *ca, *cd;

    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(dstor, Ratlas_Matrix, dmat);

    if (amat->ncol != dmat->nrow)
        rb_raise(rb_eArgError, "Dimension mismatch.");
    if (amat->type != dmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");
    
    nrow = dmat->nrow;
    ncol = amat->ncol;
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rd = (double *) dmat->data;
            for (i = 0; i < nrow; i++)
                ra[i*nrow + i] += rd[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cd = (Ratlas_Complex *) dmat->data;
            for (i = 0; i < nrow; i++) {
                cursor = i*nrow + i; 
                ca[cursor].r += cd[i].r;
                ca[cursor].i += cd[i].i;
            }
            break;
        default:
            rb_raise(rb_eArgError, "Unknown type.");
    }
    return astor;
}




/*
 * A <-- A-D, where D is a diagonal matrix with the diagonal stored as a vector
 * (1-column matrix) */
VALUE ratlas_dsub_bang(VALUE self, VALUE astor, VALUE dstor)
{
    Ratlas_Matrix *dmat, *amat;
    int i, nrow, ncol, cursor;
    double *ra, *rd;
    Ratlas_Complex *ca, *cd;
    
    Data_Get_Struct(astor, Ratlas_Matrix, amat);
    Data_Get_Struct(dstor, Ratlas_Matrix, dmat);

    if (amat->ncol != dmat->nrow)
        rb_raise(rb_eArgError, "Dimension mismatch.");
    if (amat->type != dmat->type)
        rb_raise(rb_eArgError, "Matrices must have same type.");
    
    nrow = dmat->nrow;
    ncol = amat->ncol;
    switch (amat->type) {
        case RATLAS_DFLOAT:
            ra = (double *) amat->data;
            rd = (double *) dmat->data;
            for (i = 0; i < nrow; i++)
                ra[i*nrow + i] -= rd[i];
            break;
        case RATLAS_DCOMPLEX:
            ca = (Ratlas_Complex *) amat->data;
            cd = (Ratlas_Complex *) dmat->data;
            for (i = 0; i < nrow; i++) {
                cursor = i*nrow + i; 
                ca[cursor].r -= cd[i].r;
                ca[cursor].i -= cd[i].i;
            }
            break;
        default:
            rb_raise(rb_eArgError, "Unknown type.");
    }
    return astor;
}




/**
 * pivots columns of a m x n matrix.  pivoting is given by array piv, which
 * must have n integer elements. */
VALUE ratlas_colpiv_bang(VALUE self, VALUE stor, VALUE piv)
{
    Ratlas_Matrix *mat;
/*     struct RArray *pivarr; */
    double *work, *r;
    Ratlas_Complex *c, *cwork;
    int i, colstride, colbytes, pivlen;

    Data_Get_Struct(stor, Ratlas_Matrix, mat);
/*     pivarr = RARRAY(piv); */
    Check_Type(piv, T_ARRAY);
    pivlen = RARRAY_LEN(piv);
    if (pivlen != mat->ncol)
        rb_raise(rb_eArgError, "Number of pivots must equal columns.");
    
    colstride = mat->nrow;
    colbytes = colstride * ratlas_sizeof[mat->type];
    switch (mat->type) {
        case RATLAS_DFLOAT:
            r = (double *) mat->data;
            work = ALLOCA_N(double, colstride);
            for (i = 0; i < pivlen; i++) {
                if ( i == NUM2INT(RARRAY_PTR(piv)[i]))
                    continue;
                memcpy(work, r+i*colstride, colbytes);
            }
            break;
        case RATLAS_DCOMPLEX:
            c = (Ratlas_Complex *) mat->data;
            cwork = ALLOCA_N(Ratlas_Complex, colstride);
            for (i = 0; i < pivlen; i++) {
                if ( i == NUM2INT(RARRAY_PTR(piv)[i]))
                    continue;
                memcpy(cwork, c+i*colstride, colbytes);
            }
            break;
    }
    return stor;
}



/**
 * allocates and returns a new matrix C = transpose(A) */
VALUE ratlas_transpose(VALUE self, VALUE stor)
{
    Ratlas_Matrix *orig, *trans;
    unsigned long int i, nelem, idx;
    int m, n;
    double *rt, *r;
    
    Data_Get_Struct(stor, Ratlas_Matrix, orig);
    m = orig->nrow;
    n = orig->ncol;

    if (m==1 || n==1) {
        trans = ratlas_matrix_alloc(orig->type, orig->matrixtype, n, m, 
                orig->data);
    } else {
        trans = ratlas_matrix_alloc(orig->type, orig->matrixtype, n, m, NULL);
        nelem = m*n;
        switch (orig->type) {
            case RATLAS_DFLOAT:
                r = (double *) orig->data;
                rt = (double *) trans->data;
                for (i = 0; i < nelem; i++) {
                    idx = (i%m)*n + i/m;
                    rt[idx] = r[i];
                }
                break;
            case RATLAS_DCOMPLEX:
                rb_raise(rb_eArgError, "No complex support yet.");
                break;
        }
    }
    return Data_Wrap_Struct(ratlas_storage_class, 0,
                                ratlas_matrix_free, trans);
}



/**
 * Expands a upper or lower triangular symmetric matrix to the full symmetric
 * matrix.  The matrix argument has memory for holding the full matrix, but
 * symmetric operations only reference the upper or lower triangular part.
 * Results need to be copied over if one is to do non-symmetric operations
 * later. */
void ratlas_matrix_sym2dense(int uplo, Ratlas_Matrix *mat)
{
    int i, j, ncol;
    unsigned long int iup, ilo;
    double *rcur;
    Ratlas_Complex *ccur;
    
    if (mat->ncol != mat->nrow)
        rb_raise(rb_eArgError, "Matrix must be square.");
    
    ncol = mat->ncol;
    switch (uplo) {
        case CblasLower:
            switch (mat->type) {
                case RATLAS_DFLOAT:
                    rcur = (double *) mat->data;
                    for (j = 1; j < ncol; j++) {
                        for (i = 0; i < ncol; i++) {
                            if (i == j) break;
                            iup = i+j*ncol;
                            ilo = i*ncol+j;
                            rcur[iup] = rcur[ilo];
                        }
                    }
                    break;
                case RATLAS_DCOMPLEX:
                    /* Note: make a Hermitian */
                    ccur = (Ratlas_Complex *) mat->data;
                    for (j = 1; j < ncol; j++) {
                        for (i = 0; i < ncol; i++) {
                            if (i == j) break;
                            iup = i+j*ncol;
                            ilo = i*ncol+j;
                            ccur[iup].r = ccur[ilo].r;
                            ccur[iup].i = - ccur[ilo].i;
                        }
                    }
                    break;
            }
            break;
        case CblasUpper:
            switch (mat->type) {
                case RATLAS_DFLOAT:
                    rcur = (double *) mat->data;
                    for (j = 1; j < ncol; j++) {
                        for (i = 0; i < ncol; i++) {
                            if (i == j) break;
                            iup = i+j*ncol;
                            ilo = i*ncol+j;
                            rcur[ilo] = rcur[iup];
                        }
                    }
                    break;
                case RATLAS_DCOMPLEX:
                    /* Note: make a Hermitian */
                    ccur = (Ratlas_Complex *) mat->data;
                    for (j = 1; j < ncol; j++) {
                        for (i = 0; i < ncol; i++) {
                            if (i == j) break;
                            iup = i+j*ncol;
                            ilo = i*ncol+j;
                            ccur[ilo].r = ccur[iup].r;
                            ccur[ilo].i = - ccur[iup].i;
                        }
                    }
                    break;
            }
            break;
    }
}




VALUE ratlas_storage_sym2dense_bang(VALUE self, VALUE uplo, VALUE stor)
{
    Ratlas_Matrix *mat;
    int uploflag;
    
    uploflag = FIX2INT(uplo);
    Data_Get_Struct(stor, Ratlas_Matrix, mat);
    ratlas_matrix_sym2dense(uploflag, mat);
    return stor;
}




/**
 * For lower triangular matrix, set upper triagonal to zero.
 * For upper triangular matrix, set lower triagonal to zero.
 * If diagonal is to be unit, set it to unit. */
void ratlas_matrix_uplo_nullify(int unitdiag, int uplo, Ratlas_Matrix *mat)
{
    int i, j, ndiag;
    unsigned long int idx;
    double *r;
    Ratlas_Complex *c;
    
    switch (uplo) {
        case CblasLower:
            switch (mat->type) {
                case RATLAS_DFLOAT:
                    r = (double *) mat->data;
                    for (j = 1; j < mat->ncol; j++) {
                        for (i = 0; i < mat->nrow; i++) {
                            if (i == j) break;
                            idx = i+j*mat->nrow;
                            r[idx] = 0;
                        }
                    }
                    break;
                case RATLAS_DCOMPLEX:
                    c = (Ratlas_Complex *) mat->data;
                    for (j = 1; j < mat->ncol; j++) {
                        for (i = 0; i < mat->nrow; i++) {
                            if (i == j) break;
                            idx = i+j*mat->nrow;
                            c[idx].r = 0;
                            c[idx].i = 0;
                        }
                    }
                    break;
            }
            break;
        case CblasUpper:
            switch (mat->type) {
                case RATLAS_DFLOAT:
                    r = (double *) mat->data;
                    for (i = 1; i < mat->nrow; i++) {
                        for (j = 0; j < mat->ncol; j++) {
                            if (i == j) break;
                            idx = i+j*mat->nrow;
                            r[idx] = 0;
                        }
                    }
                    break;
                case RATLAS_DCOMPLEX:
                    c = (Ratlas_Complex *) mat->data;
                    for (i = 1; i < mat->nrow; i++) {
                        for (j = 0; j < mat->ncol; j++) {
                            if (i == j) break;
                            idx = i+j*mat->nrow;
                            c[idx].r = 0;
                            c[idx].i = 0;
                        }
                    }
                    break;
            }
            break;

    }

    if (unitdiag == CblasUnit)
    {
        ndiag = ratlas_imin(mat->nrow, mat->ncol);
        switch (mat->type) {
            case RATLAS_DFLOAT:
                r = (double *) mat->data;
                for (i = 0; i < ndiag; i++)
                    r[mat->nrow*i + i] = 1.0;
                break;
            case RATLAS_DCOMPLEX:
                c = (Ratlas_Complex *) mat->data;
                for (i = 0; i < ndiag; i++) {
                    idx = mat->nrow*i + i;
                    c[idx].r = 1.0;
                    c[idx].i = 0.0;
                }
                break;

        }
    }
}




/**
 * Expands a diagonal matrix, stored as a 1-column matrix, to a full diagonal
 * dense matrix */
VALUE ratlas_storage_diag2dense(VALUE self, VALUE dstor)
{
    int i, nrow;
    unsigned long int nelem;
    double *zeros, *rm, *rd;
    Ratlas_Complex *cm, *cd;
    Ratlas_Matrix *mat, *diag;
    VALUE retval=Qnil;
    
    Data_Get_Struct(dstor, Ratlas_Matrix, diag);
    nrow = diag->nrow;
    nelem = nrow*nrow;

    switch (diag->type) {
        case RATLAS_DFLOAT:
            zeros = xcalloc(nelem, ratlas_sizeof[RATLAS_DFLOAT]);
            mat = ratlas_matrix_alloc(RATLAS_DFLOAT, RATLAS_DENSE,
                    nrow, nrow, zeros);
            free(zeros);
            rm = (double *) mat->data;
            rd = (double *) diag->data;
            for (i = 0; i < nelem; i++)
                rm[i+i*nrow] = rd[i];    
            retval = Data_Wrap_Struct(ratlas_storage_class, 0, 
                    ratlas_matrix_free, mat);
            break;
        case RATLAS_DCOMPLEX:
            zeros = xcalloc(nelem, ratlas_sizeof[RATLAS_DCOMPLEX]);
            mat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE,
                    nrow, nrow, zeros);
            free(zeros);
            cm = (Ratlas_Complex *) mat->data;
            cd = (Ratlas_Complex *) diag->data;
            for (i = 0; i < nelem; i++) {
                cm[i+i*nrow].r = cd[i].r;
                cm[i+i*nrow].i = cd[i].i;
            }
            retval = Data_Wrap_Struct(ratlas_storage_class, 0,
                    ratlas_matrix_free, mat);
            break;
    }
    return retval;
}





/**
 * Casts a real matrix to complex matrix. */
VALUE ratlas_storage_real2complex(VALUE self, VALUE restor)
{
    Ratlas_Matrix *remat;
    Ratlas_Matrix *cplxmat;
    double *r;
    Ratlas_Complex *c;
    unsigned long int i, nelem;
    
    Data_Get_Struct(restor, Ratlas_Matrix, remat);
    cplxmat = ratlas_matrix_alloc(RATLAS_DCOMPLEX, RATLAS_DENSE,
            remat->nrow, remat->ncol, NULL);
    
    nelem = remat->nrow * remat->ncol;
    r = (double *) remat->data;
    c = (Ratlas_Complex *) cplxmat->data;
    for (i = 0; i < nelem; i++) {
        c[i].r = r[i];
        c[i].i = 0.0;
    }

    return Data_Wrap_Struct(ratlas_storage_class, 0, ratlas_matrix_free,
            cplxmat);
}
