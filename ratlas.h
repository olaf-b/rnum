/*
 *  ratlas.h
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#ifndef RATLAS_H
#define RATLAS_H

#include <stddef.h>
#include <stdio.h>
#include <ruby.h>
#include "ratlas_config.h"
#include <math.h>
/* Ratlas error codes */
#define RATLAS_ERROR     -1
#define PROPAGATESYM 1 /* undef when full support for symmetric matrices is
              implemented */

/* Matrix element types */
#define    RATLAS_DFLOAT    0
#define    RATLAS_DCOMPLEX  1
#define    RATLAS_NTYPES    2

/* Matrix types */
#define    RATLAS_DENSE        0xa0
#define    RATLAS_BAND         0xa1
#define    RATLAS_SYM          0xa2
#define RATLAS_SYMBAND      0xa3
#define RATLAS_SYMPACK      0xa4
#define RATLAS_TRI          0xa5
#define RATLAS_TRIBAND      0xa6
#define RATLAS_TRIPACK      0xa7

typedef struct Ratlas_Matrix Ratlas_Matrix;
typedef struct Ratlas_Complex Ratlas_Complex;

extern ID id_Complex;
extern ID id_Range;
extern ID id_atre;
extern ID id_atim;
extern ID id_imag;
extern ID id_call;
extern ID id_first;
extern ID id_last;
extern ID id_exclude_end;

extern VALUE complex_class;
extern VALUE range_class;
extern VALUE ratlas_storage_class;

struct Ratlas_Matrix {
    unsigned int nrow;         // Number of rows.
    unsigned int ncol;         // Number of columns.
    void *data;                // Pointer to data.
    int type;                  // Ratlas_Types
    int matrixtype;            // Ratlas_Matrix_Types
    unsigned long int memsize; // Number of bytes in memory for data;
};

struct Ratlas_Complex {
    double r, i;
};

#include "ratlas_func.h"
#include "ratlas_cblas.h"
#include "ratlas_lapack.h"

#endif /* ifndef RATLAS_H */
