#include <stdio.h>
#include <gsl/gsl_blas.h>

int
main (void)
{
  double a[] = { 0.11, 0.12, 0.13,
                 0.21, 0.22, 0.23 };

  double b[] = { 1011, 1012,
                 1021, 1022,
                 1031, 1032 };

  double c[] = { 0.00, 0.00,
                 0.00, 0.00 };

  gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
  gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
  gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);

  /* Compute C = A B */

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

  printf ("[ %g, %g\n", c[0], c[1]);
  printf ("  %g, %g ]\n", c[2], c[3]);

  return 0;  
}
//Function: int gsl_blas_dger (double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A)
//void gsl_matrix_set_all (gsl_matrix * m, double x)//set all element to zero
//void gsl_matrix_set_identity (gsl_matrix * m)//make identity matrix
//gsl_vector_view gsl_matrix_diagonal (gsl_matrix * m) //return a vector view of the diagonal of the matrix m
Function: int gsl_blas_sgemv (CBLAS_TRANSPOSE_t TransA, float alpha, const gsl_matrix_float * A, const gsl_vector_float * x, float beta, gsl_vector_float * y)
Function: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
Function: int gsl_blas_cgemv (CBLAS_TRANSPOSE_t TransA, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_vector_complex_float * x, const gsl_complex_float beta, gsl_vector_complex_float * y)
Function: int gsl_blas_zgemv (CBLAS_TRANSPOSE_t TransA, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_vector_complex * x, const gsl_complex beta, gsl_vector_complex * y)

Function: int gsl_blas_sdsdot (float alpha, const gsl_vector_float * x, const gsl_vector_float * y, float * result)
This function computes the sum \alpha + x^T y for the vectors x and y, returning the result in result.

Function: int gsl_blas_sdot (const gsl_vector_float * x, const gsl_vector_float * y, float * result)
Function: int gsl_blas_dsdot (const gsl_vector_float * x, const gsl_vector_float * y, double * result)
Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
These functions compute the scalar product x^T y for the vectors x and y, returning the result in result.

gsl_matrix * my_diag_alloc(gsl_vector * X)
{
    gsl_matrix * mat = gsl_matrix_alloc(X->size, X->size);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_matrix_set_all(mat, 0.0); //or whatever number you like
    gsl_vector_memcpy(&diag.vector, X);
    return mat;
}
