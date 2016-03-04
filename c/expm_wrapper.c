#include <math.h>
#include <libhades.h>
#include <libhades/expm.h>

void expm_(complex_t *A, int *dim_p);

void expm_(complex_t *A, int *dim_p)
{
    const int dim = *dim_p;
    matrix_complex_t M;

    M.rows    = dim;
    M.columns = dim;
    M.min     = dim;
    M.size    = dim*dim;
    M.type    = 0;
    M.view    = 0;
    M.M       = A;

    matrix_complex_expm(&M);
}
