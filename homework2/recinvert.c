#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <time.h>

#define N 10000
#define P 100

float* recinvert(float* X, int rows, int columns, int ldx, float* res) {
  if(columns == 2) {
    // Create a new array
    float array[4];

    // a = n
    array[0] = (float)rows;

    // b = sum(X[,2])
    // d = sum(X[,2]^2)
    array[1] = 0.0;
    array[3] = 0.0;
    for(int i = 0; i < rows; i++){
      array[1] += X[ldx*i+1];
      array[3] += X[ldx*i+1]*X[ldx*i+1];
    }
    array[2] = array[1];

    // Calculate determinant
    float det = array[0]*array[3] - array[1]*array[2];

    // Inverse the matrix
    //float* invarray = malloc(4*sizeof(float));
    float* invarray = res;
    invarray[0*ldx] = array[3]/det;
    invarray[0*ldx + 1] = -array[1]/det;
    invarray[1*ldx] = invarray[1];
    invarray[1*ldx+1] = array[0]/det;

    return invarray;
  }

  int prevrow = rows-1;
  int prevcol = columns-1;
  
  float B[rows-1];
  cblas_sgemv(CblasRowMajor, CblasTrans, rows, prevcol, 1.0, X, ldx, &X[prevcol], ldx, 0, B, 1);

  /*
  printf("** B **\n");
  for(int i = 0; i < prevrow; i++){
    printf("[%f] ", B[i]);
  }
  printf("\n\n");
  */

  float d = cblas_sdsdot(rows, 0, &X[prevcol], ldx, &X[prevcol], ldx);

  //printf("[%f]", d);

  float* invA = recinvert(X, rows, prevcol, ldx, res);

  /*
  printf("**INV A**\n");
  for(int i = 0; i < columns - 1; i++){
    for(int j = 0; j < columns -1; j++){
      printf("%f ", invA[i*ldx + j]);
    }
    printf("\n");
  }
  printf("\n");
  */

  float Bt_invA[prevcol];
  cblas_sgemv(CblasRowMajor, CblasTrans, prevcol, prevcol, 1.0, invA, ldx, B, 1, 0, Bt_invA, 1);

  /*
  printf("** BT INV A **\n");
  for(int i = 0; i < columns -1; i++){
    printf("%f ", Bt_invA[i]);
  }
  printf("\n\n");
  */

  float* invA_B = Bt_invA;
  //float invA_B[columns-1];// = Bt_invA;
  //cblas_sgemv(CblasRowMajor, CblasNoTrans, columns-1, columns-1, 1.0, invA, columns-1, B, 1, 0, invA_B, 1);

  /*
  printf("** INV A B **\n");
  for(int i = 0; i < columns -1; i++){
    printf("%f ", invA_B[i]);
  }
  printf("\n\n");
  */

  float Bt_invA_B = cblas_sdsdot(prevcol, 0, Bt_invA, 1, B, 1);

  float k = d - Bt_invA_B;

  //printf("[%f]\n", k);

  float invA_B_Bt_invA[prevcol * prevcol];
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, prevcol, prevcol, 1, 1.0/k, invA_B, 1, Bt_invA, columns-1, 0.0, invA_B_Bt_invA, prevcol);

  /*
  for(int i = 0; i < columns - 1; i++){
    for(int j = 0; j < columns -1; j++){
      printf("%f ", invA_B_Bt_invA[i*(columns-1) + j]);
    }
    printf("\n");
  }
  printf("\n");
  */

  cblas_sscal(prevcol, -1.0/k, invA_B, 1);

  float* result = res;
  for(int i = 0; i < prevcol; i++){
    for(int j = 0; j < prevcol; j++){
      result[i*ldx + j] += invA_B_Bt_invA[i*(prevcol) + j];
    }
  }

  for(int i = 0; i < prevcol; i++){
    result[i*ldx + (prevcol)] = invA_B[i];
    result[(prevcol)*ldx + i] = Bt_invA[i];
  }

  result[(prevcol)*ldx + prevcol] = 1.0/k;

  return result;
}

int
main (void)
{
  srand(1);

  int ldx = P;

  float X[N*P];// = { 1, 2, 5, 9, 1, 3, 8, 12, 1, 4, 7, 4, 1, 3, 4, 2 };
  
  for(int i = 0; i < N*P; i++) {
    X[i] = ((float)rand()/(float)(RAND_MAX));
  }
  /*
  for(int i = 0; i < N; i++){
    for(int j = 0; j < P; j++){
      printf("%f ", X[j + i*ldx]);
    }
    printf("\n");
  }
  printf("\n");
  */

  float res[P*P];


  float XtX[P*P];
  clock_t begin = clock(); 
  //for(int i = 0; i < 50; i++){
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, P, P, N, 1.0, X, ldx, X, ldx, 0, XtX, P);
    lapack_int ipiv[N];
    lapack_int info;
    info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, P, P, XtX, P, ipiv);
    info = LAPACKE_sgetri(LAPACK_ROW_MAJOR, P, XtX, P, ipiv);
  //}
  clock_t end = clock();
  double time_spent_inv = (double)(end - begin) / CLOCKS_PER_SEC;

  begin = clock();
  //for(int i = 0; i < 50; i++){
    float* invX = recinvert(X, N, P, P, res);
  //}
  end = clock();
  double time_spent_rec = (double)(end - begin) / CLOCKS_PER_SEC;


/*
  for(int i = 0; i < P; i++){
    for(int j = 0; j < P; j++){
      printf("%f ", invX[j + i*ldx]);
    }
    printf("\n");
  }
*/
  printf("Total time recursive: %f\n", time_spent_rec);
  printf("Total time LAPACK LU inverse: %f\n", time_spent_inv);
  return 0;  
}

