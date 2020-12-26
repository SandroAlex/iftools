/* 
===============================================================================
Basic functions written in C programming language.
===============================================================================
*/

// Load libraries:
#include <stdio.h>
#include <math.h>

// 1d array size macros.
#define array1d_size1(a1d) (sizeof(a1d) / sizeof(a1d[0])) // Array size.

// 2d array size macros.
#define array2d_size1(a2d) (sizeof(a2d) / sizeof(a2D[0]))       // Number of rows.
#define array2d_size2(a2d) (sizeof(a2d[0]) / sizeof(a2D[0][0])) // Number of columns. 

// 3d array size macros.
#define array3d_size1(a3d) (sizeof(a3d) / sizeof(a3d[0]))             // Number of matrices. 
#define array3d_size2(a3d) (sizeof(a3d[0]) / sizeof(a3d[0][0]))       // Number of rows.
#define array3d_size3(a3d) (sizeof(a3d[0][0]) / sizeof(a3d[0][0][0])) // Number of columns. 

// Functions prototypes.
////////////////////////
float transfer_entropy_summation(
    int nbins_xnn, int nbins_xn, int nbins_yn,
    float p_xnn_xn_yn[nbins_xnn][nbins_xn][nbins_yn], 
    float p_xnn_xn[nbins_xnn][nbins_xn], 
    float p_xn_yn[nbins_xn][nbins_yn], 
    float p_xn[nbins_xn]);

/////////////////////////
int isclose(float a, float b);

// Main program.
////////////////
int main(){

    return 0;
}

// Functions definitions.
/////////////////////////
float transfer_entropy_summation(
    int nbins_xnn, int nbins_xn, int nbins_yn,
    float p_xnn_xn_yn[nbins_xnn][nbins_xn][nbins_yn], 
    float p_xnn_xn[nbins_xnn][nbins_xn], 
    float p_xn_yn[nbins_xn][nbins_yn], 
    float p_xn[nbins_xn]){
    
    // Variables definitions.    
    float transfer_entropy=0.0;
    float zero=0.0, arg;

    // Calculate transfer entropy.
    for (int i=0; i < nbins_xnn; i++){       // xnn.
        for(int j=0; j < nbins_xn; j++){     // xn.
            for(int k=0; k < nbins_yn; k++){ // yn.

                // We are going to add only non zero terms.
                // Assumptions:  0*log2(a/0)=0.
                //               0*log2(0) = 0.
                // Here we need to further investigate these assumptions.
                if (!isclose(p_xnn_xn_yn[i][j][k], zero) && 
                    !isclose(p_xn[j], zero) &&
                    !isclose(p_xn_yn[j][k], zero) &&
                    !isclose(p_xnn_xn[i][j], zero)){
                    
                        // Argument for log2 function.
                        arg = p_xnn_xn_yn[i][j][k] * p_xn[j] /
                              p_xn_yn[j][k] / p_xnn_xn[i][j];
                       
                        // Accumulate for transfer entropy.
                        transfer_entropy += p_xnn_xn_yn[i][j][k] * log2f(arg);
                }
            }
        }
    }
    return transfer_entropy;
}

/////////////////////////
int isclose(float a, float b){

    // Definitions of variables.
    int close=1, not_close=0;
    float atol=1.00e-8, rtol=1.00e-5; // Absolute and relative tolerances. 
    float left, right;  // Left and right sides of the following equation:
                        // absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    left = fabs(a - b);
    right = atol + rtol * fabs(b);

    // They are close. Return true.
    if (left <= right){
        return close;
    }

    // Not close. Return false.
    else{
        return not_close;
    }
}