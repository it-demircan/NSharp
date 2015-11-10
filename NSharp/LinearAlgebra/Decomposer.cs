using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.LinearAlgebra
{
    public enum DecomposerType
    {
        LU
    }

    /// <summary>
    /// This Class implements algorithms to decompose matrices
    /// </summary>
    public class Decomposer
    {

        /// <summary>
        /// Public method to access the algorithms, which decompose the passed matrix by using the selected method.
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public static Matrix[] Decompose(Matrix mat, DecomposerType type){
            switch (type)
            {
                case DecomposerType.LU:
                    return DecomposeUsingLU(mat);
                default:
                    return DecomposeUsingLU(mat);
            }
        }

        /// <summary>
        /// Decompse a matrix which has full rank into a lower and a upper triangular matrix L*U = A.
        /// According to Crout's Method.
        /// </summary>
        /// <param name="mat">Matrix A</param>
        /// <returns>Matrix Array with lower triangular matrix [0] and upper triagular matrix U [1]</returns>
        private static Matrix[] DecomposeUsingLU(Matrix mat)
        {
            int N = mat.NoColumns;
            if(N != mat.NoRows)
                throw new Exception("Not a square matrix.");
            Matrix[] decomposition = new Matrix[2];
            Matrix lowerMatrix = new Matrix(N,N);
            Matrix upperMatrix = new Matrix(N,N);
            double temp;


            //Crout Algorithm
            for (int i = 0; i < N; i++)
            {
                lowerMatrix[i, i] = 1.0;
            }
         
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i <= j; i++)
                {              
                    temp = 0.0;
                    for (int k = 0; k < i; k++)
                    {
                        temp += lowerMatrix[i, k] * upperMatrix[k, j];
                    }
                    upperMatrix[i, j] = mat[i, j] - temp;
                }

                for (int i = j + 1; i < N; i++)
                {
                    temp = 0.0;
                    for(int k = 0 ; k < j; k++){
                        temp -= lowerMatrix[i,k]*upperMatrix[k,j];
                    }
                    temp += mat[i,j];
                    lowerMatrix[i, j] = temp / upperMatrix[j, j];
                }
            }
            decomposition[0] = lowerMatrix;
            decomposition[1] = upperMatrix;

            return decomposition;
        }
    }
}
