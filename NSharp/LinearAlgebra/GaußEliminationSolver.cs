using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.LinearAlgebra
{
    public class SolverException : Exception
    {
        public SolverException()
        {
        }

        public SolverException(string message)
            : base(message)
        {
        }

        public SolverException(string message, Exception inner)
            : base(message, inner)
        {
        }
    }

    /// <summary>
    /// Equation Solver using Gauß Jordan Elimination
    /// Author: Muhammed Demircan
    /// Date: 16.06.15
    /// </summary>
    public class GaußEliminationSolver
    {
        /// <summary>
        /// Solves the equation Ax = b by using Gauß Elimination Algorithm for a square matrix A
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vector Solve(Matrix A, Vector b)
        {
            int N = A.NoRows;

            if (A.NoRows != A.NoColumns)
                throw new SolverException("Matrix A is not square.");


            //Transform to lower triangular matrix
            for (int i = 0; i < N; i++)
            {
                int pivotIdx = A.GetPivotRow(i, i);
                if (A[pivotIdx,i] == 0.00)
                    throw new SolverException("Matrix A is not singular");

                A.SwapRow(i, pivotIdx);
                b.SwapRow(i, pivotIdx);

                for (int row = i + 1; row < N; row++)
                {
                    for (int col = i + 1; col < N; col++)
                    {
                        A[row, col] = A[row, col] - A[i, col] * ((A[row, i] / A[i, i]));
                        
                    }

                    //Update vecctor b
                    b[row, 0] = b[row, 0] - b[i, 0] * ((A[row, i] / A[i, i]));
                    
                    //Set lower Triangular Matrix with zeros
                    A[row, i] = 0.0;
                }
            }

            //Transform lower triangular matrix to diagnol matrix
            for (int i = N-1; i >= 0; i--)
            {
                //Normalize
                b[i,0] = b[i,0]/A[i,i];
                A[i, i] = 1.0;
                for (int row = i - 1; row >= 0; row--)
                {       
                    b[row, 0] = b[row, 0] - b[i, 0] * (A[row, i]);
                    A[row, i] = 0.0;
                }
            }

            return b;
        }
    }
}
