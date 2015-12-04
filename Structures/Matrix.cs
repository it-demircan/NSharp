using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Structures
{
    public enum MatrixType
    {
        ONES, EYE
    }

    /// <summary>
    /// Exception Class, which is used for error and exceptions thrown in class "Matrix"
    /// </summary>
    public class MatrixException : Exception
    {
        public MatrixException()
        {
        }

        public MatrixException(string message)
            : base(message)
        {
        }

        public MatrixException(string message, Exception inner)
            : base(message, inner)
        {
        }

    }

    /* Structure for a n x m Matrix.
     * Autor : Muhammed Demircan
     * */
    public class Matrix : ICloneable
    {
        private double[,] values;
        private int noRows;
        private int noColumns;

        /// <summary>
        /// Constructor 
        /// </summary>
        /// <param name="_n">Number of rows</param>
        /// <param name="_m">Number of columns</param>
        public Matrix(int _n, int _m)
        {
            noRows = _n;
            noColumns = _m;
            values = new double[noRows, noColumns];
        }

        /// <summary>
        /// Constructs a matrix with a special structure, like an "Eye" Matrix (see Matlab)
        /// </summary>
        /// <param name="_n">Number of rows</param>
        /// <param name="_m">Number of columns</param>
        /// <param name="structure">Structure Type </param>
        public Matrix(int _n, int _m, MatrixType type)
        {
            noRows = _n;
            noColumns = _m;
            values = new double[noRows, noColumns];

            switch (type)
            {
                case MatrixType.EYE:
                    for (int i = 0; i < noRows; i++)
                        for (int j = 0; j < noColumns; j++)
                            if (i == j)
                                this[i, j] = 1;
                    break;
                case MatrixType.ONES:
                    for (int i = 0; i < noRows; i++)
                        for (int j = 0; j < noColumns; j++)
                            this[i, j] = 1;
                    break;
            }
        }

        /// <summary>
        /// Creates a diagonal matrix using values from input vector for diagonal entries.
        /// </summary>
        /// <param name="_diag">Vector, which should be mapped as a diagonal matrix.</param>
        public Matrix(Vector _diag)
        {
            noRows = _diag.Length;
            noColumns = noRows;
            values = new double[noRows, noColumns];

            for (int i = 0; i < noRows; i++)
                values[i, i] = _diag[i];
        }

        /// <summary>
        /// Multiplicates two matrices A,B
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C = A*B</returns>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            int n_mat1 = A.NoRows;
            int m_mat1 = A.NoColumns;

            int n_mat2 = B.NoRows;
            int m_mat2 = B.noColumns;

            //Dimension missmatch
            if (m_mat1 != n_mat2)
            {
                throw new Exception("Dimension missmatch");
            }

            Matrix resultMatrix = new Matrix(n_mat1, m_mat2);

            for (int i = 0; i < n_mat1; i++)
            {
                for (int j = 0; j < m_mat2; j++)
                {
                    for (int k = 0; k < n_mat2; k++)
                    {
                        resultMatrix[i, j] += A[i, k] * B[k, j];
                    }
                }
            }
            return resultMatrix;
        }

        /// <summary>
        /// Sum of two matrices A,B
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C= A+B</returns>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            int n_mat1 = A.NoRows;
            int m_mat1 = A.NoColumns;

            int n_mat2 = B.NoRows;
            int m_mat2 = B.noColumns;

            //Dimension missmatch
            if (m_mat1 != m_mat2 || n_mat1 != n_mat2)
            {
                throw new Exception("Dimension missmatch");
            }

            Matrix resultMatrix = new Matrix(n_mat1, m_mat2);

            for (int i = 0; i < n_mat1; i++)
            {
                for (int j = 0; j < m_mat2; j++)
                {
                    resultMatrix[i, j] += A[i, j] + B[i, j];
                }
            }
            return resultMatrix;
        }

        /// <summary>
        /// Subtract matrices A,B
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C= A-B</returns>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            int n_mat1 = A.NoRows;
            int m_mat1 = A.NoColumns;

            int n_mat2 = B.NoRows;
            int m_mat2 = B.noColumns;

            //Dimension missmatch
            if (m_mat1 != m_mat2 || n_mat1 != n_mat2)
            {
                throw new Exception("Dimension missmatch");
            }

            Matrix resultMatrix = new Matrix(n_mat1, m_mat2);

            for (int i = 0; i < n_mat1; i++)
            {
                for (int j = 0; j < m_mat2; j++)
                {
                    resultMatrix[i, j] += A[i, j] - B[i, j];
                }
            }
            return resultMatrix;
        }

        /// <summary>
        /// Matrix Vector Multiplication
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="B">Vektor B</param>
        /// <returns>Matrix C = A*B</returns>
        public static Vector operator *(Matrix A, Vector B)
        {
            if (A.NoColumns != B.NoRows)
                throw new Exception("Dimension missmatch.");

            Vector result = new Vector(A.noRows);

            for (int i = 0; i < A.noRows; i++)
            {
                for (int j = 0; j < A.noColumns; j++)
                {
                    result[i] += A[i, j] * B[j];
                }
            }
            return result;
        }

        /// <summary>
        /// Vector Matrix Multiplication
        /// </summary>
        /// <param name="A">Vektor A[transponiert]</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Vektor C = A*B</returns>
        public static Vector operator *(Vector A, Matrix B)
        {
            if (A.NoRows != B.noRows)
                throw new Exception("Dimension missmatch.");

            Vector result = new Vector(B.NoColumns);


            for (int i = 0; i < B.NoColumns; i++)
                for (int j = 0; j < B.noRows; j++)
                {
                    result[i] += A[j] * B[j, i];
                }
            return result;
        }


        /// <summary>
        /// Skalar x Matrix Multiplication
        /// </summary>
        /// <param name="A">Skalar A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C = a*B</returns>
        public static Matrix operator *(double a, Matrix B)
        {

            Matrix result = new Matrix(B.NoRows, B.NoColumns);
            for (int i = 0; i < B.NoColumns; i++)
                for (int j = 0; j < B.noRows; j++)
                {
                    result[i, j] = a * B[i, j];
                }
            return result;
        }

        /// <summary>
        /// Skalar x Matrix Multiplication
        /// </summary>
        /// <param name="A">Skalar A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C = a*B</returns>
        public static Matrix operator *(Matrix B, double a)
        {

            Matrix result = new Matrix(B.NoRows, B.NoColumns);
            for (int i = 0; i < B.NoColumns; i++)
                for (int j = 0; j < B.noRows; j++)
                {
                    result[i, j] = a * B[i, j];
                }
            return result;
        }

        /// <summary>
        /// Transpose matrix A
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <returns>A ^{T}</returns>
        public static Matrix operator !(Matrix A)
        {
            Matrix result = new Matrix(A.noColumns, A.noRows);
            for (int j = 0; j < A.noRows; j++)
            {
                for (int i = 0; i < A.noColumns; i++)
                {
                    result[i, j] = A[j, i];
                }
            }
            return result;
        }

        /// <summary>
        /// Swaps two rows at given position "sourceRow" and "destinationRow"
        /// </summary>
        /// <param name="sourceRow">Row Index between 0 and NoRows-1</param>
        /// <param name="destinationRow">Row Index between 0 and NoRows-1</param>
        public void SwapRow(int sourceRow, int destinationRow)
        {
            if (sourceRow == destinationRow)
                return;
            if (sourceRow > this.NoRows || destinationRow > this.NoRows || sourceRow < 0 || destinationRow < 0)
            {
                throw new MatrixException("[IndexOutOfBound] - Can't swap rows for given indices.");
            }

            double val;
            for (int i = 0; i < this.NoColumns; i++)
            {
                val = this[sourceRow, i];
                this[sourceRow, i] = this[destinationRow, i];
                this[destinationRow, i] = val;
            }
        }

        /// <summary>
        /// Swaps two columns at given position "sourceColumn" and "destinationColumn"
        /// </summary>
        /// <param name="sourceRow">Row Index between 0 and NoColumn-1</param>
        /// <param name="destinationRow">Row Index between 0 and NoColumn-1</param>
        public void SwapColumn(int sourceColumn, int destinationColumn)
        {
            if (sourceColumn == destinationColumn)
                return;
            if (sourceColumn > this.NoColumns || destinationColumn > this.NoColumns || sourceColumn < 0 || destinationColumn < 0)
            {
                throw new MatrixException("[IndexOutOfBound] - Can't swap columns for given indices.");
            }

            double val;
            for (int i = 0; i < this.NoRows; i++)
            {
                val = this[i,sourceColumn];
                this[sourceColumn, i] = this[i,destinationColumn];
                this[i,destinationColumn] = val;
            }
        }


        /// <summary>
        /// Finds a pivot element for a specific column beginnend at a given row index until last row.
        /// </summary>
        /// <param name="downWardsRow">Start row</param>
        /// <param name="columnIdx">Decisive column index </param>
        /// <returns>Row index of pivot element</returns>
        public int GetPivotRow(int downWardsRow, int columnIdx)
        {
            int pivotRow = downWardsRow;
            for (int i = downWardsRow; i < this.NoRows; i++)
            {
                if (Math.Abs(this[pivotRow, columnIdx]) < Math.Abs(this[i, columnIdx]))
                {
                    pivotRow = i;
                }
            }

            return pivotRow;
        }

        /// <summary>
        /// Getter for a specific column
        /// </summary>
        /// <param name="columnIdx">Column Index</param>
        /// <returns>Column as a Vector</returns>
        public Vector GetColumn(int columnIdx)
        {
            double[] colValues = new double[this.NoRows];

            for (int i = 0; i < this.NoRows; i++)
                colValues[i] = this[i, columnIdx];

            return new Vector(colValues);
        }

        /// <summary>
        /// Accessor for an element in matrix A at position [x,y]
        /// </summary>
        /// <param name="x">Row index, starts at 0 </param>
        /// <param name="y">Column index, start at 0</param>
        /// <returns>Value at [x,y] </returns>
        public double this[int x, int y]
        {
            get { return values[x, y]; }
            set { values[x, y] = value; }
        }


        /// <summary>
        /// Number of rows
        /// </summary>
        public int NoRows
        {
            get { return noRows; }
        }
        /// <summary>
        /// Number of columns
        /// </summary>
        public int NoColumns
        {
            get { return noColumns; }
        }

        public void InjectMatrixAtPosition(Matrix matrixToBeInjected, int rowPosition, int columnPosition)
        {
            for (int i = 0; i < matrixToBeInjected.NoRows; i++)
            {
                for (int k = 0; k < matrixToBeInjected.NoColumns; k++)
                {
                    this[rowPosition + i, columnPosition + k] = matrixToBeInjected[i, k];
                }
            }
        }

        /// <summary>
        /// Create an Outputstring containing all values in recent matrix object.
        /// </summary>
        /// <param name="decimalRound">Number of digits (rounding)</param>
        /// <returns>Matrix values in an formatted string</returns>
        public string toString(int decimalRound)
        {
            string outputString = "";

            for (int i = 0; i < noRows; i++)
            {
                for (int j = 0; j < noColumns; j++)
                {
                    outputString += " [" + Math.Round(this.values[i, j], decimalRound) + "] ";
                }
                outputString += "\r\n";
            }
            return outputString;
        }

        public object Clone()
        {
            Matrix clone = new Matrix(this.NoRows, this.NoColumns);

            for (int i = 0; i < this.NoRows; i++)
            {
                for (int j = 0; j < this.NoColumns; j++)
                {
                    clone[i, j] = this[i, j];
                }
            }
            return clone;
        }
    }
}
