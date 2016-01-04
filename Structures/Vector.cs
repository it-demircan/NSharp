using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Structures
{

    public enum VectorType
    {
        ONES
    }

    /// <summary>
    /// Diese Klasse spiegelt Vektoren wieder und ist eine Unterklasse von Matrix. Die üblichen Operatoren sind implementiert.
    /// 
    /// Autor : Muhammed Demircan
    /// Datum : 25.07.2014
    /// </summary>
    public class Vector : Matrix
    {
        //Standardkonstruktor
        public Vector(int _n): base(_n, 1)
        {
            //Mit Nullen gefüllter Vektor
        }

        //Konstruktor, der unmittelbar ein Double Array entgegen nimmt und daraus eine Vektor Objekt konstruiert.
        public Vector(double[] _values):base(_values.Length, 1)
        {
            for (int i = 0; i < base.NoRows; i++)
            {
                this[i] = _values[i];
            }
        }

        //Spezielle Vektoren erzeugen.
        public Vector(int _n, VectorType type): base(_n, 1)
        {
            switch (type)
            {
                case VectorType.ONES:
                    for (int i = 0; i < base.NoRows; i++)
                    {
                        this[i] = 1.0;
                    }
                    break;
            }
        }

        // Einheitsvektor mit einer Eins an der Position gegeben durch "index".
        public Vector(int _n, int index)
            : base(_n, 1)
        {
            this[index] = 1.0;
        }

        /// <summary>
        /// Addiert zwei Vektoren A,B
        /// </summary>
        /// <param name="_mat1">Vektor A</param>
        /// <param name="_mat2">Vektor B</param>
        /// <returns>Summe C = A*B</returns>
        public static Vector operator +(Vector A, Vector B)
        {
            
            //Dimension missmatch
            if (A.NoRows != B.NoRows)
            {
                throw new Exception("Dimension missmatch");
            }

            Vector result = new Vector(A.NoRows);
            for (int i = 0; i < A.NoRows; i++)
            {
                result[i] = A[i, 0] + B[i, 0];
            }
            return result;
        }

        /// <summary>
        /// Subtrahiert zwei Vektoren A,B
        /// </summary>
        /// <param name="_mat1">Vektor A</param>
        /// <param name="_mat2">Vektor B</param>
        /// <returns>Summe C = A*B</returns>
        public static Vector operator -(Vector A, Vector B)
        {

            //Dimension missmatch
            if (A.NoRows != B.NoRows)
            {
                throw new Exception("Dimension missmatch");
            }

            Vector result = new Vector(A.NoRows);
            for (int i = 0; i < A.NoRows; i++)
            {
                result[i] = A[i, 0] - B[i, 0];
            }
            return result;
        }

        /// <summary>
        /// Multipliziert zwei Vektoren A,B
        /// </summary>
        /// <param name="_mat1">Vektor A</param>
        /// <param name="_mat2">Vektor B</param>
        /// <returns>Skalarprodukt C = A*B</returns>
        public static double operator *(Vector A, Vector B)
        {
            double innerProd = 0.0;
            //Dimension missmatch
            if (A.NoRows != B.NoRows)
            {
                throw new Exception("Dimension missmatch");
            }
            for (int i = 0; i < A.NoRows;i++)
            {
                innerProd += A[i, 0] * B[i, 0];
            }
            return innerProd;
        }


        /// <summary>
        /// Verkettet zwei Vektoren A,B
        /// </summary>
        /// <param name="_mat1">Vektor A</param>
        /// <param name="_mat2">Vektor B</param>
        /// <returns>Verkettzng A->B</returns>
        public static Vector operator &(Vector A, Vector B)
        {
            Vector result = new Vector(A.Length + B.Length);

            for (int i = 0; i < A.Length; i++)
            {
                result[i] = A[i];
            }

            for (int i = A.Length; i < A.Length + B.Length; i++)
            {
                result[i] = B[i - A.Length];
            }
            return result;
        }

        /// <summary>
        /// Skalar x Vektor Multiplikation
        /// </summary>
        /// <param name="A">Skalar A</param>
        /// <param name="B">Matrix B</param>
        /// <returns>Matrix C = a*B</returns>
        public static Matrix operator *(double a, Vector B)
        {
            Vector result = new Vector(B.Length);
            for (int i = 0; i < B.Length; i++)
                    result[i] = a * B[i];
            return result;
        }

        /// <summary>
        /// Skalar x Vektor Multiplikation
        /// </summary>
        /// <param name="A">Skalar A</param>
        /// <param name="B">Vektor B</param>
        /// <returns>Vektor C = a*B</returns>
        public static Matrix operator *(Vector B, double a)
        {
            Vector result = new Vector(B.Length);
            for (int i = 0; i < B.Length; i++)
                result[i] = a * B[i];
            return result;
        }

        //Zugriff auf Vektor Element
        public double this[int x]
        {
            get { return base[x,0]; }
            set { base[x, 0] = value; }
        }

        public int ContainsValue(double x)
        {
            for(int i = 0; i < this.Length; i++)
            {
                if (GeneralHelper.isXAlmostEqualToY(this[i], x))
                    return i;
            }
            return -1;
        }

        public double GetMaxAbsValue()
        {
            double tempMax = Math.Abs(this[0]);
            for(int i = 0; i < this.Length; i++)
            {
                if(tempMax < Math.Abs(this[0]))
                {
                    tempMax = Math.Abs(this[0]);
                }
            }
            return tempMax;
        }

        //Gibt die Länge des Vektors wieder
        public int Length{
            get { return base.NoRows; }
        }
    }
}
