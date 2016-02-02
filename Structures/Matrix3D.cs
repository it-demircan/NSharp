using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Structures
{
    public class Matrix3D
    {
        Matrix[] matrices = null;
        public int DimZ{get;set;}
        public int DimX { get; set; }
        public int DimY { get; set; }


        public Matrix3D(int dimN, int dimM, int dimZ)
        {
            matrices = new Matrix[dimZ];
            this.DimZ = dimZ;
            this.DimX = dimN;
            this.DimY = dimM;
            for (int i = 0; i < dimZ; i++)
            {
                matrices[i] = new Matrix(dimN, dimM);
            }
        }

        public Vector this[int x, int y]
        {
            get { return GetZValues(x, y); }
            set { SetZValues(x, y,value); }
        }

        public double this[int x, int y, int z]
        {
            get { return matrices[z][x, y]; }
            set { matrices[z][x, y] = value; }
        }

        private Vector GetZValues(int x, int y)
        {
            Vector zValues = new Vector(DimZ);

            for(int i = 0; i < DimZ; i++)
            {
                zValues[i] = matrices[i][x, y];
            }

            return zValues;
        }

        private void SetZValues(int x, int y, Vector zValues)
        {
            for (int i = 0; i < DimZ; i++)
            {
                matrices[i][x, y] = zValues[i];
            }
        }

        public Matrix this[int z]
        {
            get { return matrices[z]; }
            set { matrices[z] = value; }
        }
    }
}
