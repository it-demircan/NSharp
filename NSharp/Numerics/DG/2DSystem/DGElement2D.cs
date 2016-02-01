using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._2DSystem
{
    class DGElement2D
    {
        public Matrix3D Solution { get; set; }
        public Matrix3D FluxF { get; set; }
        public Matrix3D FluxG { get; set; }

        public int N { get; set; }
        public int SystemDimension { get; set; }

        public double LeftXBorder{get;set;}
        public double RightXBorder{ get; set; }

        public double TopYBorder { get; set; }
        public double BottomYBorder { get; set; }

        private Vector Nodes;
        private Vector IntegrationWeights;

        private double J { get; set; }
        public Matrix3D BottomBorderValues { get; set; }
        public Matrix3D LeftBorderValues { get; set; }
        public Matrix3D RightBorderValues { get; set; }
        public Matrix3D TopBorderValues { get; set; }


        Func<Vector, Vector> InitialFunction;
        Func<Vector, Vector> FluxF_Function;
        Func<Vector, Vector> FluxG_Function;
        Func<Vector,Vector, Vector> NumFluxFunctionF;
        Func<Vector, Vector, Vector> NumFluxFunctionG;


        Matrix DifferentialMatrix, DifferentialMatrixTransposed;
        Matrix S;
        LagrangeInterpolator interpolator;


        private void InitializeDGElement()
        {
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out Nodes, out IntegrationWeights);
            InitializeStartSolution();
            BottomBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            TopBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            LeftBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            RightBorderValues = new Matrix3D(N + 1, 1, SystemDimension);

            this.UpdateBorderValues();

            S = new Matrix(N + 1, N + 1);
            S[0, 0] = 1.0 / IntegrationWeights[0];
            S[N, N] = -1.0 / IntegrationWeights[N];

            interpolator = new LagrangeInterpolator(Nodes);
            DifferentialMatrix = interpolator.computeLagrangePolynomeDerivativeMatrix();
            DifferentialMatrixTransposed = !DifferentialMatrix;

            J = ((RightXBorder - LeftXBorder) * (TopYBorder - BottomYBorder)) / 4.0;
        }


        public void InitializeStartSolution()
        {
            Matrix3D originNodes = GetOriginNodes();

            for(int i = 0; i < originNodes.DimX; i++)
            {
                for(int k = 0; k < originNodes.DimY; k++)
                {
                    Solution[i, k] = InitialFunction(originNodes[i, k]);
                }
            }
        }


        public DGElement2D Left, Right, Top, Bottom;
        public Matrix3D ComputeNumericalFluxFAndScale()
        {
            Matrix3D NumFluxF = new Matrix3D(N + 1, N + 1, SystemDimension);

            for(int i = 0; i < N+1; i++)
            {
                NumFluxF[0, i] = (Vector) (((TopYBorder-BottomYBorder)/2.0) * NumFluxFunctionF(this.LeftBorderValues[i, 0], Left.RightBorderValues[i, 0]));
                NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], Right.LeftBorderValues[i, 0]));
            }

            return NumFluxF;
        }

        public Matrix3D ComputeNumericalFluxGAndScale()
        {
            Matrix3D NumFluxG = new Matrix3D(N + 1, N + 1, SystemDimension);

            for (int i = 0; i < N + 1; i++)
            {
                NumFluxG[i,0] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.BottomBorderValues[i, 0], Bottom.TopBorderValues[i, 0]));
                NumFluxG[i,N] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.TopBorderValues[i, 0], Top.BottomBorderValues[i, 0]));
            }
            return NumFluxG;
        }

        public void UpdateFluxesAndScale()
        {
            for(int i = 0; i <= N; i++)
            {
                for(int k = 0; k <= N; k++)
                {
                    FluxF[i, k] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) *   FluxF_Function(Solution[i, k]));
                    FluxG[i, k] = (Vector)(((RightXBorder-LeftXBorder)/2.0) *       FluxG_Function(Solution[i, k]));
                }
            }
        }

        public void UpdateSolution(Matrix3D solution)
        {
            this.Solution = solution;
        }

        public void UpdateBorderValues()
        {
            for(int i = 0; i< Nodes.Length; i++)
            {
                BottomBorderValues[i, 0] = Solution[i, 0];
                TopBorderValues[i, 0] = Solution[i, N];
                LeftBorderValues[i, 0] = Solution[0, i];
                RightBorderValues[i, 0] = Solution[N, i];
            }
        }

        public Matrix3D GetOriginNodes()
        {
            Matrix3D originNodes = new Matrix3D(Nodes.Length, Nodes.Length, 2);
            for (int x = 0; x < Nodes.Length; x++)
                for (int y = 0; y < Nodes.Length; y++)
                {
                    originNodes[x, y, 0] = MapToOriginX(Nodes[x]);
                    originNodes[x, y, 1] = MapToOriginY(Nodes[y]);
                }
            return originNodes;
        }

        public double MapToOriginX(double x)
        {
            return LeftXBorder + ((x + 1.0) / 2.0) * (RightXBorder - LeftXBorder);
        }

        public double MapToOriginY(double y)
        {
            return BottomYBorder + ((y + 1.0) / 2.0) * (TopYBorder - BottomYBorder);
        }
    }
}
