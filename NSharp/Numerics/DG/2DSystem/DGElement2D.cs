using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._2DSystem
{
    public class DGElement2D
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
        Func<Vector,double, Vector> ExactFunction;


        Matrix DifferentialMatrix, DifferentialMatrixTransposed;
        Matrix S;
        LagrangeInterpolator interpolator;
        public DGElement2D Left, Right, Top, Bottom;



        public DGElement2D(int N, int sysDim, double leftX, double rightX, double bottomY, double topY,
                Func<Vector, Vector> InitialFunction,
                Func<Vector, Vector> FluxF_Function,
                Func<Vector, Vector> FluxG_Function,
                Func<Vector, Vector, Vector> NumFluxFunctionF,
                Func<Vector, Vector, Vector> NumFluxFunctionG,
                Func<Vector,double, Vector> ExactFunction = null
            )
        {
            this.N = N;
            this.SystemDimension = sysDim;
            this.LeftXBorder = leftX;
            this.RightXBorder = rightX;
            this.TopYBorder = topY;
            this.BottomYBorder = bottomY;

            this.InitialFunction = InitialFunction;
            this.FluxF_Function = FluxF_Function;
            this.FluxG_Function = FluxG_Function;
            this.NumFluxFunctionF = NumFluxFunctionF;
            this.NumFluxFunctionG = NumFluxFunctionG;
            this.ExactFunction = ExactFunction;
            InitializeDGElement();
        }

        private void InitializeDGElement()
        {
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out Nodes, out IntegrationWeights);
            InitializeStartSolution();

            BottomBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            TopBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            LeftBorderValues = new Matrix3D(N + 1, 1, SystemDimension);
            RightBorderValues = new Matrix3D(N + 1, 1, SystemDimension);

            this.FluxF = new Matrix3D(N + 1, N + 1, SystemDimension);
            this.FluxG = new Matrix3D(N + 1, N + 1, SystemDimension);

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
            Solution = new Matrix3D(N + 1, N + 1, SystemDimension);

            for(int i = 0; i < originNodes.DimX; i++)
            {
                for(int k = 0; k < originNodes.DimY; k++)
                {
                    Solution[i, k] = InitialFunction(originNodes[i, k]);
                }
            }
        }


        public Matrix3D ComputeTimeDerivative(double time)
        {
            Matrix3D TimeDerivative = new Matrix3D(N + 1, N + 1, SystemDimension);

            Matrix3D NumFluxFEvaluated = ComputeNumericalFluxFAndScale();

            Matrix3D NumFluxGEvaluated = ComputeNumericalFluxGAndScale();

            UpdateFluxesAndScale();

            for (int i = 0; i < SystemDimension; i++)
            {
                Matrix tempResult = new Matrix(N + 1, N + 1);

                tempResult = S * (NumFluxFEvaluated[i] - FluxF[i])
                            + (NumFluxGEvaluated[i] - FluxG[i]) * S
                            - FluxG[i] * DifferentialMatrixTransposed
                            - DifferentialMatrix * FluxF[i];

                tempResult = (1.0 / J) * tempResult;
                TimeDerivative[i] = tempResult;
            }

            return TimeDerivative;
        }


        public Matrix3D ComputeNumericalFluxFAndScale()
        {
            Matrix3D NumFluxF = new Matrix3D(N + 1, N + 1, SystemDimension);

            for(int i = 0; i < N+1; i++)
            {
                NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(Left.RightBorderValues[i, 0], this.LeftBorderValues[i, 0]));
                NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], Right.LeftBorderValues[i, 0]));
            }

            return NumFluxF;
        }

        public Matrix3D ComputeNumericalFluxFAndScaleTaskThree()
        {
            Matrix3D NumFluxF = new Matrix3D(N + 1, N + 1, SystemDimension);

            for (int i = 0; i < N + 1; i++)
            {
                Vector temp;
                //Reflektierend
                if (this.LeftXBorder == 0.3 && this.TopYBorder <= 0.4)
                {
                    temp = this.LeftBorderValues[i, 0];
                    temp[1] = -1.0 * temp[1];
                    NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(temp, this.LeftBorderValues[i, 0]));
                    NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], Right.LeftBorderValues[i, 0]));
                }
                else if (this.RightXBorder == 0.7 && this.TopYBorder <= 0.4)
                {
                    temp = this.RightBorderValues[i, 0];
                    temp[1] = -1.0 * temp[1];

                    NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], temp));
                    NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(Left.RightBorderValues[i, 0], this.LeftBorderValues[i, 0]));

                }
                else if(this.RightXBorder == 1.0)
                {
                    temp = new Vector(3);
                    temp[2] = 2.0;
                    NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], temp));
                    NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(Left.RightBorderValues[i, 0], this.LeftBorderValues[i, 0]));
                }
                else if(this.LeftXBorder == 0.0)
                {
                    temp = new Vector(3);
                    temp[2] = 2.0;
                    NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(temp, this.LeftBorderValues[i, 0]));
                    NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], Right.LeftBorderValues[i, 0]));
                }
                else {
                    NumFluxF[0, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(Left.RightBorderValues[i, 0], this.LeftBorderValues[i, 0]));
                    NumFluxF[N, i] = (Vector)(((TopYBorder - BottomYBorder) / 2.0) * NumFluxFunctionF(this.RightBorderValues[i, 0], Right.LeftBorderValues[i, 0]));
                }
            }

            return NumFluxF;
        }

        public Matrix3D ComputeNumericalFluxGAndScaleTaskThree()
        {
            Matrix3D NumFluxG = new Matrix3D(N + 1, N + 1, SystemDimension);

            for (int i = 0; i < N + 1; i++)
            {
                Vector temp;
                if (this.TopYBorder == 1.0)
                {
                    temp = new Vector(3);
                    temp[2] = 2.0;

                    NumFluxG[i, N] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.TopBorderValues[i, 0], temp));
                    NumFluxG[i, 0] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(Bottom.TopBorderValues[i, 0], this.BottomBorderValues[i, 0]));

                }               
                else if(this.BottomYBorder == 0.0 || (this.BottomYBorder == 0.4 && (this.RightXBorder <= 0.3 || this.LeftXBorder >= 0.7)))
                {
                    temp = this.BottomBorderValues[i, 0];
                    temp[0] = -1.0 * temp[0];

                    NumFluxG[i, 0] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(temp, this.BottomBorderValues[i, 0]));
                    NumFluxG[i, N] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.TopBorderValues[i, 0], Top.BottomBorderValues[i, 0]));
                }
                else
                {
                    NumFluxG[i, 0] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(Bottom.TopBorderValues[i, 0], this.BottomBorderValues[i, 0]));
                    NumFluxG[i, N] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.TopBorderValues[i, 0], Top.BottomBorderValues[i, 0]));
                }
            }
            return NumFluxG;
        }

        public Matrix3D ComputeNumericalFluxGAndScale()
        {
            Matrix3D NumFluxG = new Matrix3D(N + 1, N + 1, SystemDimension);

            for (int i = 0; i < N + 1; i++)
            {
                NumFluxG[i,0] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(Bottom.TopBorderValues[i, 0],this.BottomBorderValues[i, 0]));
                NumFluxG[i,N] = (Vector)(((RightXBorder - LeftXBorder) / 2.0) * NumFluxFunctionG(this.TopBorderValues[i, 0], Top.BottomBorderValues[i, 0]));
            }
            return NumFluxG;
        }


        public Matrix3D EvaluateExactSolution(double time)
        {
            Matrix3D originNodes = GetOriginNodes();
            Matrix3D eva = new Matrix3D(N + 1, N + 1, SystemDimension);

            for (int i = 0; i < originNodes.DimX; i++)
            {
                for (int k = 0; k < originNodes.DimY; k++)
                {
                    eva[i, k] = ExactFunction(originNodes[i, k], time);
                }
            }

            return eva;
        }

        public double GetMaxErrorAtTimeForSystemEquation(double time, int sysIdx)
        {
            Matrix approx = Solution[sysIdx];
            Matrix3D ExactAll = EvaluateExactSolution(time);
            Matrix exact = ExactAll[sysIdx];

            Matrix diff = exact - approx;
            double maxError = Math.Abs(diff[0, 0]);

            for(int i = 0;i< diff.NoRows; i++)
                for(int k = 0; k < diff.NoColumns; k++)
                {
                    if (maxError < Math.Abs(diff[i, k]))
                        maxError = Math.Abs(diff[i, k]);
                }

            return maxError;
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
