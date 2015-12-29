using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._1DSystem
{
    public class DGSystemElement
    {
        public DGSystemElement LeftNeighbour, RightNeighbour;
        public double leftSpaceBoundary, rightSpaceBoundary;
        public Vector nodes, integrationWeights;

        LagrangeInterpolator interpolator;
        int SystemDimension;
        int N; //Polynomgrad

        Matrix SolutionSystem;
        Matrix FluxEvaluation;

        Matrix MassMatrix;
        Matrix InvMassMatrix;
        Matrix DifferentialMatrix;
        Matrix B;

        Matrix VolumeMatrix;
        Matrix SurfaceMatrix;

        Func<Vector, Vector, Vector> NumericalFluxFunction;
        Func<Vector, double, double, Vector> InhomogenuousFunction;
        Func<Vector, Vector> FluxFunction;
        Func<Vector, Vector> InitialFunction;

        public Vector LeftBoarderValue { get; set; }
        public Vector RightBoarderValue { get; set; }


        public DGSystemElement(double leftSpaceBoundary, double rightSpaceBoundary, int polynomOrder,int systemDimension, Func<Vector, Vector, Vector> NumericalFluxFunction, Func<Vector,double, double, Vector> InhomogenuousFunction, Func<Vector, Vector> FluxFunction, Func<Vector, Vector> InitialFunction)
        {
            this.leftSpaceBoundary = leftSpaceBoundary;
            this.rightSpaceBoundary = rightSpaceBoundary;
            this.N = polynomOrder;
            this.SystemDimension = systemDimension;
            this.FluxFunction = FluxFunction;
            this.NumericalFluxFunction = NumericalFluxFunction;
            this.InhomogenuousFunction = InhomogenuousFunction;
            this.InitialFunction = InitialFunction;

            Initialize();
        }

        public void Initialize() {
            FluxEvaluation = new Matrix(N + 1, SystemDimension);
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out nodes, out integrationWeights);
            SolutionSystem = ComputeStartSolution();
            LeftBoarderValue = new Vector(SystemDimension);
            RightBoarderValue = new Vector(SystemDimension);

            this.UpdateBorderValues();

            interpolator = new LagrangeInterpolator(nodes);
            MassMatrix = IntegrationToolbox.generateMassMatrix(integrationWeights);


            InvMassMatrix = new Matrix(N + 1, N + 1);
            for (int i = 0; i < N + 1; i++)
                InvMassMatrix[i, i] = 1.0 / MassMatrix[i, i];

            DifferentialMatrix = interpolator.computeLagrangePolynomeDerivativeMatrix();
            B = new Matrix(N + 1, N + 1);
            B[0, 0] = -1.0;
            B[N, N] = 1.0;

            VolumeMatrix = InvMassMatrix * (!DifferentialMatrix) * MassMatrix;
            SurfaceMatrix = InvMassMatrix * B;
        }

        public Matrix EvaluateTimeDerivativeGaussLobatto(double time)
        {

            Matrix SystemTimeDerivative = new Matrix(N + 1, SystemDimension);

            Vector NumFluxLeft = NumericalFluxFunction(LeftNeighbour.RightBoarderValue, this.LeftBoarderValue);
            Vector NumFluxRight = NumericalFluxFunction(this.RightBoarderValue, RightNeighbour.LeftBoarderValue);


            for(int i = 0; i < N + 1; i++)
            {
                Vector SolutionAt = SolutionSystem.GetRowAsColumnVector(i);
                Vector EvaluatedFlux = FluxFunction(SolutionAt);
                FluxEvaluation.InjectMatrixAtPosition(!EvaluatedFlux, i, 0);
            }

            Matrix Inhomogenuous = ComputeInhomogenuousPart(time);

            //WEAK
            //for (int sysIdx = 0; sysIdx < SystemDimension; sysIdx++)
            //{
            //    Vector timeDerivate;
            //    Vector volume = new Vector(N + 1);
            //    Vector surface = new Vector(N + 1);
            //    surface[0] = NumFluxLeft[sysIdx];
            //    surface[N] = NumFluxRight[sysIdx];

            //    for (int vol = 0; vol < N + 1; vol++)
            //    {
            //        volume[vol] = FluxEvaluation[vol, sysIdx];
            //    }
            //    timeDerivate = VolumeMatrix * volume - SurfaceMatrix * surface;
            //    timeDerivate = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivate);
            //    timeDerivate += Inhomogenuous.GetColumn(sysIdx);
            //    SystemTimeDerivative.InjectMatrixAtPosition(timeDerivate, 0, sysIdx);
            //}

            //Strong
            for (int sysIdx = 0; sysIdx < SystemDimension; sysIdx++)
            {
                Vector timeDerivate;
                Vector volume = new Vector(N + 1);
                Vector surface = new Vector(N + 1);
                surface[0] = NumFluxLeft[sysIdx];
                surface[N] = NumFluxRight[sysIdx];

                for (int vol = 0; vol < N + 1; vol++)
                {
                    volume[vol] = FluxEvaluation[vol, sysIdx];
                }

                surface = surface - volume;

                timeDerivate = -1.0 * InvMassMatrix * B * surface - DifferentialMatrix * volume;
                timeDerivate = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivate);
                timeDerivate += Inhomogenuous.GetColumn(sysIdx);
                
                SystemTimeDerivative.InjectMatrixAtPosition(timeDerivate, 0, sysIdx);
            }

            return SystemTimeDerivative;
        }

        public void UpdateSolutionSystem(Matrix UpdatedSolution)
        {
            this.SolutionSystem = UpdatedSolution;
        }


        public Matrix GetSolution()
        {
            return SolutionSystem;
        }
        private Matrix ComputeInhomogenuousPart(double time)
        {
            Matrix InhomogMatrix = new Matrix(N + 1, SystemDimension);
            Vector evaluationNodes = new Vector(SystemDimension);
            Vector solutionAtNode = new Vector(SystemDimension);
            
            for(int i = 0; i < N + 1; i++)
            {
                for (int m = 0; m < SystemDimension; m++)
                {
                    evaluationNodes[m] = MapToOriginSpace(nodes[i]);                  
                }
                solutionAtNode = SolutionSystem.GetRowAsColumnVector(i);
                Vector res = InhomogenuousFunction(solutionAtNode,evaluationNodes[0], time);
                InhomogMatrix.InjectMatrixAtPosition(!res, i, 0);
            }
            return InhomogMatrix;
        }

        private Matrix ComputeStartSolution()
        {
            Matrix startSolution = new Matrix(N + 1,SystemDimension);

            for (int i = 0; i <= N; i++)
            {
                Vector evaluationNodes = new Vector(SystemDimension);
                for (int m = 0; m < SystemDimension; m++)
                {
                    evaluationNodes[m] = MapToOriginSpace(nodes[i]);
                }

                Vector evaluation = InitialFunction(evaluationNodes);

                for (int k = 0; k < SystemDimension; k++)
                    startSolution[i, k] = evaluation[k];
            }
            return startSolution;
        }

        public Vector GetOriginNodes()
        {
            Vector originNodes = new Vector(nodes.Length);
            for (int i = 0; i < nodes.Length; i++)
                originNodes[i] = MapToOriginSpace(nodes[i]);
            return originNodes;
        }

        public double MapToOriginSpace(double x)
        {
            return leftSpaceBoundary + ((x + 1.0) / 2.0) * (rightSpaceBoundary - leftSpaceBoundary);
        }

        private double MapToReferenceSpace(double x)
        {
            return (2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * (x - leftSpaceBoundary) - 1.0;
        }

        public double GetSpaceLength()
        {
            return rightSpaceBoundary - leftSpaceBoundary;
        }

        public void UpdateBorderValues()
        {
            for (int i = 0; i < SystemDimension; i++)
            {
                LeftBoarderValue[i] = SolutionSystem[0, i];
                RightBoarderValue[i] = SolutionSystem[N, i];
            }
        }
    }
}
