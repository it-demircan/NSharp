using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public enum IntegrationMode
    {
        GaussLegendre, GaussLobatto
    }
    class DGElement
    {
        public DGElement LeftNeighbour, RightNeighbour;
        public double leftSpaceBoundary, rightSpaceBoundary;
        public Vector nodes, integrationWeights;
        IntegrationMode myIntegrationMode;

        LagrangeInterpolator interpolator;
        int N; //Polynomgrad
        Vector solution;
        Vector fluxEvaluation;

        Matrix massMatrix;
        Matrix invMassMatrix;
        Matrix differentialMatrix;
        Matrix B;
        Vector LeftBoarderInterpolation;
        Vector RightBoarderInterpolation;

        Matrix volumeMatrix;
        Matrix surfaceMatrix;

        double numFluxLeft, numFluxRight;

        Func<double, double, double> numericalFluxFunction;
        Func<double, double, double> inhomogenuousFunction;
        Func<double, double> fluxFunction;
        Func<double, double> initialFunction;


        public double LeftBorderValue{ get; set; }
        public double RightBorderValue { get; set; }

        public DGElement(IntegrationMode mode,double leftSpaceBoundary,double rightSpaceBoundary, int polynomOrder, Func<double,double> fluxFunction, Func<double, double, double> numericalFluxFunction, Func<double, double, double> inhomogenuousFunction, Func<double, double> initialFunction)
        {
            myIntegrationMode = mode;
            this.leftSpaceBoundary = leftSpaceBoundary;
            this.rightSpaceBoundary = rightSpaceBoundary;
            this.N = polynomOrder;
            this.fluxFunction = fluxFunction;
            this.numericalFluxFunction = numericalFluxFunction;
            this.inhomogenuousFunction = inhomogenuousFunction;
            this.initialFunction = initialFunction;

            Initialize();
        }

        private void Initialize()
        {
            if (myIntegrationMode == IntegrationMode.GaussLobatto)
                InitializeDGWithGaussLobattoNodes();
            else if (myIntegrationMode == IntegrationMode.GaussLegendre)
                InitializeDGWithGaussNodes();
        }

        public Vector EvaluateTimeDerivative(double time)
        {
            if (myIntegrationMode == IntegrationMode.GaussLobatto)
                return EvaluateTimeDerivativeGaussLobatto(time);
            else if (myIntegrationMode == IntegrationMode.GaussLegendre)
                return EvaluateTimeDerivativeGaussLegendre(time);
            return null;
        }

        private void InitializeDGWithGaussLobattoNodes()
        {
            fluxEvaluation = new Vector(N + 1);
            
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out nodes, out integrationWeights);

            solution = ComputeStartSolution() ;
            this.UpdateBorderValues();

            interpolator = new LagrangeInterpolator(nodes);
            massMatrix = IntegrationToolbox.generateMassMatrix(integrationWeights);

            invMassMatrix = new Matrix(N + 1, N+1);
            for (int i = 0; i < N + 1; i++)
                invMassMatrix[i,i] = 1.0 / massMatrix[i,i];

            differentialMatrix = interpolator.computeLagrangePolynomeDerivativeMatrix();
            B = new Matrix(N + 1, N + 1);
            B[0, 0] = -1.0;
            B[N, N] = 1.0;

            volumeMatrix = invMassMatrix*(!differentialMatrix) * massMatrix;
            surfaceMatrix = invMassMatrix * B;    
        }

        private void InitializeDGWithGaussNodes()
        {
            fluxEvaluation = new Vector(N + 1);

            LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(N, out nodes, out integrationWeights);

            solution = ComputeStartSolution();
           
            interpolator = new LagrangeInterpolator(nodes);
            massMatrix = IntegrationToolbox.generateMassMatrix(integrationWeights);
            this.UpdateBorderValues();

            invMassMatrix = new Matrix(N + 1, N + 1);
            for (int i = 0; i < N + 1; i++)
                invMassMatrix[i, i] = 1.0 / massMatrix[i, i];

            differentialMatrix = interpolator.computeLagrangePolynomeDerivativeMatrix();

            LeftBoarderInterpolation = new Vector(nodes.Length);
            RightBoarderInterpolation = new Vector(nodes.Length);

            for (int i = 0; i < nodes.Length; i++)
            {
                LeftBoarderInterpolation[i] = interpolator.evaluateLagrangePolynome(-1.0, i);
                RightBoarderInterpolation[i] = interpolator.evaluateLagrangePolynome(1.0, i);
            }

            B = new Matrix(N + 1, N + 1);
            for (int i = 0; i <= N; i++)
            {
                B[i, 0] = LeftBoarderInterpolation[i];
                B[i, N] = RightBoarderInterpolation[i];
            }

            volumeMatrix = invMassMatrix * (!differentialMatrix) * massMatrix;
            LeftBoarderInterpolation = invMassMatrix * LeftBoarderInterpolation;
            RightBoarderInterpolation = invMassMatrix * RightBoarderInterpolation;
            
        }

        public Vector EvaluateTimeDerivativeGaussLegendre(double time)
        {
            Vector timeDerivative;      
            numFluxLeft = numericalFluxFunction(LeftNeighbour.RightBorderValue, this.LeftBorderValue);
            numFluxRight = numericalFluxFunction(this.RightBorderValue, RightNeighbour.LeftBorderValue);

            Vector volume = new Vector(N+1);
            Vector surface = new Vector(N+1);

            surface[0] = numFluxLeft;
            surface[N] = numFluxRight;

            for (int i = 0; i <= N ; i++)
                volume[i] = fluxFunction(solution[i]);

            timeDerivative = volumeMatrix * volume-(Vector)(numFluxRight*RightBoarderInterpolation) + (Vector)(numFluxLeft*LeftBoarderInterpolation);
            timeDerivative = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivative);
            timeDerivative += ComputeInhomogenuousPart(time);
            return timeDerivative;
        }

        public Vector EvaluateTimeDerivativeGaussLobatto(double time)
        {     
            Vector timeDerivative;
            numFluxLeft = numericalFluxFunction(LeftNeighbour.RightBorderValue, this.LeftBorderValue);
            numFluxRight = numericalFluxFunction(this.RightBorderValue, RightNeighbour.LeftBorderValue);

            Vector volume = new Vector(N + 1);
            Vector surface = new Vector(N + 1);

            surface[0] = numFluxLeft;
            surface[N] = numFluxRight;

            for (int i = 0; i <= N; i++)
                volume[i] = fluxFunction(solution[i]);

            timeDerivative = volumeMatrix * volume - surfaceMatrix * surface;

            timeDerivative = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivative);
            timeDerivative += ComputeInhomogenuousPart(time);
            return timeDerivative;
        }



        public void UpdateSolution(Vector solution)
        {
            this.solution = solution;
        }

        public void UpdateBorderValues()
        {
            if (myIntegrationMode == IntegrationMode.GaussLobatto)
            {
                this.LeftBorderValue = solution[0];
                this.RightBorderValue = solution[N];
            }
            else
            {
                this.LeftBorderValue = interpolator.evaluateLagrangeRepresentation(-1.0, solution);
                this.RightBorderValue = interpolator.evaluateLagrangeRepresentation(1.0, solution);
            }
        }

        private Vector ComputeInhomogenuousPart(double time)
        {
            Vector inhomoPart = new Vector(N + 1);
            for (int i = 0; i < N + 1; i++)
                inhomoPart[i] = inhomogenuousFunction(nodes[i], time);

            return inhomoPart;
        }

        //Gibt M^-1 * B zurück
        public Matrix GetSMatrix()
        {
            return invMassMatrix * B;
        }

        public Matrix GetDifferentialMatrix()
        {
            return differentialMatrix;
        }

        public Matrix GetL1Matrix()
        {
            Matrix L1 = new Matrix(N + 1, N + 1);
            for (int i = 0; i <= N; i++)
                L1[0, i] = interpolator.evaluateLagrangePolynome(1.0, i);
            return L1;
        }

        public Matrix GetL2Matrix()
        {
            Matrix L2 = new Matrix(N + 1, N + 1);
            for (int i = 0; i <= N; i++)
                L2[N, i] = interpolator.evaluateLagrangePolynome(1.0, i);
            return L2;
        }

        public Matrix GetBMatrix()
        {
            return B;
        }

        public Matrix GetMassMatrix()
        {
            return massMatrix;
        }
        public Matrix GetInverseMassMatrix()
        {
            return invMassMatrix;
        }


        private Vector ComputeStartSolution()
        {
            Vector startSolution = new Vector(N + 1);
            for(int i = 0; i <= N; i++)
            {
                startSolution[i] = initialFunction(MapToOriginSpace(nodes[i]));
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

        public Vector GetSolution()
        {
            return solution;
        }
        /**
            Periodic Boundary
        **/

        public double GetSpaceLength()
        {
            return rightSpaceBoundary - leftSpaceBoundary;
        }



        private double LeftPeriodicBoundary(double time)
        {
            return LeftNeighbour.RightBorderValue;
        }

        private double RightPeriodicBoundary(double time)
        {
            return RightNeighbour.LeftBorderValue;
        }
    }
}
