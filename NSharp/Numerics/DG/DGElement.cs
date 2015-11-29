using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    class DGElement
    {
        public DGElement LeftNeighbour, RightNeighbour;
        double leftSpaceBoundary, rightSpaceBoundary;
        double spaceTransformationFactor;
        Vector gaussLobattoNodes, gaussLobattoWeights;

        LagrangeInterpolator interpolator;
        int N; //Polynom order
        Vector solution;
        Vector timeDerivative;
        Vector fluxEvaluation;

        Matrix massMatrix;
        Matrix invMassMatrix;
        Matrix differentialMatrix;
        Matrix B;

        double numFluxLeft, numFluxRight;

        Func<double, double> leftBoundaryFunction;
        Func<double, double> rightBoundaryFunction;
        Func<double, double, double> numericalFluxFunction;
        Func<double, double, double> inhomogenuousFunction;
        Func<double, double> fluxFunction;


        public double LeftBorderValue{ get; set; }
        public double RightBorderValue { get; set; }

        public DGElement(double leftSpaceBoundary,double rightSpaceBoundary, int N, Func<double,double> fluxFunction, Func<double, double> leftBoundaryFunction, Func<double, double> rightBoundaryFunction, Func<double, double, double> numericalFluxFunction, Func<double, double, double> inhomogenuousFunction)
        {
            this.leftBoundaryFunction = leftBoundaryFunction;
            this.rightBoundaryFunction = rightBoundaryFunction;
            this.numericalFluxFunction = numericalFluxFunction;
            this.inhomogenuousFunction = inhomogenuousFunction;
            this.fluxFunction = fluxFunction;
            this.leftSpaceBoundary = leftSpaceBoundary;
            this.rightSpaceBoundary = rightSpaceBoundary;
            this.N = N;
        }

        public void Initialize(Vector startSolution)
        {
            solution = startSolution;
            timeDerivative = new Vector(N + 1);
            fluxEvaluation = new Vector(N + 1);
            
            LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(N, out gaussLobattoNodes, out gaussLobattoWeights);

            interpolator = new LagrangeInterpolator(gaussLobattoNodes);

            massMatrix = IntegrationToolbox.generateMassMatrix(gaussLobattoWeights);

            invMassMatrix = new Matrix(N + 1, N+1);
            for (int i = 0; i < N + 1; i++)
                invMassMatrix[i,i] = 1.0 / massMatrix[i,i];

            differentialMatrix = interpolator.computeLagrangePolynomeDerivativeMatrix();
            
            //Randauswertung
            B = new Matrix(N + 1, N + 1);
            B[0, 0] = -1.0;
            B[N, N] = 1.0;
        }

        public void updateSolution(double t, double timeStep)
        {
            List<Func<double, double, double>> mySystem = new List<Func<double, double, double>>();
            mySystem.Add(timeDerivativeMapping);
            OrdinaryPartialEquationsSolver.OrdinaryDifferentialEquation myODE = new OrdinaryPartialEquationsSolver.OrdinaryDifferentialEquation(mySystem);

            OrdinaryPartialEquationsSolver.IODESolver odeSolver = new OrdinaryPartialEquationsSolver.RungeKuttaSolver();
            solution = odeSolver.computeSolutionForNextStep(solution, myODE, t, t + timeStep);
        }

        private void evaluateTimeDerivative(double time)
        {           
            //Auswertung an numerischen Flüssen bestimmen
            if (LeftNeighbour == null)
            {
                double left = leftBoundaryFunction(time);
                numFluxLeft = numericalFluxFunction(left, solution[0]);
            }
            else
            {
                numFluxLeft = numericalFluxFunction(LeftNeighbour.RightBorderValue, solution[0]);
            }

            if (RightNeighbour == null)
            {
                numFluxRight = numericalFluxFunction(solution[gaussLobattoNodes.Length - 1], rightBoundaryFunction(time));
            }
            else
            {
                double right = solution[gaussLobattoNodes.Length - 1];
                numFluxRight = numericalFluxFunction(right, RightNeighbour.LeftBorderValue);
            }


            for (int i = 0; i < gaussLobattoNodes.Length; i++)
            {
                fluxEvaluation[i] = fluxFunction(solution[i]);
            }
            Vector volume = new Vector(N+1);
            Vector surface = new Vector(N+1);

            volume[0] = fluxEvaluation[0];//numFluxLeft;
            surface[0] = numFluxLeft;

            volume[N] = fluxEvaluation[N];// numFluxRight;
            surface[N] = numFluxRight;


            for (int i = 1; i < N ; i++)
                volume[i] = fluxEvaluation[i];

            timeDerivative = invMassMatrix * ((((!differentialMatrix) * massMatrix) * volume) - (B * surface));
            timeDerivative = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivative);
            timeDerivative += computeInhomogenuousPart(time);
        }

        private Vector computeInhomogenuousPart(double time)
        {
            Vector inhomoPart = new Vector(N + 1);
            for (int i = 0; i < N + 1; i++)
                inhomoPart[i] = inhomogenuousFunction(gaussLobattoNodes[i], time);

            return inhomoPart;
        }

        public Vector getOriginNodes()
        {
            Vector originNodes = new Vector(gaussLobattoNodes.Length);
            for (int i = 0; i < gaussLobattoNodes.Length; i++)
                originNodes[i] = MapToOriginSpace(gaussLobattoNodes[i]);
            return originNodes;
        }

        private double MapToOriginSpace(double x)
        {
            return leftSpaceBoundary + ((x + 1.0) / 2.0) * (rightSpaceBoundary - leftSpaceBoundary);
        }

        private double MapToReferenceSpace(double x)
        {
            return (2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * (x - leftSpaceBoundary) - 1.0;
        }

        private double timeDerivativeMapping(double u, double time)
        {
            return u;
        }
    }
}
