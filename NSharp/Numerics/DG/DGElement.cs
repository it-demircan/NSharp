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
        Vector nodes, integrationWeights;

        LagrangeInterpolator interpolator;
        int N; //Polynomgrad
        Vector solution;
        Vector fluxEvaluation;

        Matrix massMatrix;
        Matrix invMassMatrix;
        Matrix differentialMatrix;
        Matrix B;

        Matrix volumeMatrix;
        Matrix surfaceMatrix;

        double numFluxLeft, numFluxRight;

        Func<double, double, double> numericalFluxFunction;
        Func<double, double, double> inhomogenuousFunction;
        Func<double, double> fluxFunction;
        Func<double, double> initialFunction;


        public double LeftBorderValue{ get; set; }
        public double RightBorderValue { get; set; }

        public DGElement(double leftSpaceBoundary,double rightSpaceBoundary, int polynomOrder, Func<double,double> fluxFunction, Func<double, double, double> numericalFluxFunction, Func<double, double, double> inhomogenuousFunction, Func<double, double> initialFunction)
        {
            this.leftSpaceBoundary = leftSpaceBoundary;
            this.rightSpaceBoundary = rightSpaceBoundary;
            this.N = polynomOrder;
            this.fluxFunction = fluxFunction;
            this.numericalFluxFunction = numericalFluxFunction;
            this.inhomogenuousFunction = inhomogenuousFunction;
            this.initialFunction = initialFunction;

            InitializeDGWithGaussLobattoNodes();
        }

        private void InitializeDGWithGaussLobattoNodes()
        {
            fluxEvaluation = new Vector(N + 1);
            
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out nodes, out integrationWeights);

            solution = ComputeStartSolution() ;
            Console.Write(solution.toString(15));
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
           
        public void UpdateSolution(Vector solution)
        {
            this.solution = solution;
        }     
        public void UpdateBorderValues()
        {
            this.LeftBorderValue = solution[0];
            this.RightBorderValue = solution[N];
        }

        public Vector EvaluateTimeDerivative(double time)
        {
            Vector timeDerivative;      
            numFluxLeft = numericalFluxFunction(LeftNeighbour.RightBorderValue, solution[0]);
            numFluxRight = numericalFluxFunction(solution[nodes.Length - 1], RightNeighbour.LeftBorderValue);

            Vector volume = new Vector(N+1);
            Vector surface = new Vector(N+1);

            surface[0] = numFluxLeft;
            surface[N] = numFluxRight;

            for (int i = 0; i <= N ; i++)
                volume[i] = fluxFunction(solution[i]);

            //timeDerivative = invMassMatrix * ( (((!differentialMatrix) * massMatrix) * volume) - (B * surface));

            //timeDerivative = invMassMatrix * ((volumeMatrix* volume) - (B * surface));
            timeDerivative = volumeMatrix * volume - surfaceMatrix * surface;

            timeDerivative = (Vector)((2.0 / (rightSpaceBoundary - leftSpaceBoundary)) * timeDerivative);
            timeDerivative += ComputeInhomogenuousPart(time);
            return timeDerivative;
        }

        private Vector ComputeInhomogenuousPart(double time)
        {
            Vector inhomoPart = new Vector(N + 1);
            for (int i = 0; i < N + 1; i++)
                inhomoPart[i] = inhomogenuousFunction(nodes[i], time);

            return inhomoPart;
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

        private double MapToOriginSpace(double x)
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
