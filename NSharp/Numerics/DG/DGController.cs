using NSharp.Numerics.Interpolation;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public class DGController
    {

        //Variablen für Runga Kutta
        double[] A = { 0.0,
                -567301805773.0/1357537059087.0,
                -2404267990393.0/ 2016746695238.0,
                - 3550918686646.0/ 2091501179385.0,
                -1275806237668.0/ 842570457699.0};

        double[] B = { 0.0,
            1432997174477.0 / 9575080441755.0,
            2526269341429.0/ 6820363962896.0,
            2006345519317.0/ 3224310063776.0,
            2802321613138.0/ 2924317926251.0 };

        double[] C = { 1432997174477.0/ 9575080441755.0,
            5161836677717.0 / 13612068292357.0,
            1720146321549.0 / 2090206949498.0,
            3134564353537.0 / 4481467310338.0,
            2277821191437.0 / 14882151754819.0};

        DGElement[] elements;
        int polynomOrder; //N
        double veloCity = 2.0;

        IntegrationMode myMode;
        double spaceLengthInElements;

        public void createDGElements(int numberOfDGElements,IntegrationMode mode, int polynomOrder, double leftBoundary, double rightBoundary)
        {
            this.polynomOrder = polynomOrder;
            myMode= mode;
            spaceLengthInElements = (rightBoundary-leftBoundary)/(double)numberOfDGElements;
            elements = new DGElement[numberOfDGElements];
            for (int i = 0; i<numberOfDGElements; i++)
            {
                double leftSpaceBorder = leftBoundary + (double)i * (rightBoundary - leftBoundary) / (double)numberOfDGElements;
                double rightSpaceBorder = leftBoundary + (double)(i + 1) * (rightBoundary-leftBoundary) / (double)numberOfDGElements;
                elements[i] = new DGElement(mode,leftSpaceBorder, rightSpaceBorder, polynomOrder, FluxFunction, NumFlux, InhomogenuousPart, InitialFunction);

                if (i > 0)
                {
                    elements[i].LeftNeighbour = elements[i - 1];
                    elements[i - 1].RightNeighbour = elements[i];
                }
            }

            //Periodic Boundary Condition
            if (numberOfDGElements > 1)
            {
                elements[0].LeftNeighbour = elements[numberOfDGElements - 1];
                elements[numberOfDGElements - 1].RightNeighbour = elements[0];
            }
            else
            {
                elements[0].LeftNeighbour = elements[0];
                elements[0].RightNeighbour = elements[0];
            }
        }


        public double computeSolution(double endTime, double timeStep)
        {
            double recentTime = 0.0;
            Vector[] tempTimeDerivaties = new Vector[elements.Length];

            while (recentTime < endTime)
            {
                if (recentTime + timeStep > endTime)
                    timeStep = endTime - recentTime;
                //Vektoren zurüksetzen
                for (int i = 0; i < elements.Length; i++)
                    tempTimeDerivaties[i] = new Vector(polynomOrder+1);

               //Runga Kutta
                for (int k = 0; k < 5; k++)
                {
                    //Runga Kutta für jedes Element, damit Zeitsynchron berechnet wird
                    for (int i = 0; i < elements.Length; i++)
                    {                        
                        Vector solution = elements[i].GetSolution();
                        double nextTimeStep = recentTime + B[k] * timeStep;
                        tempTimeDerivaties[i] = (Vector)(A[k] * tempTimeDerivaties[i]) + elements[i].EvaluateTimeDerivative(nextTimeStep);
                        solution = solution + (Vector)(C[k] * timeStep * tempTimeDerivaties[i]);
                        elements[i].UpdateSolution(solution);
                    }
                    //Rand updaten
                    for (int i = 0; i < elements.Length; i++)
                    {
                        elements[i].UpdateBorderValues();
                    }
                }

                recentTime += timeStep;
            }

            //for (int k = 0; k < elements.Length; k++)
            //{
            //    for (int i = 0; i < elements[k].GetSolution().Length; i++)
            //        Console.WriteLine(elements[k].GetSolution()[i]);
            //}
            return computeError(endTime);
        }


        public double ComputeTimeStep(double cfl){
            return myMode == IntegrationMode.GaussLobatto ? cfl * ComputeSpaceLengthInElement()  / veloCity : cfl * ComputeSpaceLengthInElement() / veloCity;
        }

        public Matrix ConstructDGLMatrix()
        {
            return myMode == IntegrationMode.GaussLobatto ? ConstructDGLMatrixWithGaussLobattoNodes() : ConstructDGLMatrixWithGaussNodes();
        }

        public Matrix ConstructDGLMatrixWithGaussLobattoNodes()
        {
            Matrix A = new Matrix((polynomOrder+1) * elements.Length, (polynomOrder+1)  * elements.Length);
            Matrix S = elements[0].GetSMatrix();
            Matrix D = elements[0].GetDifferentialMatrix();
            Matrix E_0N = new Matrix((polynomOrder + 1), (polynomOrder + 1));
            E_0N[0, polynomOrder] = 1.0;
            Matrix E_NN = new Matrix((polynomOrder + 1), (polynomOrder + 1));
            E_0N[polynomOrder, polynomOrder] = 1.0;

            Matrix MainDiagonal = -1.0*S + D + S * E_NN;
            Matrix LowerDiagonal = S * E_0N;

            for (int i = 0; i < elements.Length-1; i++)
            {
                A.InjectMatrixAtPosition(MainDiagonal, i * (polynomOrder + 1), i * (polynomOrder + 1));
                A.InjectMatrixAtPosition(LowerDiagonal, (i + 1) * (polynomOrder + 1), i * (polynomOrder + 1));
            }
            A.InjectMatrixAtPosition(MainDiagonal, (elements.Length - 1) * (polynomOrder + 1), (elements.Length - 1) * (polynomOrder + 1));
            A.InjectMatrixAtPosition(LowerDiagonal, 0, (elements.Length - 1) * (polynomOrder + 1));

            A = -1.0 * veloCity / (elements[0].GetSpaceLength() / 2.0) * A;
            return A;
        }

        public Matrix ConstructDGLMatrixWithGaussNodes()
        {
            Matrix A = new Matrix((polynomOrder + 1) * elements.Length, (polynomOrder + 1) * elements.Length);
            Matrix invMass = elements[0].GetInverseMassMatrix();
            Matrix massMatrix = elements[0].GetMassMatrix();
            Matrix D = elements[0].GetDifferentialMatrix();

            Matrix L2 = elements[0].GetL2Matrix();
            Matrix L1 = elements[0].GetL1Matrix();
            Matrix B = elements[0].GetBMatrix();

            Matrix MainDiagonal = -1.0*massMatrix * (!D) * massMatrix + invMass * B * L2;
            Matrix LowerDiagonal = invMass * B * L1;

            for (int i = 0; i < elements.Length - 1; i++)
            {
                A.InjectMatrixAtPosition(MainDiagonal, i * (polynomOrder + 1), i * (polynomOrder + 1));
                A.InjectMatrixAtPosition(LowerDiagonal, (i + 1) * (polynomOrder + 1), i * (polynomOrder + 1));
            }
            A.InjectMatrixAtPosition(MainDiagonal, (elements.Length - 1) * (polynomOrder + 1), (elements.Length - 1) * (polynomOrder + 1));
            A.InjectMatrixAtPosition(LowerDiagonal, 0, (elements.Length - 1) * (polynomOrder + 1));

            A = -1.0 * veloCity / (elements[0].GetSpaceLength() / 2.0) * A;

            return A;
        }


        public Vector getCompleteSolution()
        {
            Vector solution = new Vector(elements.Length * (polynomOrder + 1));
            for (int i = 0; i < elements.Length; i++)
            {
                for (int k = 0; k <= polynomOrder; k++)
                    solution[i * (polynomOrder + 1) + k] = elements[i].GetSolution()[k];
            }

            return solution;
        }
        
        public Vector getOriginSpace()
        {
            Vector space = new Vector(elements.Length * (polynomOrder + 1));
            for (int i = 0; i < elements.Length; i++)
            {
                for (int k = 0; k <= polynomOrder; k++)
                    space[i * (polynomOrder + 1) + k] = elements[i].GetOriginNodes()[k];
            }

            return space;
        }
        public double ComputeSpaceLengthInElement()
        {
            return myMode == IntegrationMode.GaussLobatto ? spaceLengthInElements / (double)polynomOrder : spaceLengthInElements / ((double)polynomOrder+2);
        }
        //Berechnet den L2 Fehler indem die Lösung mit 2N Polynomen approximiert wird.
        public double computeError(double endTime)
        {
            Vector errorNodes, errorWeights;
            LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(2*polynomOrder, out errorNodes, out errorWeights);
            LagrangeInterpolator interpolator = new LagrangeInterpolator(elements[0].nodes);

            double error = 0.0;
            double tempErr = 0.0;

            for (int i = 0; i < elements.Length; i++)
            {
                for (int m = 0; m < errorNodes.Length; m++)
                {
                    tempErr = interpolator.evaluateLagrangeRepresentation(errorNodes[m], elements[i].GetSolution());
                    double trafoSpace = elements[i].MapToOriginSpace(errorNodes[m]);
                    tempErr = (ExactSolution(trafoSpace, endTime) - tempErr) * errorWeights[m]* (ExactSolution(trafoSpace, endTime) - tempErr) * ((elements[i].rightSpaceBoundary - elements[i].leftSpaceBoundary) / 2.0);
                    error += tempErr;
                }         
            }
            return Math.Sqrt(error);
        }

        public Vector computeExactSolution(Vector nodes, double time)
        {
            Vector exactSolution = new Vector(nodes.Length);
            for (int i = 0; i < nodes.Length; i++)
                exactSolution[i] = ExactSolution(nodes[i], time);
            return exactSolution;
        }

        private double InitialFunction(double space)
        {
            return Math.Sin(veloCity * Math.PI * space);
            //if (space >= 0.3 && space <= 0.7)
            //    return 1.0;
            //else
            //    return 0.0;
        }

        private double InhomogenuousPart(double space, double time)
        {
            return 0.0;
        }

        private double FluxFunction(double u)
        {
            return veloCity * u;
        }

        private double NumFlux(double left, double right)
        {
            return FluxFunction(left);
            //return 1.0 / 2.0 * (FluxFunction(left) + FluxFunction(right)) - 1.0 / 12.0 *((right - left) * (right - left));
        }


        private double ExactSolution(double space, double time) 
        {
            return Math.Sin(2.0 * Math.PI * (space - veloCity * time));
        }

    }
}
