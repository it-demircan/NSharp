using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._1DSystem
{
    public class DGSystemController
    {
        int polynomOrder;//N
        int systemDimension;
        double spaceLengthInElements;
       
        DGSystemElement[] elements;


        public void createDGElements(int numberOfDGElements, int polynomOrder, double leftBoundary, double rightBoundary, int systemDimension)
        {
            this.polynomOrder = polynomOrder;
            this.systemDimension = systemDimension;
            spaceLengthInElements = (rightBoundary - leftBoundary) / (double)numberOfDGElements;
            elements = new DGSystemElement[numberOfDGElements];
            for (int i = 0; i < numberOfDGElements; i++)
            {
                double leftSpaceBorder = leftBoundary + (double)i * (rightBoundary - leftBoundary) / (double)numberOfDGElements;
                double rightSpaceBorder = leftBoundary + (double)(i + 1) * (rightBoundary - leftBoundary) / (double)numberOfDGElements;
                elements[i] = new DGSystemElement(leftSpaceBorder, rightSpaceBorder, polynomOrder,systemDimension, NumFlux, InhomogenuousPart, FluxFunction, InitialFunction);

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


        public Matrix GetSolution()
        {
            Matrix SolutionSystem = new Matrix(elements.Length * (polynomOrder + 1), systemDimension);

            for(int i = 0; i < elements.Length; i++)
            {
                Matrix cellSolution = elements[i].GetSolution();
                SolutionSystem.InjectMatrixAtPosition(cellSolution, i * (polynomOrder + 1), 0);
            }
            return SolutionSystem;
        }

        //timestep == 0.0 => automatisch
        public void ComputeSolution(double endTime, double timeStep = 0.0)
        {
            double recentTime = 0.0;
            Matrix[] recentTimeDerivatives = new Matrix[elements.Length];
            double recentTimeStep = timeStep;


            while (recentTime < endTime)
            {
                if (timeStep == 0.0)
                {
                    double lambdaMax = GetMaximumLambdaOverall();
                    recentTimeStep = ComputeTimeStep(0.5, lambdaMax);
                }
                //Console.WriteLine("RecentTime: " + recentTime);

                if (recentTime + recentTimeStep > endTime)
                    recentTimeStep = endTime - recentTime;
                //Vektoren zurüksetzen
                for (int i = 0; i < elements.Length; i++)
                    recentTimeDerivatives[i] = new Matrix(polynomOrder + 1, systemDimension);

                //Runge Kutte
                for (int k = 0; k < 5; k++)
                {
                    for (int i = 0; i < elements.Length; i++)
                    {
                        Matrix solutionSystem = elements[i].GetSolution();
                        double nextTimeStep = recentTime + B[k] * recentTimeStep;
                        Matrix EvaluatedTimeDerivative = elements[i].EvaluateTimeDerivativeGaussLobatto(nextTimeStep);

                        for (int sysIdx = 0; sysIdx < systemDimension; sysIdx++)
                        {
                            Vector tempTimeDerivative = (Vector)(A[k] * recentTimeDerivatives[i].GetColumn(sysIdx)) + EvaluatedTimeDerivative.GetColumn(sysIdx);
                            recentTimeDerivatives[i].InjectMatrixAtPosition(tempTimeDerivative, 0, sysIdx);
                            Vector solution = solutionSystem.GetColumn(sysIdx) + (Vector)(C[k] * recentTimeStep * recentTimeDerivatives[i].GetColumn(sysIdx));
                            solutionSystem.InjectMatrixAtPosition(solution, 0, sysIdx);
                        }

                        elements[i].UpdateSolutionSystem(solutionSystem);
                    }
                    //Rand updaten
                    for (int i = 0; i < elements.Length; i++)
                    {
                        elements[i].UpdateBorderValues();
                    }
                }

                recentTime += recentTimeStep;
            }
        }

        public double ComputeTimeStep(double cfl)
        {
            //TODO
            return cfl * ComputeSpaceLengthInElement() / 2.0;
        }

        public double ComputeSpaceLengthInElement()
        {
            return spaceLengthInElements / (double)(polynomOrder + 1);
        }

        public double ComputeTimeStep(double cfl, double lambdaMax)
        {
            return cfl * ComputeSpaceLengthInElement() / lambdaMax;
        }


        /**
          Problem Abhängige Funktionen
        **/
        const double GRAVITATION = 1.0;

        public Vector InitialFunction(Vector nodes)
        {
            Vector result = new Vector(systemDimension);
            //Test
            result[0] = Math.Sin(nodes[0]) + 2.0;//Math.Sin(2*Math.PI*nodes[0])+2.0;
            result[1] = result[0];

            //Aufgabe 1 - A)
            //result[0] = Math.Sin(2.0*Math.PI*nodes[0]) + 2.0;//Math.Sin(2*Math.PI*nodes[0])+2.0;
            //result[1] = result[0];

            return result;
        }

        public Vector FluxFunction(Vector solution)
        {
            Vector result = new Vector(systemDimension);

            result[0] = solution[1];
            result[1] = ( (solution[1]/solution[0]) * (solution[1] / solution[0]) * solution[0]) + (1.0 / 2.0) * GRAVITATION * solution[0] * solution[0];
            return result;
        }

        public Vector NumFlux(Vector left,Vector right)
        {
            Vector temp = FluxFunction(left) + FluxFunction(right);
            return (Vector)((1.0 / 2.0) * temp) - (Vector)((ComputeMaxEigenWert(left, right) / 2.0) * (right - left));
        }

        public Vector InhomogenuousPart(Vector solution,double node, double time)
        {
            Vector inhomo = new Vector(systemDimension);
            //Test
            inhomo[0] = 0;
            inhomo[1] = -Math.Cos(node - time) + Math.Cos(time - node) * (GRAVITATION * (-Math.Sin(time - node)) + 2.0 * GRAVITATION + 1.0);

            //Aufgabe 1 - A)
            //inhomo[0] = 0;
            //inhomo[1] = -GRAVITATION * solution[0] * 2.0 * Math.PI * Math.Cos(2.0 * Math.PI * (node - time));

            ////Residuum
            //inhomo[1] += 8.0 * GRAVITATION * Math.PI* Math.Cos(2.0 * Math.PI * (node - time)) + 4.0 * GRAVITATION *Math.PI* Math.Sin(2.0 * Math.PI * (node - time)) * Math.Cos(2.0 * Math.PI * (node - time));
            //End Aufgabe 1 - A)

            return inhomo;
        }
        
        public Matrix ComputeExactSolution(double time)
        {
            Matrix ExactSolution = new Matrix((polynomOrder + 1) * elements.Length, systemDimension);

            for(int i = 0; i < elements.Length; i++)
            {
                Vector nodes;
                for(int k = 0; k < (nodes=elements[i].GetOriginNodes()).Length; k++)
                {
                    //TEST 
                    ExactSolution[i * (nodes.Length) + k, 0] = 2.0 + Math.Sin((nodes[k] - time));
                    ExactSolution[i * (nodes.Length) + k, 1] = 1.0;

                    //Aufgabe 1 - A
                    //ExactSolution[i * (nodes.Length) + k, 0] = 2.0 + Math.Sin(2.0 * Math.PI * (nodes[k] - time));
                    //ExactSolution[i * (nodes.Length) + k, 1] = 1.0;
                    //End Aufgabe 1 - A
                }
            }
            return ExactSolution;
        }
        private double ComputeMaxEigenWert(Vector left, Vector right)
        {
            List<Vector> boundaryValues = new List<Vector>() { ComputeEigenwert(left), ComputeEigenwert(right) };
            double maxEigenvalue = ComputeMaxAbsolutEigenvalue(boundaryValues);
            return maxEigenvalue;
        }

        private Vector ComputeEigenwert(Vector evaluated)
        {
            double h = evaluated[0];
            double hv = evaluated[1];
            Vector Eigenvalues = new Vector(2);
            Eigenvalues[0] = hv / h + Math.Sqrt(GRAVITATION * h);
            Eigenvalues[1] = hv / h - Math.Sqrt(GRAVITATION * h);
            return Eigenvalues;
        }
        public double ComputeMaxAbsolutEigenvalue(List<Vector> vectorsWithEigenvalues)
        {
            double maxEigenValue = 0.0;

            for (int i = 0; i < vectorsWithEigenvalues.Count; i++)
            {
                Vector tmp = vectorsWithEigenvalues.ElementAt(i);
                for(int m = 0; m < tmp.Length; m++) {
                    if(Math.Abs(tmp[m])> maxEigenValue)
                    {
                        maxEigenValue = Math.Abs(tmp[m]);
                    }
                }
            }

            return maxEigenValue;
        }

        public double GetMaximumLambdaOverall()
        {
            List<Vector> AllEvaluations = new List<Vector>();
            double maxEigenvalue;

            for (int i = 0; i < elements.Length; i++)
            {
                DGSystemElement temp = elements[i];
                for(int m = 0; m < temp.GetSolution().NoRows;m++)
                {
                    AllEvaluations.Add(temp.GetSolution().GetRowAsColumnVector(m));
                }
            }

            maxEigenvalue = ComputeMaxAbsolutEigenvalue(AllEvaluations);
            return maxEigenvalue;
        }

        #region Runge Kutte 
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

    #endregion 
    }
}
