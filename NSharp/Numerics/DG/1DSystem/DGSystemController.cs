using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._1DSystem
{
    public enum CalculationMode
    {
        EOC, Energy, WellBalanced, Stone
    }

    public enum Task
    {
        TaskOne, TaskTwo, TaskThree
    }

    public class DGSystemController
    {
        int polynomOrder;//N
        int systemDimension;
        double spaceLengthInElements;
        CalculationMode calcMode;
        Task myTask;

        const double CFL = 1.45;
       
        DGSystemElement[] elements;


        public void createDGElements(int numberOfDGElements, int polynomOrder, double leftBoundary, double rightBoundary, int systemDimension)
        {
            calcMode = CalculationMode.Stone;
            myTask = Task.TaskTwo;

            this.polynomOrder = polynomOrder;
            this.systemDimension = systemDimension;
            spaceLengthInElements = (rightBoundary - leftBoundary) / (double)numberOfDGElements;
            elements = new DGSystemElement[numberOfDGElements];
            for (int i = 0; i < numberOfDGElements; i++)
            {
                double leftSpaceBorder = leftBoundary + (double)i * (rightBoundary - leftBoundary) / (double)numberOfDGElements;
                double rightSpaceBorder = leftBoundary + (double)(i + 1) * (rightBoundary - leftBoundary) / (double)numberOfDGElements;

                if (myTask ==Task.TaskOne)
                    elements[i] = new DGSystemElement(DGMODE.STRONG, calcMode, leftSpaceBorder, rightSpaceBorder, polynomOrder,systemDimension, NumFlux, InhomogenuousPart, FluxFunction, InitialFunction);
                else if(myTask == Task.TaskTwo)
                    elements[i] = new DGSystemElement(DGMODE.ECON, calcMode,leftSpaceBorder, rightSpaceBorder, polynomOrder, systemDimension, NumFluxECON, InhomogenuousPart, FluxFunction, InitialFunction);
                else if(myTask == Task.TaskThree)
                    elements[i] = new DGSystemElement(DGMODE.ECON, calcMode, leftSpaceBorder, rightSpaceBorder, polynomOrder, systemDimension, NumFluxES, InhomogenuousPart, FluxFunction, InitialFunction);

                if (i > 0)
                {
                    elements[i].LeftNeighbour = elements[i - 1];
                    elements[i - 1].RightNeighbour = elements[i];
                }
            }

           
            if (calcMode != CalculationMode.Stone)
            {   //Periodic Boundary Condition
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
            else
            {
                //Dirichlet 
                DGSystemElement dummy = new DGSystemElement(DGMODE.ECON, calcMode, 0.0, 0.0, polynomOrder, systemDimension, NumFluxES, InhomogenuousPart, FluxFunction, InitialFunction);
                Vector leftDirichletBoundary = new Vector(2);
                leftDirichletBoundary[0] = 1.1;
                leftDirichletBoundary[1] = leftDirichletBoundary[0] * -0.525;
                
                Vector rightDirichletBoundary = new Vector(2);
                rightDirichletBoundary[0] = 1.1;
                rightDirichletBoundary[1] = rightDirichletBoundary[0] * 0.525;

                dummy.RightBoarderValue = leftDirichletBoundary;
                dummy.LeftBoarderValue = rightDirichletBoundary;

                elements[0].LeftNeighbour = dummy;
                elements[numberOfDGElements - 1].RightNeighbour = dummy;
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
        public Dictionary<double,double> ComputeSolution(double endTime, double timeStep = 0.0)
        {
            double recentTime = 0.0;
            Matrix[] recentTimeDerivatives = new Matrix[elements.Length];
            double recentTimeStep = timeStep;

            Dictionary<double,double> energys = new Dictionary<double,double>();

            while (recentTime < endTime)
            {
                if (timeStep == 0.0)
                {
                    double lambdaMax = GetMaximumLambdaOverall();
                    recentTimeStep = ComputeTimeStep(CFL, lambdaMax);
                }
                //Console.WriteLine("RecentTime: " + recentTime);

                if (recentTime + recentTimeStep > endTime)
                    recentTimeStep = endTime - recentTime;
                //Vektoren zurüksetzen
                for (int i = 0; i < elements.Length; i++)
                    recentTimeDerivatives[i] = new Matrix(polynomOrder + 1, systemDimension);
                Console.WriteLine("TIME: " + recentTime + "\t\t ENERGY: " + this.ComputeEnergy());
                energys.Add(recentTime, this.ComputeEnergy());
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

            return energys;
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


        public Vector GetOriginNodes()
        {
            Vector constant = new Vector((polynomOrder + 1) * elements.Length);
            for (int i = 0; i < elements.Length; i++)
            {
                constant.InjectMatrixAtPosition(elements[i].GetOriginNodes(), i * (polynomOrder + 1), 0);
            }
            return constant;
        }

        /**
          Problem Abhängige Funktionen
        **/

        public static double GRAVITATION = 1.0;

        public double ComputeEnergy()
        {
            double energy = 0.0;
            int counter = 1;
            foreach (DGSystemElement element in elements)
            {
                double inte= element.ComputeEnergy(GRAVITATION, Ground);
                //Console.WriteLine("Integral CellNo:"+ counter +" Int:" + inte);
                energy +=  element.ComputeEnergy(GRAVITATION, Ground);
                counter++;
            }
            return energy;
        }

        public Vector ComputeMass()
        {
            Vector massVector = new Vector(systemDimension);
            foreach(DGSystemElement element in elements)
            {
                massVector = massVector + element.ComputeMass();
            }
            return massVector;
        }

        
        public Vector ComputeConstant()
        {
            Vector constant = new Vector((polynomOrder + 1) * elements.Length);
            for(int i = 0; i < elements.Length; i++)
            {
                constant.InjectMatrixAtPosition(elements[i].GetConstantSolution(Ground), i * (polynomOrder+1), 0);
            }
            return constant;
        }

        public Vector InitialFunction(Vector nodes, bool fromLeft = true)
        {
            Vector result = new Vector(systemDimension);

            //Aufgabe 1D)
            if(calcMode == CalculationMode.WellBalanced)
            {
                result[0] = 3.0 - Ground(nodes[0]);
                result[1] = 0.0;
            }
            else if(calcMode == CalculationMode.Energy)
            {
                result[0] = nodes[0] <= 10.0 ? 3.0 - Ground(nodes[0]) : 2.5 - Ground(nodes[0]);
                result[1] = 0.0;
            }
            else if(calcMode ==CalculationMode.EOC)
            {   //EOC
                result[0] = Math.Sin(nodes[0]) + 2.0;//Math.Sin(2*Math.PI*nodes[0])+2.0;
                result[1] = result[0];
            }
            else if(calcMode == CalculationMode.Stone)
            {
                result[0] = 1.1 - Ground(nodes[0]);
                result[1] = (nodes[0] <= 5.0 && fromLeft ? -0.525 : 0.525)* result[0];
            }
            return result;
        }

        public double Ground(double node)
        {
            //A1A
            if(calcMode == CalculationMode.EOC)
                return 0.0;
            else if (calcMode == CalculationMode.Energy || calcMode == CalculationMode.Stone)
            {
                return (Math.Abs(node - 10.0) <= 2 ? Math.Sin((Math.PI / 4.0) * node) : 0.0);
            }
            else//WellBalanced
            {
                return Math.Sin((Math.PI / 4.0) * node);
            }   
        }

        public Vector InhomogenuousPart(Vector solution,double node, double time)
        {
            Vector inhomo = new Vector(systemDimension);
            if (calcMode == CalculationMode.WellBalanced)
            {
                //Aufgabe 1D
                inhomo[0] = 0;
                inhomo[1] = -1.0 * GRAVITATION * solution[0] * (Math.PI / 4.0) * Math.Cos((Math.PI / 4.0) * node);
            }else if(calcMode == CalculationMode.Energy || calcMode == CalculationMode.Stone)
            {   //ENERGY
                inhomo[0] = 0.0;
                inhomo[1] = -1.0 * GRAVITATION * solution[0];//* (Math.Abs(node - 10.0) <= 2.0 ? (Math.PI / 4.0) * Math.Cos((Math.PI / 4.0) * node) : 0.0);
            }
            else
            {   //EOC
                inhomo[0] = 0;
                inhomo[1] = -Math.Cos(node - time) + Math.Cos(time - node) * (GRAVITATION * (-Math.Sin(time - node)) + 2.0 * GRAVITATION + 1.0);

            }

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

        public Vector FluxFunction(Vector solution)
        {
            Vector result = new Vector(systemDimension);

            result[0] = solution[1];
            result[1] = ((solution[1] / solution[0]) * (solution[1] / solution[0]) * solution[0]) + (1.0 / 2.0) * GRAVITATION * solution[0] * solution[0];
            return result;
        }

        public Vector NumFlux(Vector left, Vector right, double node)
        {
            Vector temp = FluxFunction(left) + FluxFunction(right);
            return (Vector)((1.0 / 2.0) * temp) - (Vector)((ComputeMaxEigenWert(left, right) / 2.0) * (right - left));
        }

        public Vector NumFluxECON(Vector left, Vector right, double node)
        {
            Vector numFlux = new Vector(2);
            double vLeft = left[1] / left[0];
            double vRight = right[1] / right[0];

            numFlux[0] = ((vLeft + vRight) / 2.0) * ((left[0] + right[0]) / 2.0);
            numFlux[1] = ((vLeft + vRight) / 2.0) * ((vLeft + vRight) / 2.0) * ((left[0] + right[0]) / 2.0)
                + (1.0 / 2.0) * GRAVITATION * ((left[0] * left[0] + right[0] * right[0]) / 2.0);
            return numFlux;

        }

        public Vector NumFluxES(Vector left, Vector right, double node)
        {
            Vector result = new Vector(2);
            Vector econ = NumFluxECON(left, right, node);

            //Mitterlwert berechnen fuer A und R
            //Matrix R = (1.0 / 2.0) * (ComputeRMatrix(right) + ComputeRMatrix(left));
            //Vector eigenValues = (Vector)((1.0/2.0)* (ComputeEigenwert(right) + ComputeEigenwert(left)));
            //Matrix A = new Matrix(2, 2);
            //A[0, 0] = eigenValues[1];
            //A[1, 1] = eigenValues[0];

            //Vector test = (Vector)((1.0 / 2.0) * (right + left));
            //Matrix R = ComputeRMatrix(test);

            //Matrix A = new Matrix(2, 2);
            //A[0, 0] = Math.Abs(eigenValues[1]);
            //A[1, 1] = Math.Abs(eigenValues[0]);

            Vector eigenValues = ComputeEigenwert((Vector)((1.0 / 2.0) * (right + left)));
            double negEW = eigenValues[1];
            double posEW = eigenValues[0];

            Vector deriEnergy = ComputeEnergyDerivative(right, node) - ComputeEnergyDerivative(left, node);
            Vector esVector = new Vector(2);
            esVector[0] = (1.0 / (4.0 * GRAVITATION)) * (((Math.Abs(negEW) + Math.Abs(posEW)) * deriEnergy[0]) + (Math.Abs(negEW) * negEW + Math.Abs(posEW) * posEW) * deriEnergy[1]);
            esVector[1] = (1.0 / (4.0 * GRAVITATION)) * ((Math.Abs(negEW) * negEW + Math.Abs(posEW) * posEW) * deriEnergy[0] + (Math.Abs(negEW) * negEW * negEW + Math.Abs(posEW) * posEW * posEW) * deriEnergy[1]);

            result = econ - esVector;
            return result;
        }

        private Matrix ComputeRMatrix(Vector vector)
        {
            double h = vector[0];
            double hv = vector[1];
            Vector eigenVectorNegativ = new Vector(2);
            Vector eigenVectorPositiv = new Vector(2);
            eigenVectorNegativ[0] = -(h / ((Math.Sqrt(GRAVITATION * h * h * h)) - hv));
            eigenVectorNegativ[1] = 1.0;

            eigenVectorPositiv[0] = (h / ((Math.Sqrt(GRAVITATION * h * h * h)) + hv));
            eigenVectorPositiv[1] = 1.0;

            Matrix R = new Matrix(2, 2);
            R.InjectMatrixAtPosition(eigenVectorNegativ, 0, 0);
            R.InjectMatrixAtPosition(eigenVectorPositiv, 0, 1);

            return R;
        }

        private Vector ComputeEnergyDerivative(Vector vector, double node)
        {
            Vector result = new Vector(2);
            result[0] = -(1.0 / 2.0) * ((vector[1] * vector[1]) / (vector[0] * vector[0])) + GRAVITATION * vector[0] + GRAVITATION * Ground(node);
            result[1] = vector[1] / vector[0];
            return result;
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
            Eigenvalues[0] = hv/h + Math.Sqrt(GRAVITATION * h);
            Eigenvalues[1] = hv/h - Math.Sqrt(GRAVITATION * h);
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

            List<Vector> EigenValues = new List<Vector>();
            for (int i = 0; i < AllEvaluations.Count; i++)
                EigenValues.Add(ComputeEigenwert(AllEvaluations.ElementAt(i)));

            maxEigenvalue = ComputeMaxAbsolutEigenvalue(EigenValues);
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
