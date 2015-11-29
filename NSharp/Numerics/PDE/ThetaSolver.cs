using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.PDE
{
    public class ThetaSolver
    {
        Func<double, double> initialFunction;
        Func<double, double> leftBoundaryFunction;
        Func<double, double> rightBoundaryFunction;

        //Space 
        double leftBoundary;
        double rightBoundary;

        const double startTime = 0.0;
        double timeBoundary;

        double theta;
        double alpha;

        int M; //Anzahl der Zerlegungen von Zeitraum
        int N; //Anzahl der Zerlegungen von Raum


        public Vector SolvePDE(Func<double, double> initialFunction,
                                Func<double, double> leftBoundaryFunction,
                                Func<double, double> rightBoundaryFunction, double leftBoundary,
                                double rightBoundary,double timeBoundary, int M , int N, double theta)
        {
            //Initialisiere Werte
            Vector result;
            this.M = M;
            this.N = N;
            this.theta = theta;
            this.leftBoundary = leftBoundary;
            this.rightBoundary = rightBoundary;
            this.timeBoundary = timeBoundary;
            this.initialFunction = initialFunction;
            this.leftBoundaryFunction = leftBoundaryFunction;
            this.rightBoundaryFunction = rightBoundaryFunction;

           
            double spaceStepLength = computeSpaceStepLength();
            double timeStepLength = computeTimeStepLength();
            this.alpha = computeAlpha(timeStepLength, spaceStepLength);

            Matrix leftMatrix = composeLeftMatrix();
            Matrix rightMatrix = composeRightMatrix();
            Vector initalSolutionVector = evaluateInitialFunction(spaceStepLength);

            
            //Hier könnte man die LU Zerlegung verwendung, da die Linke Matrix sich nicht verändert.
            //Matrix[] decompsedMatrix = NSharp.LinearAlgebra.Decomposer.Decompose(leftMatrix, NSharp.LinearAlgebra.DecomposerType.LU);

            //#1 Berechnung mit Anfangsbedingung
            Vector tempVector = rightMatrix * initalSolutionVector;
            tempVector = tempVector + (computeRightSideBoundaryVector(startTime) - computeLeftSideBoundaryVector(startTime + timeStepLength));

            result = NSharp.LinearAlgebra.GaußEliminationSolver.Solve((Matrix)leftMatrix.Clone(), tempVector);

            //Restliche Berechnungen
            for (int k = 1; k < M; k++)
            {
                tempVector = rightMatrix * result;
                tempVector = tempVector + (computeRightSideBoundaryVector(startTime+timeStepLength*k) - computeLeftSideBoundaryVector(startTime + timeStepLength*(k+1.0)));
                result = NSharp.LinearAlgebra.GaußEliminationSolver.Solve((Matrix)leftMatrix.Clone(), tempVector);
            }
            return result;
        }
        private double computeSpaceStepLength()
        {
            return (rightBoundary + Math.Abs(leftBoundary)) / (double)N;
        }

        private double computeTimeStepLength()
        {
            return timeBoundary / (double)M;
        }

        private double computeAlpha(double timeStepLength, double spaceStepLength)
        {
            return timeStepLength / (spaceStepLength*spaceStepLength);
        }

        private Matrix composeLeftMatrix()
        {
            Matrix leftMatrix = new Matrix(N - 1, N - 1);

            if(N == 2)
            {
                leftMatrix[0, 0] = 2.0 * alpha * theta + 1;
                return leftMatrix;
            }

            for(int i = 1; i< leftMatrix.NoRows-1; i++)
            {
                leftMatrix[i,i] = 2.0 * alpha * theta + 1;
                leftMatrix[i, i-1] = -alpha * theta;
                leftMatrix[i, i+1] = leftMatrix[i, i-1];
            }

            leftMatrix[0, 0] = leftMatrix[1,1];
            leftMatrix[0, 1] = leftMatrix[1,2];

            leftMatrix[leftMatrix.NoRows - 1, leftMatrix.NoRows - 1] = leftMatrix[0, 0];
            leftMatrix[leftMatrix.NoRows - 1, leftMatrix.NoRows - 2] = leftMatrix[0, 1];

            return leftMatrix;
        }
        
        private Matrix composeRightMatrix()
        {
            Matrix rightMatrix = new Matrix(N - 1, N - 1);

            if (N == 2)
            {
                rightMatrix[0, 0] = -(2.0 * alpha * (1.0 - theta) - 1);
                return rightMatrix;
            }

            for (int i = 1; i < rightMatrix.NoRows - 1; i++)
            {
                rightMatrix[i, i] = -(2.0*alpha*(1.0-theta)-1);
                rightMatrix[i, i - 1] = alpha * (1.0-theta);
                rightMatrix[i, i + 1] = rightMatrix[i, i - 1];
            }

            rightMatrix[0, 0] = rightMatrix[1,1];
            rightMatrix[0, 1] = rightMatrix[1,2];

            rightMatrix[rightMatrix.NoRows - 1, rightMatrix.NoRows - 1] = rightMatrix[0, 0];
            rightMatrix[rightMatrix.NoRows - 1, rightMatrix.NoRows - 2] = rightMatrix[0, 1];

            return rightMatrix;
        }

        private Vector computeLeftSideBoundaryVector(double time)
        {
            Vector leftBoundaryVector = new Vector(N - 1);
            leftBoundaryVector[0] = -alpha * theta * leftBoundaryFunction(time);
            leftBoundaryVector[N - 2] = -alpha * theta * rightBoundaryFunction(time);
            return leftBoundaryVector;
        }

        private Vector computeRightSideBoundaryVector(double time)
        {
            Vector rightSideBoundaryVector = new Vector(N - 1);
            rightSideBoundaryVector[0] = alpha*(1.0- theta) * leftBoundaryFunction(time);
            rightSideBoundaryVector[N - 2] = alpha * (1.0 - theta) * rightBoundaryFunction(time);
            return rightSideBoundaryVector;
        }
       
        private Vector evaluateInitialFunction(double spaceStepLength)
        {
            Vector evaluation = new Vector(N - 1);

            for (int i = 0; i < evaluation.Length; i++)
                evaluation[i] = initialFunction(leftBoundary + (i + 1) * spaceStepLength);

            return evaluation;
        }
    }
}
