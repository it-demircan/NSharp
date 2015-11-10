using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Modelling.LocationOptimization
{
    /// <summary>
    /// This Class implements an Multi Facility Solver according to 
    /// 
    /// [1]    C. Charalambous, An iterative Algortihm for the Multi Facility Minimax Location Problem with euclidean distances
    ///        Naval Research Logistics Quarterly, Volume 28, Issue 2, pages 325-337, 1981.
    /// </summary>
    public class MultiFacilityProblemSolver
    {
        /// <summary>
        /// Example of Use
        /// </summary>
        /// <returns>Optimized Location</returns>
        public static Matrix StartTest()
        {
            const int MAX_ITER = 100;
            const double MAX_ERR = 0.000000001;

            //Position of existing Facilites 
            /*
             * (0,1)
             * 
             * (0,0)    (1,0) 
             * */
            Matrix oldFacilities = new Matrix(3, 2);
            oldFacilities[0, 0] = 0.0;
            oldFacilities[0, 1] = 0.0;
            oldFacilities[1, 0] = 1.0;
            oldFacilities[1, 1] = 0.0;
            oldFacilities[2, 0] = 0.0;
            oldFacilities[2, 1] = 1.0;

            //Weights between new and existing Facilities
            Matrix quantifierOldLocationNewLocation = new Matrix(2, 3);
            quantifierOldLocationNewLocation[0, 0] = 1;
            quantifierOldLocationNewLocation[0, 1] = 1;
            quantifierOldLocationNewLocation[0, 2] = 0;
            quantifierOldLocationNewLocation[1, 0] = 0;
            quantifierOldLocationNewLocation[1, 1] = 5;
            quantifierOldLocationNewLocation[1, 2] = 2;

            //Upper Triangular Matrix without diagonal elements
            Matrix quantifierNewLocationNewLocation = new Matrix(2, 2);
            quantifierNewLocationNewLocation[0, 1] = 0;

            Matrix result = SolveMFP(oldFacilities, quantifierOldLocationNewLocation, quantifierNewLocationNewLocation, MAX_ITER, MAX_ERR);
            return result;
        }



        /// <summary>
        /// Solve the Multi Facility Problem according to [1]. 
        /// </summary>
        /// <param name="oldFacilityPositions">Position of existing facilities</param>
        /// <param name="quantifierOldLocationNewLocation">Weights between new and existing facilities</param>
        /// <param name="quantifierNewLocationNewLocation">Weights between new facilties themselves</param>
        /// <param name="maxIteration">Number of Iterations</param>
        /// <param name="maxError">Stop criterion</param>
        /// <returns></returns>
        public static Matrix SolveMFP(Matrix oldFacilityPositions, Matrix quantifierOldLocationNewLocation, Matrix quantifierNewLocationNewLocation, int maxIteration = 100, double maxError = 0.001)
        {
            Matrix newFacilityPosition = null;
            int n = quantifierOldLocationNewLocation.NoRows; //# New Facilities //W
            int m = quantifierOldLocationNewLocation.NoColumns;//# Old Facilites //V

            Matrix lambda = new Matrix(n, m, MatrixType.ONES);

            //For stop criterion
            double minMaxDistance = 0.0, maxCost = 0.0;
            double relError = 0.0;

            //Construct µ (upper triangular matrix with ones but without diagonal elements, which are set to 0)
            Matrix mu = new Matrix(n, n);
            for (int l = 0; l < n - 1; l++)
            {
                for (int k = l + 1; k < n; k++)
                {
                    mu[l, k] = 1.0;
                }
            }

            for (int i = 0; i < maxIteration; i++)
            {
                newFacilityPosition = SolveQuadraticFunction(oldFacilityPositions, quantifierOldLocationNewLocation, quantifierNewLocationNewLocation, lambda, mu);
                double S = ComputeS(newFacilityPosition, oldFacilityPositions, lambda, mu, quantifierOldLocationNewLocation, quantifierNewLocationNewLocation);
                lambda = UpdateLambda(newFacilityPosition, oldFacilityPositions, lambda, quantifierOldLocationNewLocation, S);
                mu = UpdateMu(newFacilityPosition, mu, quantifierNewLocationNewLocation, S);

                minMaxDistance = ComputeS(newFacilityPosition, oldFacilityPositions, lambda, mu, quantifierOldLocationNewLocation, quantifierNewLocationNewLocation);
                maxCost = ComputeMaximumCost(oldFacilityPositions, newFacilityPosition, quantifierOldLocationNewLocation, quantifierNewLocationNewLocation);

                //Stop Criterion
                if ((relError = Math.Abs((minMaxDistance - maxCost) / maxCost)) < maxError)
                {
                    i = maxIteration;
                }
                
            }
            return newFacilityPosition;
        }



        /// <summary>
        /// Compute the maximum cost based on (2) in [1]
        /// </summary>
        /// <returns>Maximum Cost</returns>
        private static double ComputeMaximumCost(Matrix oldFacilityPositions, Matrix newFacilityPositions, Matrix quantifierOldLocationNewLocation, Matrix quantifierNewLocationNewLocation)
        {
            double maxCost = 0.0;
            int n = newFacilityPositions.NoRows;
            int m = oldFacilityPositions.NoRows;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    double tempCost = quantifierOldLocationNewLocation[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityPositions[i, 0], newFacilityPositions[i, 1], oldFacilityPositions[j, 0], oldFacilityPositions[j, 1]);
                    if (tempCost >= maxCost)
                    {
                        maxCost = tempCost;
                    }
                }
            }

            for (int i = 0; i < n - 1; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    double tempCost = quantifierNewLocationNewLocation[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityPositions[i, 0], newFacilityPositions[i, 1], newFacilityPositions[j, 0], newFacilityPositions[j, 1]);
                    if (tempCost >= maxCost)
                    {
                        maxCost = tempCost;
                    }
                }
            }

            return maxCost;
        }

        /// <summary>
        /// Compute S according to [1] at chapter 3 step 3
        /// </summary>
        /// <returns>Returns S</returns>
        private static double ComputeS(Matrix newFacilityPosition, Matrix oldFacilityPosition, Matrix lambda, Matrix mu, Matrix w, Matrix v)
        {
            double s = 0.0;
            int n = newFacilityPosition.NoRows;
            int m = oldFacilityPosition.NoRows;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    s += lambda[i, j] * w[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityPosition[i, 0], newFacilityPosition[i, 1], oldFacilityPosition[j, 0], oldFacilityPosition[j, 1]);
                }
            }

            for (int i = 0; i < n - 1; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    s += mu[i, j] * v[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityPosition[i, 0], newFacilityPosition[i, 1], newFacilityPosition[j, 0], newFacilityPosition[j, 1]);
                }
            }

            return s;
        }

        /// <summary>
        /// Update Lambda^{r-1} to Lambda^{r} according to [1] at chapter 3 step 3
        /// </summary>
        /// <returns>Lambda^{r}</returns>
        private static Matrix UpdateLambda(Matrix newFacilityPosition, Matrix oldFacilityPosition, Matrix lambda, Matrix w, double S)
        {
            int n = newFacilityPosition.NoRows;
            int m = oldFacilityPosition.NoRows;
            Matrix updatedLambda = new Matrix(n, m);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    updatedLambda[i, j] = (lambda[i, j] * w[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityPosition[i, 0], newFacilityPosition[i, 1], oldFacilityPosition[j, 0], oldFacilityPosition[j, 1])) / S;
                }
            }

            return updatedLambda;
        }

        /// <summary>
        /// Update Mu^{r-1} to Mu^{r} according to [1] at chapter 3 step 3
        /// </summary>
        /// <returns>Mu^{r}</returns>
        private static Matrix UpdateMu(Matrix newFacilityLocation, Matrix mu, Matrix v, double S)
        {
            int n = newFacilityLocation.NoRows;
            Matrix updatedMu = new Matrix(n, n);

            for (int i = 0; i < n - 1; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    updatedMu[i, j] = mu[i, j] * v[i, j] * Measures.MeasureFunctions.CalculateEuclidianDistance(newFacilityLocation[i, 0], newFacilityLocation[i, 1], newFacilityLocation[j, 0], newFacilityLocation[j, 1]);
                }
            }

            return updatedMu;
        }

        #region Optimum Solution of the Quadratic Function

        /// <summary>
        /// Computes the optimum solution for x,y according to [1] chapter 3.1
        /// </summary>
        /// <param name="oldFacilityPositions">Position of existing locations </param>
        /// <param name="quantifierOldLocationNewLocation">Weight Matrix W</param>
        /// <param name="quantifierNewLocationNewLocation">Weigth Matrix V</param>
        /// <param name="lambda">Minimax Multifiers Lambda</param>
        /// <param name="mu">Minimax Multifiers Mu</param>
        /// <returns>Optimum Solution of X,Y as an Matrix (X|Y)</returns>
        private static Matrix SolveQuadraticFunction(Matrix oldFacilityPositions, Matrix quantifierOldLocationNewLocation, Matrix quantifierNewLocationNewLocation, Matrix lambda, Matrix mu)
        {
            int n = quantifierOldLocationNewLocation.NoRows;
            Matrix scaledW = ComputeScaledWeight(quantifierOldLocationNewLocation, lambda);
            Matrix scaledV = ComputeScaledWeight(quantifierNewLocationNewLocation, mu);
            //See [1] at (10a,b)
            Vector aTilde = scaledW * oldFacilityPositions.GetColumn(0); // X Position of locations
            Vector bTilde = scaledW * oldFacilityPositions.GetColumn(1); // Y Position of locations

            //Create Matrix A, see [1] at (11)
            Matrix A = new Matrix(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                        A[i, j] = ComputeBeta(scaledW, scaledV, i);
                    else if (i < j)
                    {
                        A[i, j] = (-1) * scaledV[i, j];
                    }
                    else
                        A[i, j] = (-1) * scaledV[j, i];
                }
            }
            //Solve Ax = aTilde, see [1] 13a
            Vector X = LinearAlgebra.GaußEliminationSolver.Solve((Matrix)A.Clone(), aTilde);
            //Solve Ax = bTilde, see [1] 13b
            Vector Y = LinearAlgebra.GaußEliminationSolver.Solve((Matrix)A.Clone(), bTilde);

            Matrix result = new Matrix(X.NoRows, 2);

            for (int i = 0; i < X.NoRows; i++)
            {
                result[i, 0] = X[i];
                result[i, 1] = Y[i];
            }

            return result;
        }

        /// <summary>
        /// Calculate scaled minimax multipliers \tilde{\lambda}} and \tilde{\mu}
        /// </summary>
        /// <param name="miniMaxMultiplier">Lambda oder mu - Minimax Multiplier Matrix </param>
        /// <param name="weights">Weight Matrix w or v</param>
        /// <returns></returns>
        private static Matrix ComputeScaledWeight(Matrix miniMaxMultiplier, Matrix weights)
        {
            Matrix result = new Matrix(weights.NoRows, weights.NoColumns);
            int i, k;
            for (i = 0; i < weights.NoRows; i++)
            {
                for (k = 0; k < weights.NoColumns; k++)
                {
                    result[i, k] = miniMaxMultiplier[i, k] * weights[i, k] * weights[i, k];
                }
            }
            return result;
        }

        /// <summary>
        /// Computes Beta for Matrix A, see [1] at (12) for more Information
        /// </summary>
        /// <param name="scaledWeightW">Scaled Weight W Matrix</param>
        /// <param name="scaledWeightV">Scaled Weight V Matrix</param>
        /// <param name="index">Index for Beta</param>
        /// <returns>Returns Beta_{index}</returns>
        private static double ComputeBeta(Matrix scaledWeightW, Matrix scaledWeightV, int index)
        {
            double result = 0.0;
            int n = scaledWeightW.NoRows;
            int m = scaledWeightW.NoColumns;


            for (int j = 0; j < m; j++)
            {
                result += scaledWeightW[index, j];
            }

            for (int j = 0; j < n; j++)
            {
                if (j != index)
                    result += scaledWeightV[index, j];
            }

            return result;
        }
        #endregion

    }
}
