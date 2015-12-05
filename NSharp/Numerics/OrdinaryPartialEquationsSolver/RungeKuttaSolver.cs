using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.OrdinaryPartialEquationsSolver
{
    public class RungeKuttaSolver : IODESolver
    {
        static readonly double[] A = { 0.0,
            -567301805773.0/1357537059087.0,
         -2404267990393.0/ 2016746695238.0,
         - 3550918686646.0/ 2091501179385.0,
         -1275806237668.0/ 842570457699.0};

        static readonly double[] B = { 0.0,
        1432997174477.0 / 9575080441755.0,
        2526269341429.0/ 6820363962896.0,
        2006345519317.0/ 3224310063776.0,
        2802321613138.0/ 2924317926251.0 };

        static readonly double[] C = { 1432997174477.0/ 9575080441755.0,
        5161836677717.0 / 13612068292357.0,
        1720146321549.0 / 2090206949498.0,
        3134564353537.0 / 4481467310338.0,
        2277821191437.0 / 14882151754819.0};

        public Vector computeSolutionForNextStep(Vector initial,OrdinaryDifferentialEquation ode, double startTime, double endTime)
        {
            Vector tempSolution = initial;
            Vector evaluatedODE = ode.evaluateOrdinaryDifferentialEquation(tempSolution, startTime);
            double deltaTime = endTime - startTime;

            for(int i = 0; i < 5; i++)
            {
                double nextTimeStep = startTime + B[i] * deltaTime;
                evaluatedODE = (Vector) (A[i] * evaluatedODE) + ode.evaluateOrdinaryDifferentialEquation(tempSolution, nextTimeStep);
                tempSolution = tempSolution + (Vector)(C[i] * deltaTime * evaluatedODE);
            }

            return tempSolution;
        }


        public Vector computeSolutionWithMultipleSteps(Vector initial, OrdinaryDifferentialEquation ode, double startTime, double endTime, double step)
        {
            double tempTime = startTime;
            int N =  Convert.ToInt32((endTime - startTime)/ step);
            Vector tempSolution = initial;

            for(int i = 1; i <= N; i++)
            {
                tempSolution = computeSolutionForNextStep(tempSolution, ode, tempTime, tempTime + step);
                tempTime += step;
            }

            if (!GeneralHelper.isXAlmostEqualToY(startTime + N * step, endTime))
            {
                tempTime = startTime + N * step;
                tempSolution = computeSolutionForNextStep(tempSolution, ode, tempTime, endTime);
            }

            return tempSolution;
        }


        public Vector computeSolutionVectorWithMultipleSteps(Vector initial, OrdinaryDifferentialEquation ode, double startTime, double endTime, double step)
        {

            double tempTime = startTime;
            int N = Convert.ToInt32((endTime - startTime) / step) + 1;
            Vector tempSolution = initial;

            Vector solution = new Vector(N);
            solution[0] = initial[0];

            Vector temp = new Vector(1);

            for (int i = 1; i < N; i++)
            {

                temp[0] = solution[i - 1];
                solution[i] = computeSolutionForNextStep(temp, ode, tempTime, tempTime + step)[0];
                tempTime += step;
            }

            if (!GeneralHelper.isXAlmostEqualToY(startTime + (N-1) * step, endTime))
            {
                tempTime = startTime + N * step;
                tempSolution = computeSolutionForNextStep(tempSolution, ode, tempTime, endTime);
            }

            return solution;
        }
    }
}
