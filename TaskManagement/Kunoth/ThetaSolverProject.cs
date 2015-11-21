using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.PDE;
using Structures;

namespace TaskManagement.Kunoth
{
    public class ThetaSolverProject
    {
        public void evaluate()
        {
            //ThetaSolver ts = new ThetaSolver();
            //double leftBoudary = -Math.PI;
            //double rightBoundary = -leftBoudary;
            //double timeBoundary = 1.0;
            //int M = 32;
            //int N = 32;
            //double theta = 0.50;

            //Vector res = ts.SolvePDE(initalFunction, leftBoundaryFunction, rightBoundaryFunction, leftBoudary, rightBoundary, timeBoundary, M, N, theta);
            //Vector exact = evaluateExactSolution(timeBoundary, leftBoudary, rightBoundary, N);
            ////Console.WriteLine(res.toString(15));

            //double absErr =  NSharp.Measures.MeasureFunctions.CalculateDiscreteNorm((res - exact), (res - exact));
            //double relErr = absErr/NSharp.Measures.MeasureFunctions.CalculateDiscreteNorm(exact, exact);
            //Console.WriteLine(relErr);

            double[] relError = computeRelativeError(5);
            double[] EOC = computeEOC(relError);
            Console.Write("Rel Error:");
            for(int i = 0; i < relError.Length; i++)
            {
                Console.WriteLine(relError[i]);
            }
            Console.Write("EOC:");
            for (int i = 0; i < EOC.Length; i++)
            {
                Console.WriteLine(EOC[i]);
            }
            Console.ReadKey();
        }

        public double[] computeRelativeError(int maxPow)
        {
                double[] EOC = new double[maxPow - 1];
                double[] relError = new double[maxPow];
                ThetaSolver ts = new ThetaSolver();
                double leftBoudary = -Math.PI;
                double rightBoundary = -leftBoudary;
                double timeBoundary = 1.0;
                int M;
                int N;
                double theta = 0.50;
                Vector exact;
                Vector res;

                for (int j = 1; j <= maxPow; j++) {
                    M = (int)Math.Pow(2.0, j);
                    N = M;
                    //Berechne exakte und numerische Lösung
                    res = ts.SolvePDE(initalFunction, leftBoundaryFunction, rightBoundaryFunction, leftBoudary, rightBoundary, timeBoundary, M, N, theta);
                    exact = evaluateExactSolution(timeBoundary, leftBoudary, rightBoundary, N);

                    //Fehler berechnen
                    double absErr = NSharp.Measures.MeasureFunctions.CalculateDiscreteNorm((res - exact), (res - exact));
                    relError[j - 1] = absErr / NSharp.Measures.MeasureFunctions.CalculateDiscreteNorm(exact, exact);
                }

            return relError;
        }

        private double[] computeEOC(double[] relError)
        {
            double[] EOC = new double[relError.Length - 1];
            for (int j = 0; j < EOC.Length - 1; j++)
            {
                EOC[j] = Math.Log(relError[j] / relError[j + 1]) / Math.Log(2.0);
            }
            return EOC;
        }

        private static Vector evaluateExactSolution(double time, double leftSpaceBoundary, double rightSpaceBoundary, int N)
        {
            Vector evaluation = new Vector(N - 1);
            double spaceStep = (Math.Abs(leftSpaceBoundary) + Math.Abs(leftSpaceBoundary)) / (double)N;

            for(int i = 0; i< N-1; i++)
            {
                double spaceLocation = leftSpaceBoundary + (i + 1) * spaceStep;
                evaluation[i] = exactSolution(time, spaceLocation);
            }

            return evaluation;
        }

        private static double exactSolution(double time, double space)
        {
            return Math.Exp(-time) * Math.Sin(space);
        }

        private static double initalFunction(double x)
        {
            return Math.Sin(x);
        }

        private static double leftBoundaryFunction(double t)
        {
            return 0.0;
        }

        private static double rightBoundaryFunction(double t)
        {
            return 0.0;
        }
    }
}
