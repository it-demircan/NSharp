using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public class IntegrationToolbox
    {

        public static Matrix generateMassMatrix(Vector weights)
        {
            return new Matrix(weights);
        }

        public static double computeGaussianIntegrationWithGaussNodesAndWeights(Func<double,double> myFunction, int N)
        {
            Vector nodes, weights;
            LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(N, out nodes, out weights);
            double result = computeIntegralSummation(myFunction, nodes, weights);
            return result;
        }

        public static double computeGaussianIntegrationWithGaussLobattoNodesAndWeights(Func<double, double> myFunction, int N)
        {
            Vector nodes, weights;
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out nodes, out weights);
            double result = computeIntegralSummation(myFunction, nodes, weights);
            return result;
        }

        private static double computeIntegralSummation(Func<double, double> myFunction, Vector nodes, Vector weights)
        {
            double evaluation = 0.0;
            for (int i = 0; i <nodes.Length; i++)
            {
                evaluation += myFunction(nodes[i]) * weights[i];
            }
            return evaluation;
        }
    }
}
