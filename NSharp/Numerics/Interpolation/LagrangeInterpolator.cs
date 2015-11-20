using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.Interpolation
{
    class LagrangeInterpolator : IInterpolater
    {
        Vector nodes;
        Vector barycentricWeights;

        public LagrangeInterpolator(Vector nodes)
        {
            this.nodes = nodes;
            barycentricWeights = computeBarycentricWeights(this.nodes);
        }
        public Vector computeBarycentricWeights(Vector nodes)
        {
            double tempProd = 0.0;
            int N = nodes.Length;
            Vector baryWeights = new Vector(N);

            for (int i = 0; i < N; i++)
            {
                tempProd = 1.0;
                for (int j = 0; j < N; j++)
                {
                    if (i != j)
                        tempProd *= (nodes[i] - nodes[j]);
                }
                baryWeights[i] = 1.0 / tempProd;
            }

            return baryWeights;
        }

        public Matrix computeLagrangePolynomeDerivativeMatrix()
        {
            Matrix D = new Matrix(nodes.Length, nodes.Length);
            Vector baryWeights = barycentricWeights;

            for (int i = 0; i < nodes.Length; i++)
            {
                D[i, i] = 0.0;
                for (int j = 0; j < nodes.Length; j++)
                {
                    if (i != j)
                    {
                        D[i, j] = baryWeights[j] / baryWeights[i] * (1.0 / (nodes[i] - nodes[j]));
                        D[i, i] -= D[i, j];
                    }
                }
            }

            return D;
        }

        public double evaluateLagrangePolynome(double x, int j)
        {
            int idx;
            //Gibt die Position zurück, wenn x einer Stützerstelle entspricht.
            if ((idx = nodes.ContainsValue(x)) != -1)
            {
                if (idx == j)
                    return 1.0;
                else
                    return 0.0;
            }

            double tempProd = 1.0;

            for (int i = 0; i < nodes.Length; i++)
            {
                if (i != j)
                    tempProd *= (x - nodes[i]);
            }

            tempProd *= barycentricWeights[j];

            return tempProd;
        }

        public double evaluateLagrangeRepresentation(double x, Vector functionValues)
        {
            int idx;
            //Gibt die Position zurück, wenn x einer Stützerstelle entspricht.
            if ((idx = nodes.ContainsValue(x)) != -1)
                return functionValues[idx];
           
            double tempQuot = 0.0;
            double denominator = 0.0;
            double numerator = 0.0;

            for (int i = 0; i < nodes.Length; i++)
            {
                tempQuot = barycentricWeights[i] / (x - nodes[i]);
                numerator += functionValues[i] * tempQuot;
                denominator += tempQuot;
            }

            return numerator / denominator;
        }

        public double evaluateInterpolation(double x, Vector nodes, Vector functionValues)
        {
            this.nodes = nodes;
            barycentricWeights = computeBarycentricWeights(nodes);
            return evaluateLagrangeRepresentation(x, functionValues);
        }
    }
}
