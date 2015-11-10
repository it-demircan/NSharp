using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public class InterpolationToolbox
    {

        public static double evaluateLagrangePolynome(double x, Vector nodes, int j)
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

            Vector barycentricWeights = computeBarycentricWeights(nodes);
            double tempProd = 1.0;

            for(int i = 0; i < nodes.Length; i++)
            {
                if (i != j)
                    tempProd *= (x - nodes[i]);
            }

            tempProd *= barycentricWeights[j];

            return tempProd;
        }
        public static double evaluateLagrangeRepresentation(double x, Vector nodes, Vector functionValues)
        {
            int idx;
            //Gibt die Position zurück, wenn x einer Stützerstelle entspricht.
            if( (idx = nodes.ContainsValue(x)) != -1)
                return functionValues[idx];

            Vector barycentricWeights = computeBarycentricWeights(nodes);
            double tempQuot = 0.0;
            double denominator = 0.0;
            double numerator = 0.0;

            for(int i = 0; i < nodes.Length; i++)
            {
                tempQuot = barycentricWeights[i] / (x - nodes[i]);
                numerator += functionValues[i] * tempQuot;
                denominator += tempQuot;
            }

            return numerator / denominator;
        }

        public static Vector computeBarycentricWeights(Vector nodes)
        {
            double tempProd = 0.0;
            int N = nodes.Length;
            Vector baryWeights = new Vector(N);
            
            for(int i = 0; i < N; i++)
            {
                tempProd = 1.0;
                for(int j = 0; j < N; j++)
                { 
                    if(i != j)
                        tempProd *= (nodes[i] - nodes[j]);
                }
                baryWeights[i] = 1.0/tempProd;
            }

            return baryWeights;
        }

        public static Matrix computeLagrangePolynomeDerivativeMatrix(Vector nodes)
        {   
            Matrix D = new Matrix(nodes.Length, nodes.Length);
            Vector baryWeights = computeBarycentricWeights(nodes);

            for (int i = 0; i < nodes.Length; i++)
            {
                D[i, i] = 0.0;
                for(int j = 0; j < nodes.Length; j++)
                {
                    if(i != j)
                    {
                        D[i, j] = baryWeights[j] / baryWeights[i] * (1.0/(nodes[i] - nodes[j])) ;
                        D[i, i] -= D[i, j];
                    }
                }
            }

            return D;
        }

        public static void testLagrangeEvaluation()
        {
            Vector nodes = new Vector(2);
            Vector functionValues = new Vector(2);

            nodes[0] = 0.0;
            nodes[1] = 1.0;

            functionValues[0] = 0.0;
            functionValues[1] = 1.0;
            double x = evaluateLagrangeRepresentation(Math.PI, nodes, functionValues);
            bool res = GeneralHelper.isXAlmostEqualToY(x, Math.PI);
            Console.Write(res);
            Console.ReadKey();
        }    
        public static void testBarycentricWeights()
        {
            Vector test = new Vector(3);
            test[0] = 1.0;
            test[1] = 2.0;
            test[2] = 3.0;
            Vector res = computeBarycentricWeights(test);

            Console.Write(res.toString(4));
            Console.ReadKey();

        }
    }
}
