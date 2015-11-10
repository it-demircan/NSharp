using NSharp.Numerics.DG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;
using NSharp;


namespace TestProject
{
    class Program
    {
        
        static void Main(string[] args)
        {
            //DateTime StartZeit = DateTime.Now;
            ////Hier die Funktion einfügen deren Zeit gemessen werden soll


            //double result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(firstFunction, N);
            //DateTime EndZeit = DateTime.Now;
            //Console.WriteLine(result);
            //TimeSpan GemessendeZeit = EndZeit - StartZeit;
            //Console.Write(GemessendeZeit);

            //Vector nodes, weights;

            ////SBP-Eigenschaft bei Gauß-Lobatto
            //LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N,out nodes,out weights);
            //Matrix massMatrix = IntegrationToolbox.generateMassMatrix(weights);
            //Matrix diffMatrix = InterpolationToolbox.computeLagrangePolynomeDerivativeMatrix(nodes);

            //Matrix left = massMatrix * diffMatrix;
            //Matrix B = left + !left;

            //Console.WriteLine(B.toString(4));

            ////SBP-Eigenschaft funktioniert bei Legendre Stützstellen nicht
            //LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(N, out nodes, out weights);
            //massMatrix = IntegrationToolbox.generateMassMatrix(weights);
            //diffMatrix = InterpolationToolbox.computeLagrangePolynomeDerivativeMatrix(nodes);

            //left = massMatrix * diffMatrix;
            //B = left + !left;

            //Console.WriteLine(B.toString(4));

            //Console.ReadKey();
        }


    }
}
