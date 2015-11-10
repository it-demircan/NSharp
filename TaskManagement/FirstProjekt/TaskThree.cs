using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;
using NSharp.Numerics.DG;

namespace TaskManagement.FirstProjekt
{
    //Aufgabe 3 - SBP Eigenschaft numerisch zeigen für GL-Stützstellen und widerlegen, dass diese nicht bei Legendre Stützstellen funktioniert.
    class TaskThree
    {
        static int[] elements = {1,2,3,4};
        public void evaluate()
        {
            Matrix gaußLobatto;
            Matrix gaußLegendre;
            foreach (int i in elements)
            {
                gaußLobatto = computeSBPMatrixByGaußLobatto(i);
                gaußLegendre = computeSBPMatrixByGaußLegendre(i);
            }
        }


        private Matrix computeSBPMatrixByGaußLobatto(int N)
        {
            Vector nodes, weights;
            //SBP-Eigenschaft bei Gauß-Lobatto
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out nodes, out weights);

            Matrix B = computeSBPMatrix(nodes, weights);
            return B;
        }

        private Matrix computeSBPMatrixByGaußLegendre(int N)
        {
            Vector nodes, weights;
            //SBP-Eigenschaft bei Gauß-Lobatto
            LegendrePolynomEvaluator.computeLegendreGaussNodesAndWeights(N, out nodes, out weights);

            Matrix B = computeSBPMatrix(nodes, weights);
            return B;
        }

        private Matrix computeSBPMatrix(Vector nodes, Vector weights)
        {
            Matrix massMatrix = IntegrationToolbox.generateMassMatrix(weights);
            Matrix diffMatrix = InterpolationToolbox.computeLagrangePolynomeDerivativeMatrix(nodes);

            Matrix left = massMatrix * diffMatrix;
            Matrix B = left + !left;
            return B;
        }
    }
}
