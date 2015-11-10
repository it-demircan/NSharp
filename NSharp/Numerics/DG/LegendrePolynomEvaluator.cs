using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public class LegendrePolynomEvaluator
    {
        public const int NEWTON_ITERATOR = 4;

        /// <summary>
        /// Berechne die Auswertung der Legendre Polynome und deren Ableitung vom Grad N an der Stelle x.
        /// </summary>
        /// <param name="N">Grad N </param>
        /// <param name="x">Auswertungsstelle</param>
        /// <param name="evaPolynomial">Ausgabe der Evaluation der Legendre Polynome</param>
        /// <param name="evaDerivative">Ausgabe der Evaluation der Ableitung der Legendre Polynome</param>
        public static void evaluateLegendrePolynomialAndDerivative(int N, double x, out double evaPolynomial, out double evaDerivative)
        {
            if(N == 0)
            {
                evaPolynomial = 1.0;
                evaDerivative = 0.0;
            }
            else if(N == 1)
            {
                evaPolynomial = x;
                evaDerivative = 1.0;
            }
            else
            {
                evaPolynomial = 0.0; evaDerivative = 0.0;
                double tempEvaPolynom_N_2 = 1.0, tempEvaPolynom_N_1 = x;
                double tempEvaDerivative_N_2 = 0.0, tempEvaDerivative_N_1 = 1.0;

                for(int k = 2; k <= N; k++)
                {
                    evaPolynomial = ((2.0 * k - 1.0) / (double)k) * x * tempEvaPolynom_N_1 - (k - 1.0) / (double)k * tempEvaPolynom_N_2;
                    evaDerivative = tempEvaDerivative_N_2 + (2.0 * k - 1.0) * tempEvaPolynom_N_1;

                    //Update rekursiv Variablen
                    tempEvaPolynom_N_2 = tempEvaPolynom_N_1;
                    tempEvaPolynom_N_1 = evaPolynomial;
                    tempEvaDerivative_N_2 = tempEvaDerivative_N_1;
                    tempEvaDerivative_N_1 = evaDerivative;
                }
            }
        }

        /// <summary>
        /// Berechnet die N+1 Gauß Legendre Punkte und die zugehörigen Gewichte.
        /// </summary>
        /// <param name="N">Grad vom Legendre Polynom</param>
        /// <param name="nodes">Ausgabe der Stützstellen</param>
        /// <param name="weights">Ausgabe der Gewichte</param>
        public static void computeLegendreGaussNodesAndWeights(int N, out Vector nodes, out Vector weights)
        {
            Vector gaussNodes = new Vector(N + 1);
            Vector gaussWeights = new Vector(N + 1);

            if (N == 0)
            {
                gaussNodes[0] = 0.0;
                gaussWeights[0] = 2.0;
            }
            else if (N == 1)
            {
                gaussNodes[0] = -Math.Sqrt(1.0 / 3.0);
                gaussNodes[1] = -gaussNodes[0];

                gaussWeights[0] = 1.0;
                gaussWeights[1] = gaussWeights[0];
            }
            else
            {
                int upperBound = (N + 1) / 2 - 1;
                double evaLegendre, evaLegendreDerivative, delta;

                for (int j = 0; j <= upperBound; j++)
                {
                    gaussNodes[j] = -Math.Cos((2.0 * j + 1) / (2.0 * N + 2.0) * Math.PI);
                    for (int k = 0; k <= NEWTON_ITERATOR; k++)
                    {
                        evaluateLegendrePolynomialAndDerivative(N + 1, gaussNodes[j], out evaLegendre, out evaLegendreDerivative);
                        delta = -evaLegendre / evaLegendreDerivative;
                        gaussNodes[j] += delta;
                        if (Math.Abs(delta) <= 2 * GeneralHelper.EPSILON * Math.Abs(gaussNodes[j]))
                            k = NEWTON_ITERATOR + 1;
                    }

                    evaluateLegendrePolynomialAndDerivative(N + 1, gaussNodes[j], out evaLegendre, out evaLegendreDerivative);

                    gaussNodes[N - j] = -gaussNodes[j];
                    gaussWeights[j] = 2.0 / ((1 - gaussNodes[j] * gaussNodes[j]) * evaLegendreDerivative * evaLegendreDerivative);
                    gaussWeights[N - j] = gaussWeights[j];
                }

                if (N % 2 == 0)
                {
                    evaluateLegendrePolynomialAndDerivative(N + 1, 0.0, out evaLegendre, out evaLegendreDerivative);
                    gaussNodes[N / 2] = 0.0;
                    gaussWeights[N / 2] = 2.0 / (evaLegendreDerivative * evaLegendreDerivative);
                }
            }

            nodes = gaussNodes;
            weights = gaussWeights;
        }

        /// <summary>
        /// Berechne die Auswertung der Legendrepolynome an der Stelle x und inneren Stützpolynomen q und dessen Ableitung.
        /// </summary>
        /// <param name="N">Grad N der Legendre Polynome</param>
        /// <param name="x">Auswertungsstelle</param>
        /// <param name="q">Innere Stützpolynome</param>
        /// <param name="qDerivative">Ableitung der inneren Stützpolynome</param>
        /// <param name="evaPolynomial">Auswertung der Legendre Polynome</param>
        public static void evaluateQandL(int N, double x, out double q, out double qDerivative, out double evaPolynomial)
        {
            int k;
            double evaDerivative;
            double tempEvaPolynom_N_2 = 1.0, tempEvaPolynom_N_1 = x, tempEvaPolynom_N_Plus_1;
            double tempEvaDerivative_N_2 = 0.0, tempEvaDerivative_N_1 = 1.0, tempEvaDerivative_N_Plus_1;

            for(k = 2; k < N; k++)
            {
                evaPolynomial = ((2.0 * k - 1.0) / (double)k) * x * tempEvaPolynom_N_1 - (k - 1.0) / (double)k * tempEvaPolynom_N_2;
                evaDerivative = tempEvaDerivative_N_2 + (2.0 * k - 1.0) * tempEvaPolynom_N_1;

                tempEvaPolynom_N_2 = tempEvaPolynom_N_1;
                tempEvaPolynom_N_1 = evaPolynomial;
                tempEvaDerivative_N_2 = tempEvaDerivative_N_1;
                tempEvaDerivative_N_1 = evaDerivative;
            }

            //k == N - Lauf
            evaPolynomial = ((2.0 * k - 1.0) / (double)k) * x * tempEvaPolynom_N_1 - (k - 1.0) / (double)k * tempEvaPolynom_N_2;
            evaDerivative = tempEvaDerivative_N_2 + (2.0 * k - 1.0) * tempEvaPolynom_N_1;

            k = N + 1;
            tempEvaPolynom_N_Plus_1 = ((2.0 * k - 1.0) / (double)k) * x * evaPolynomial - ((k - 1.0) / (double)k) * tempEvaPolynom_N_1;
            tempEvaDerivative_N_Plus_1 = tempEvaDerivative_N_1 + (2.0 * k - 1.0) * evaPolynomial;

            q = tempEvaPolynom_N_Plus_1 - tempEvaPolynom_N_1;
            qDerivative = tempEvaDerivative_N_Plus_1 - tempEvaDerivative_N_1;
        }

        /// <summary>
        /// Berechnet die Gauß Lobatto Stützstellen und deren zugehörige Gewichte.
        /// </summary>
        /// <param name="N">Grad N des Legendre Stützpolynoms</param>
        /// <param name="nodes">Ausgabe der GL-Stützstellen</param>
        /// <param name="weights">Ausgabe der GL-Gewichte</param>
        public static void computeGaussLobattoNodesAndWeights(int N, out Vector nodes, out Vector weights)
        {
            if( N == 0)
            {
                throw new Exception("Gauß Lobatto Nodes and Weights not defined for N = 0.");
            }
            Vector gaussNodes = new Vector(N + 1);
            Vector gaussWeights = new Vector(N + 1);
            if (N == 1)
            {
                gaussNodes[0] = -1.0;
                gaussWeights[0] = 1.0;
                gaussNodes[1] = 1.0;
                gaussWeights[1] = 1.0;
            }
            else
            {
                int upperBound = (N + 1) / 2 - 1;
                double q, qDerivative, evaPolynomial, delta;

                gaussNodes[0] = -1;
                gaussWeights[0] = 2.0 / (N * (N + 1));
                gaussNodes[N] = 1;
                gaussWeights[N] = gaussWeights[0];
                
                for(int j=1;j<= upperBound; j++)
                {
                    gaussNodes[j] = -Math.Cos( ((j + 0.25) * Math.PI / (double)N) - ((3.0 / (8.0 * N * Math.PI)) * (1.0 / (j + 0.25))) );

                    for(int k = 0; k <= NEWTON_ITERATOR; k++)
                    {
                        evaluateQandL(N, gaussNodes[j], out q,out qDerivative,out evaPolynomial);
                        delta = -q / qDerivative;
                        gaussNodes[j] += delta;
                        if (Math.Abs(delta) <= 2 * GeneralHelper.EPSILON * Math.Abs(gaussNodes[j]))
                            k = NEWTON_ITERATOR + 1;
                    }
                    evaluateQandL(N, gaussNodes[j], out q, out qDerivative, out evaPolynomial);
                    gaussNodes[N - j] = -gaussNodes[j];
                    gaussWeights[j] = 2.0 / (N * (N + 1) * evaPolynomial * evaPolynomial);
                    gaussWeights[N - j] = gaussWeights[j];
                }

                if(N%2 == 0)
                {
                    evaluateQandL(N, 0.0, out q, out qDerivative, out evaPolynomial);
                    gaussNodes[N / 2] = 0.0;
                    gaussWeights[N / 2] = 2.0 / (N * (N + 1) * evaPolynomial * evaPolynomial);
                }
            }

            nodes = gaussNodes;
            weights = gaussWeights;
        }
    }
}
