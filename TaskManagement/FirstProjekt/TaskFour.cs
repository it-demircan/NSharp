using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp;
using NSharp.Numerics.DG;
using NSharp.Converter;

namespace TaskManagement.FirstProjekt
{
    /// <summary>
    /// Aufgabe 4 - Interpolation und Visualisierung
    /// </summary>
    class TaskFour
    {
        const int N = 10;
        const int N_OUT = 1000;
        const double LEFT_BOARDER = -1.0;
        const double RIGHT_BOARDER = 1.0;

        /// <summary>
        /// Berechnet das Ergebnis der Visualisierungsmatrix und der maßgebenden Auswertungsvektoren und konvertiert dies in ein MatLab String
        /// </summary>
        public void evaluate()
        {
            evaluateFunctionAndGenerateMatLabPlotString(firstFunction);
            evaluateFunctionAndGenerateMatLabPlotString(secondFunction);
        }

        private void evaluateFunctionAndGenerateMatLabPlotString(Func<double,double> function)
        {
            Vector gaussLobattoNodes, weights, nodes; 
            LegendrePolynomEvaluator.computeGaussLobattoNodesAndWeights(N, out gaussLobattoNodes, out weights);
            Vector evaluation = evaluateFunctionAtNodes(function, gaussLobattoNodes);
            nodes = computeAquiDistantNodes(LEFT_BOARDER, RIGHT_BOARDER, N_OUT);
            Matrix visualizationMatrix = computeVisualizationMatrix(evaluation, gaussLobattoNodes, nodes);

            Vector visualizedEvaluation             = visualizeFunction(visualizationMatrix, evaluation);
            Vector visualizedDerivativeEvaluation   = visualizeFunctionDerivative(visualizationMatrix, evaluation, gaussLobattoNodes);

            String evaluationMatLabString = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(nodes, visualizedEvaluation, "X", "Pn - Lagrange Darstellung", "Visualisierung der Lagrange Interpolation anhand cos(x)");
            String derivativeEvaluationMatLabString = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(nodes, visualizedDerivativeEvaluation, "X", "(d/dx)Pn - Ableitung Lagrange Darstellung", "Visualisierung der Ableitung der Lagrange Interpolation anhand 1/(1+x^2)");
        }

        private Vector visualizeFunction(Matrix visualizationMatrix, Vector evaluation)
        {
            Vector visualizedEvaluation = visualizationMatrix * evaluation;
            return visualizedEvaluation;
        }

        private Vector visualizeFunctionDerivative(Matrix visualizationMatrix, Vector evaluation, Vector nodes)
        {
            Matrix diffMatrix = InterpolationToolbox.computeLagrangePolynomeDerivativeMatrix(nodes);
            Vector differentialFunctionEvaluation = diffMatrix * evaluation;
            Vector visualizedDerivativeEvaluation = visualizationMatrix * differentialFunctionEvaluation;
            return visualizedDerivativeEvaluation;
        }

        /// <summary>
        /// Aquidistante Stützstellen berechnen
        /// </summary>
        /// <param name="leftBorder">Linker Rand</param>
        /// <param name="rightBorder">Rechter Rand</param>
        /// <param name="N">Zerlegungsparameter</param>
        /// <returns>Vektor mit aquidistanten Stützstellen</returns>
        private Vector computeAquiDistantNodes(double leftBorder, double rightBorder, int N)
        {
            Vector aquiDistantNodes = new Vector(N + 1);

            for(int i = 0; i <= N; i++)
            {
                aquiDistantNodes[i] = leftBorder + ((rightBorder - leftBorder) / (double)N) * (double)i;
            }
            return aquiDistantNodes;
        }

        /// <summary>
        /// Berechnet die Visualisierungsmatrix
        /// </summary>
        /// <param name="lagrangeEvaluation">Auswertung der Funktion an den erzeugenden Stützstellen</param>
        /// <param name="lagrangeNodes">Lagrange Polynome erzeugende Stützstellen </param>
        /// <param name="nodes">Auswertungsstellen der Interpolation</param>
        /// <returns></returns>
        private Matrix computeVisualizationMatrix(Vector lagrangeEvaluation, Vector lagrangeNodes, Vector nodes)
        {
            Matrix visualMatrix = new Matrix(nodes.Length, lagrangeNodes.Length);

            for(int i = 0; i < nodes.Length; i++)
            {
                for(int k = 0; k < lagrangeNodes.Length; k++)
                {
                    visualMatrix[i, k] = InterpolationToolbox.evaluateLagrangePolynome(nodes[i], lagrangeNodes, k);
                }
            }
            return visualMatrix;
        }

        /// <summary>
        /// Evaluiert eine Funktion an gegebenen Stützstellen.
        /// </summary>
        /// <param name="function">Delegate an die Funktion</param>
        /// <param name="nodes">Die Stützstellen</param>
        /// <returns>Funktionsauswertung an den Stützstellen</returns>
        private Vector evaluateFunctionAtNodes(Func<double, double> function, Vector nodes)
        {
            Vector evaluation = new Vector(nodes.Length);
            for (int i = 0; i < nodes.Length; i++)
                evaluation[i] = function(nodes[i]);

            return evaluation;

        }

        private double firstFunction(double x)
        {
            return Math.Cos(x);
        }

        private double secondFunction(double x)
        {
            return 1.0 / (1.0 + x * x);
        }

    }
}
