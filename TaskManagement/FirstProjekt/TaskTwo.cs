using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.DG;

namespace TaskManagement.FirstProjekt
{
    /// <summary>
    /// Aufgabe 2 - Vergleich der numerischen Integration mit GL und Legendre Stützstellen
    /// </summary>
    class TaskTwo
    {
        static int N = 5;
        static int[] elements = {5, 10, 20};

        public void evaluate()
        {
            Console.WriteLine("############################");
            Console.WriteLine("Starte Aufgabe 2....");
            foreach(int i in elements)
            {
                N = i;
                Console.WriteLine("Gauß-Legendre with N = " + i);
                evaluateGaussLobattoIntegration();
                Console.WriteLine("Gauß-Lobatto with N = " + i);
                evaluateGaussLegendreIntegration();
            }
            Console.WriteLine("Aufgabe 2 abgeschlossen.");
            Console.WriteLine("############################");
        }

        private void evaluateGaussLegendreIntegration()
        {
            double result;
            Console.WriteLine("Erste Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(firstFunction, N);
            Console.WriteLine(result);
            
            Console.WriteLine("Zweite Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(secondFunction, N);
            Console.WriteLine(result);
            
            Console.WriteLine("Dritte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(thirdFunction, N);
            Console.WriteLine(result);
            
            Console.WriteLine("Vierte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(fourthFunction, N);
            Console.WriteLine(result);

            Console.WriteLine("Fünfte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussNodesAndWeights(fifthFunction, N);
            Console.WriteLine(result);
        }


        private void evaluateGaussLobattoIntegration()
        {
            double result;
            Console.WriteLine("Erste Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussLobattoNodesAndWeights(firstFunction, N);
            Console.WriteLine(result);

            Console.WriteLine("Zweite Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussLobattoNodesAndWeights(secondFunction, N);
            Console.WriteLine(result);

            Console.WriteLine("Dritte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussLobattoNodesAndWeights(thirdFunction, N);
            Console.WriteLine(result);

            Console.WriteLine("Vierte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussLobattoNodesAndWeights(fourthFunction, N);
            Console.WriteLine(result);

            Console.WriteLine("Fünfte Funktion mit N = " + N);
            result = IntegrationToolbox.computeGaussianIntegrationWithGaussLobattoNodesAndWeights(fifthFunction, N);
            Console.WriteLine(result);
        }


        public double constFunction(double x)
        {
            return 1;
        }
        private double firstFunction(double x)
        {
            return Math.Cos(x);
        }

        private double secondFunction(double x)
        {
            return 1.0 / (1.0 + x * x);
        }

        private double thirdFunction(double x)
        {
            return Math.Pow(x, 2.0 * N - 2.0);
        }

        private double fourthFunction(double x)
        {
            return Math.Pow(x, 2 * N);
        }


        private double fifthFunction(double x)
        {
            return Math.Pow(x, 2 * N + 2);
        }
    }
}
