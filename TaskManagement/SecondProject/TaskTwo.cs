using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;
using NSharp.Numerics.OrdinaryPartialEquationsSolver;

namespace TaskManagement.SecondProject
{
    class TaskTwo
    {
        public void evaluate()
        {
            //Anfangsbedingung
            Vector initial = new Vector(1);
            initial[0] = -1.0;
            List<Func<double, double, double>> mySystem = new List<Func<double, double, double>>();
            mySystem.Add(derivativeFunction);
            OrdinaryDifferentialEquation testEquation = new OrdinaryDifferentialEquation(mySystem);
            IODESolver odeSolver = new RungeKuttaSolver();

            double[] errorList = new double[10];

            for (int i = 1; i <= 10; i++)
            {
                double timeStep = Math.Pow(2.0, -i);
                Vector res = odeSolver.computeSolutionVectorWithMultipleSteps(initial, testEquation, 1.0, 2.0, timeStep);
                //Vector exact = null; //hier die exacte auswertung
                //double error = (res - exact) * (res - exact);
            }
            
        }

        /// <summary>
        /// Diese Funktion bildet die DGL ab, mit u ausgewertet und t als Ortsvariable
        /// Beispiel: u'(x) = u(x) * x
        /// Somit wäre evaluatedFunction u(x), also bekannt
        /// und x wäre time.
        /// Die Lösung hängt insbesondere von unserer Initial Vektoren ab.
        private double derivativeFunction(double evaluatedFunction, double time){
            return evaluatedFunction* evaluatedFunction;
        }
    }
}
