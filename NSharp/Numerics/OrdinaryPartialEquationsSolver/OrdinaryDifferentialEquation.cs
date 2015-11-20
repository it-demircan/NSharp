using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.OrdinaryPartialEquationsSolver
{
    public class OrdinaryDifferentialEquation
    {
        List<Func<double,double, double>> myOrdinaryFunctionList;

        public OrdinaryDifferentialEquation(List<Func<double,double, double>> ordinaryFunctionList)
        {
            myOrdinaryFunctionList = ordinaryFunctionList;
        }

        public Vector evaluateOrdinaryDifferentialEquation(Vector startSolution, double time)
        {
            Vector evaluatedODE = new Vector(myOrdinaryFunctionList.Count);

            for(int i = 0; i < myOrdinaryFunctionList.Count; i++)
            {
                evaluatedODE[i] = myOrdinaryFunctionList.ElementAt(i)(startSolution[i], time);
            }
            return evaluatedODE;
        }
    }
}
