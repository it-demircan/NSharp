using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.OrdinaryPartialEquationsSolver
{
    public interface IODESolver
    {
        Vector computeSolutionForNextStep(Vector initial, OrdinaryDifferentialEquation ode, double startTime, double endTime);
        Vector computeSolutionWithMultipleSteps(Vector initial, OrdinaryDifferentialEquation ode, double startTime, double endTime, double step);

        Vector computeSolutionVectorWithMultipleSteps(Vector initial, OrdinaryDifferentialEquation ode, double startTime, double endTime, double step);
    }
}
