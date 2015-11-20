using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Numerics.Interpolation
{
    public interface IInterpolater
    {
        double evaluateInterpolation(double x,Vector nodes, Vector functionValues);

    }
}
