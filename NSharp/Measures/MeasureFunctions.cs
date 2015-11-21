using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Structures;

namespace NSharp.Measures
{
    public class MeasureFunctions
    {

        /// <summary>
        /// Computes the euclidian distance between two points (x0,y0) and (x1,y1)
        /// </summary>
        /// <returns>Distance</returns>
        public static double CalculateEuclidianDistance(double x0, double y0, double x1, double y1)
        {
            double distance;
            distance = Math.Sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
            return distance;
        }

        public static double CalculateDiscreteNorm(Vector A, Vector B)
        {
            return Math.Sqrt(A * B);
        }
    }
}
