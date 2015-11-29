using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    class DGController
    {
        List<DGElement> elements;
        int polynomOrder; //N
        static double velocity = 2.0;

        public void createDGElements(int N, double LeftBoundary, double RightBoundary)
        {            
            for (int i = 0; i<N; i++)
            {
                elements[i] = new DGElement(LeftBoundary + (double)i * 1.0 / (double)N, LeftBoundary + (double)(i + 1) * 1.0 / (double)N, polynomOrder, fluxFunction, leftBoundaryFunction, rightBoundaryFunction, NumFlux, inhomogenuousPart);

                if(i>0)
                    elements[i].LeftNeighbour = elements[i - 1];
                if(i<N-1)
                    elements[i].RightNeighbour = elements[i + 1];
            }
        }

        private double leftBoundaryFunction(double x)
        {
            return 0.0;
        }

        private double rightBoundaryFunction(double x)
        {
            return 0.0;
        }
        private double inhomogenuousPart(double space, double time)
        {
            return 0.0;
        }

        private double fluxFunction(double x)
        {
            return velocity * x;
        }

        private double NumFlux(double left, double right)
        {
            return left;
        }

    }
}
