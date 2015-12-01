using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG
{
    public class DGController
    {

        //Variablen für Runga Kutta
        double[] A = { 0.0,
                -567301805773.0/1357537059087.0,
                -2404267990393.0/ 2016746695238.0,
                - 3550918686646.0/ 2091501179385.0,
                -1275806237668.0/ 842570457699.0};

        double[] B = { 0.0,
            1432997174477.0 / 9575080441755.0,
            2526269341429.0/ 6820363962896.0,
            2006345519317.0/ 3224310063776.0,
            2802321613138.0/ 2924317926251.0 };

        double[] C = { 1432997174477.0/ 9575080441755.0,
            5161836677717.0 / 13612068292357.0,
            1720146321549.0 / 2090206949498.0,
            3134564353537.0 / 4481467310338.0,
            2277821191437.0 / 14882151754819.0};

        DGElement[] elements;
        int polynomOrder; //N

        public void createDGElements(int numberOfDGElements, int polynomOrder, double leftBoundary, double rightBoundary)
        {
            this.polynomOrder = polynomOrder;
            elements = new DGElement[numberOfDGElements];
            for (int i = 0; i<numberOfDGElements; i++)
            {
                double leftSpaceBorder = leftBoundary + (double)i * (rightBoundary - leftBoundary) / (double)numberOfDGElements;
                double rightSpaceBorder = leftBoundary + (double)(i + 1) * (rightBoundary-leftBoundary) / (double)numberOfDGElements;
                elements[i] = new DGElement(leftSpaceBorder, rightSpaceBorder, polynomOrder, FluxFunction, NumFlux, InhomogenuousPart, InitialFunction);

                if (i > 0)
                {
                    elements[i].LeftNeighbour = elements[i - 1];
                    elements[i - 1].RightNeighbour = elements[i];
                }
            }

            //Periodic Boundary Condition
            if (numberOfDGElements > 1)
            {
                elements[0].LeftNeighbour = elements[numberOfDGElements - 1];
                elements[numberOfDGElements - 1].RightNeighbour = elements[0];
            }
            else
            {
                elements[0].LeftNeighbour = elements[0];
                elements[0].RightNeighbour = elements[0];
            }
        }


        public void computeSolution(double endTime, double timeStep)
        {
            double recentTime = 0.0;
            Vector[] tempTimeDerivaties = new Vector[elements.Length];

            do
            {
                //Vektoren zurüksetzen
                for (int i = 0; i < elements.Length; i++)
                    tempTimeDerivaties[i] = new Vector(polynomOrder+1);

               //Runga Kutta
                for (int k = 0; k < 5; k++)
                {
                    //Runga Kutta für jedes Element, damit Zeitsynchron berechnet wird
                    for (int i = 0; i < elements.Length; i++)
                    {                        
                        Vector solution = elements[i].GetSolution();
                        double nextTimeStep = recentTime + B[k] * timeStep;
                        tempTimeDerivaties[i] = (Vector)(A[k] * tempTimeDerivaties[i]) + elements[i].EvaluateTimeDerivative(nextTimeStep);
                        solution = solution + (Vector)(C[k] * timeStep * tempTimeDerivaties[i]);
                        elements[i].UpdateSolution(solution);
                    }
                    //Rand updaten
                    for (int i = 0; i < elements.Length; i++)
                    {
                        elements[i].UpdateBorderValues();
                    }
                }

                recentTime += timeStep;
            }
            while (recentTime < endTime);

            /**
                Hier koennte man das letzte Stuck noch erganzen.
            **/
            //recentTime -= timeStep;

            //timeStep = endTime - recentTime;
            ////Vektoren zurüksetzen
            //for (int i = 0; i < elements.Length; i++)
            //    tempTimeDerivaties[i] = new Vector(polynomOrder + 1);


            //for (int k = 0; k < 5; k++)
            //{

            //    for (int i = 0; i < elements.Length; i++)
            //    {
            //        Vector solution = elements[i].getSolution();
            //        double nextTimeStep = recentTime + B[k] * timeStep;
            //        tempTimeDerivaties[i] = (Vector)(A[k] * tempTimeDerivaties[i]) + elements[i].EvaluateTimeDerivative(nextTimeStep);
            //        solution = solution + (Vector)(C[k] * timeStep * tempTimeDerivaties[i]);
            //        elements[i].updateSolution(solution);
            //    }

            //    for (int i = 0; i < elements.Length; i++)
            //    {
            //        elements[i].UpdateBorderValues();
            //    }
            //}

            for (int k = 0; k < elements.Length; k++)
            {
                for (int i = 0; i < elements[k].GetSolution().Length; i++)
                    Console.WriteLine(elements[k].GetSolution()[i]);
            }
        }

        private double InitialFunction(double space)
        {
            //return Math.Exp(-Math.Log(2.0) * (space + 1.0) * (space + 1.0) / (0.2 * 0.2));
            return Math.Sin(Math.PI*space);
        }

        private double InhomogenuousPart(double space, double time)
        {
            return 0.0;
        }

        private double FluxFunction(double u)
        {
            return 1.0 * u;
        }

        private double NumFlux(double left, double right)
        {
            return 1.0 / 2.0 * (FluxFunction(left) + FluxFunction(right)) - 1.0 / 12.0 *((right - left) * (right - left));
        }

    }
}
