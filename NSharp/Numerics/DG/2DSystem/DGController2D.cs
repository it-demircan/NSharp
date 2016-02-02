using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Numerics.DG._2DSystem
{
    public enum TaskNr
    {
        TaskTwo, TaskThree
    }
    public class DGController2D
    {
        public int NQ = 8;
        public int MQ = 8;

        public double XLeft = 0.0;
        public double XRight = 1.0;
        public double YTop = 1.0;
        public double YBottom = 0.0;
        public int N = 5;
        public int SysDim = 3;

        public DGElement2D[] elements;

        public double CFL = 0.5;

        public void Init(int N, int NQ, int MQ, double CFL = 0.5)
        {
            this.N = N;
            this.NQ = NQ;
            this.MQ = MQ;
            this.CFL = CFL;

            CreateMesh();
            SetNeighboursPeriodic();
        }

        public void ComputeSolution(double endTime, double timeStep = 0.0)
        {
            double recentTime = 0.0;
            Matrix3D[] recentTimeDerivatives = new Matrix3D[elements.Length];
            double recentTimeStep = timeStep;

            while (recentTime < endTime)
            {
                if (timeStep == 0.0)
                {
                    double lambdaMax = 1.0;
                    recentTimeStep = ComputeTimeStep(CFL, lambdaMax);
                }

                if (recentTime + recentTimeStep > endTime)
                    recentTimeStep = endTime - recentTime;

                for (int i = 0; i < elements.Length; i++)
                    recentTimeDerivatives[i] = new Matrix3D(N + 1, N+1,  SysDim);

                //Runge Kutte
                for (int k = 0; k < 5; k++)
                {
                    for (int i = 0; i < elements.Length; i++)
                    {
                        Matrix3D solutionSystem = elements[i].Solution;
                        double nextTimeStep = recentTime + B[k] * recentTimeStep;
                        Matrix3D EvaluatedTimeDerivative = elements[i].ComputeTimeDerivative(nextTimeStep);

                        for (int sysIdx = 0; sysIdx < SysDim; sysIdx++)
                        {
                            Matrix tempTimeDerivative = A[k] * recentTimeDerivatives[i][sysIdx] + EvaluatedTimeDerivative[sysIdx];
                            recentTimeDerivatives[i][sysIdx] = tempTimeDerivative;
                            Matrix solution = solutionSystem[sysIdx] + C[k] * recentTimeStep * recentTimeDerivatives[i][sysIdx];
                            solutionSystem[sysIdx] = solution;
                        }

                        elements[i].Solution = solutionSystem;
                    }
                    //Rand updaten
                    for (int i = 0; i < elements.Length; i++)
                    {
                        elements[i].UpdateBorderValues();
                    }
                }
                recentTime += recentTimeStep;
            }
        }

        public void CreateMesh()
        {
            elements = new DGElement2D[NQ * MQ];
            for (int i = 0; i < NQ; i++)
            {
                for (int k = 0; k < MQ; k++)
                {
                    elements[i * MQ + k] = new DGElement2D(N,
                        SysDim,
                        (XRight - XLeft) / NQ * (double)i,
                        (XRight - XLeft) / NQ * (i + 1.0),
                        (YTop - YBottom) / MQ * (double)k,
                        (YTop - YBottom) / MQ * (k + 1.0),
                        StartSolution,
                        FluxF,
                        FluxG,
                        NumFluxF,
                        NumFluxG,
                        ExactSolution
                        );
                }
            }
        }

        public double ComputeMaxError(double time, int sysIdx)
        {
            double maxError = 0.0;

            for(int i = 0; i< elements.Length; i++)
            {
                double recentError = elements[i].GetMaxErrorAtTimeForSystemEquation(time, sysIdx);
                if (maxError < recentError)
                    maxError = recentError;
            }

            return maxError;
        }

        private void SetNeighboursPeriodic()
        {
            for (int i = 0; i < NQ; i++)
            {
                for (int k = 0; k < MQ; k++)
                {
                    DGElement2D recentCell = elements[i * MQ + k];
                    if (i == 0)
                    {
                        recentCell.Left = elements[(NQ - 1) * MQ + k];
                        recentCell.Right = elements[(i + 1) * MQ + k];
                    }
                    else if (i == NQ - 1)
                    {
                        recentCell.Left = elements[(i - 1) * MQ + k];
                        recentCell.Right = elements[k];
                    }
                    else
                    {
                        recentCell.Left = elements[(i - 1) * MQ + k];
                        recentCell.Right = elements[(i + 1) * MQ + k];
                    }

                    if (k == 0)
                    {
                        recentCell.Bottom = elements[i * MQ + (MQ - 1)];
                        recentCell.Top = elements[i * MQ + k + 1];
                    }
                    else if (k == MQ - 1)
                    {
                        recentCell.Top = elements[i * MQ];
                        recentCell.Bottom = elements[i * MQ + k - 1];
                    }
                    else
                    {
                        recentCell.Top = elements[i * MQ + k + 1];
                        recentCell.Bottom = elements[i * MQ + k - 1];
                    }
                }
            }
        }

        private Vector StartSolution(Vector spaceNodes)
        {
            //Aufgabe 2
            //return ExactSolution(spaceNodes, 0.0);

            //Aufgabe 3
            Vector eva = new Vector(3);
            double x = spaceNodes[0];
            double y = spaceNodes[1];

            eva[0] = 0.0;
            eva[1] = 0.4 <= x && x <= 0.6 && 0.1 <= y && y <= 0.4 ? 1.0 : 0.0;
            eva[2] = 0.4 <= x && x <= 0.6 && 0.1 <= y && y <= 0.4 ?
                3.0 * Math.Exp(-1.0 / 2.0 * ((((x - 0.5) * (x - 0.5)) + ((y - 0.25) * (y - 0.25))) / 0.01)) + 2 :
                2.0;

            return eva;
        }

        private Vector FluxF(Vector u)
        {
            Vector flux = new Vector(SysDim);
            flux[0] = u[2];
            flux[2] = u[0];
            return flux;
        }
        private Vector FluxG(Vector u)
        {
            Vector flux = new Vector(SysDim);
            flux[1] = u[2];
            flux[2] = u[1];
            return flux;
        }

        private Vector NumFluxF(Vector left, Vector right)
        {
            Vector numFlux = FluxF(left) + FluxF(right) + left - right;
            numFlux = (Vector)((1.0 / 2.0) * numFlux);
            return numFlux;
        }

        private Vector NumFluxG(Vector bottom, Vector top)
        {
            Vector numFlux = FluxG(bottom) + FluxG(top) + bottom - top;
            numFlux = (Vector)((1.0 / 2.0) * numFlux);
            return numFlux;
        }

        private Vector ExactSolution(Vector spaceNodes, double time)
        {
            Vector eva = new Vector(3);
            double x = spaceNodes[0];
            double y = spaceNodes[1];

            eva[0] = 2.0 + Math.Sin(2.0 * Math.PI * (x - time));
            eva[1] = 4.0 + Math.Cos(2.0 * Math.PI * (y - time));
            eva[2] = eva[0] + eva[1];

            return eva;
        }

        public double ComputeTimeStep(double cfl, double lambdaMax)
        {
            return cfl * ComputeSurfaceInElement() / lambdaMax;
        }

        public double ComputeSurfaceInElement()
        {
            return (1.0/2.0)* (((XRight - XLeft) / (double)NQ) / (double)(N + 1));
        }

        #region Runge Kutte 
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

        #endregion

    }
}
