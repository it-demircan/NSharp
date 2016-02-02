using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.DG._2DSystem;
using Structures;

namespace TaskManagement.FourthProject
{
    public class TaskTwo
    {



        public static void Test()
        {
            
            DGController2D myController = new DGController2D();

            int maxErrorCalculations = 4;
            Matrix error = new Matrix(maxErrorCalculations, 3);
            int polyOrder = 6;
            double endTime = 1.0;
            double deltaT = 0.01;

            for (int i = 0; i < maxErrorCalculations; i++)
            {
                myController.Init(polyOrder, (int)Math.Pow(2.0,i+1), (int)Math.Pow(2.0, i+1),0.5);
                myController.ComputeSolution(endTime);

                //Console.WriteLine((v0 - exact).toString(6));
                double v1Error = myController.ComputeMaxError(endTime, 0);
                double v2Error = myController.ComputeMaxError(endTime, 1);
                double pError = myController.ComputeMaxError(endTime, 2);

                Console.WriteLine("i = " + i + " V1-Error:" + v1Error);                
                Console.WriteLine("i = " + i + " V2-Error:" + v2Error);
                Console.WriteLine("i = " + i + " p-Error:" + pError);

                error[i, 0] = v1Error;
                error[i, 1] = v2Error;
                error[i, 2] = pError;
            }

            Matrix eoc = computeEOC(error);

            Console.WriteLine();
            Console.WriteLine(eoc.toString(5));
            string errorStr = error.toString(15);
            string eocStr = eoc.toString(15);
        }


        private static Matrix computeEOC(Matrix errorMatrix)
        {
            Matrix EOC = new Matrix(errorMatrix.NoRows - 1, errorMatrix.NoColumns);
            for (int k = 0; k < errorMatrix.NoColumns; k++)
            {
                for (int i = 0; i < errorMatrix.NoRows - 1; i++)
                {
                    EOC[i, k] = Math.Log(errorMatrix[i + 1, k] / errorMatrix[i, k]) / Math.Log(1.0 / 2.0);
                }
            }
            return EOC;
        }





        public static void TestStart()
        {
            int NQ = 50, MQ = 50;

            _2DElement[] elements = new _2DElement[NQ * MQ];
            double yBottom = 0.0;
            double yTop = 1.0;
            double xLeft = 0.0;
            double xRight = 1.0;
            int N = 3;
            int sysDim = 3;

            for(int i = 0;  i < NQ; i++)
            {
                for(int k = 0; k < MQ; k++)
                {
                    elements[i*MQ+k] = new _2DElement(sysDim, N, (yTop- yBottom)/MQ * (double)k , (yTop - yBottom) / MQ * (k+1.0), (xRight - xLeft) / NQ * (double)i, (xRight - xLeft) / NQ * (i+1.0));
                }
            }

            //Nachbarn setzten
            for (int i = 0; i < NQ; i++)
            {
                for (int k = 0; k < MQ; k++)
                {
                    _2DElement recentCell = elements[i * MQ + k];
                    if(i == 0)
                    {
                        recentCell.Left = elements[(NQ - 1) * MQ + k];
                        recentCell.Right = elements[(i + 1) * MQ + k];
                    }
                    else if(i == NQ-1)
                    {
                        recentCell.Left = elements[(i - 1) * MQ + k];
                        recentCell.Right = elements[k];
                    }
                    else
                    {
                        recentCell.Left =  elements[(i-1) * MQ + k];
                        recentCell.Right = elements[(i + 1) * MQ + k];
                    }

                    if(k == 0)
                    {
                        recentCell.Bottom = elements[i * MQ + (MQ - 1)];
                        recentCell.Top = elements[i * MQ + k+1];
                    }
                    else if( k == MQ - 1){
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

            Matrix[] res = elements[0].ComputeTimeEvaluation();
        }
    }
}
