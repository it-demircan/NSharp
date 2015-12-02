using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.DG;
using Structures;

namespace TaskManagement.SecondProject
{
    class TaskThree
    {
        static int[] elementNumber  = {8,16,32};
        static int[] polynomOrders = {3,7};
        static IntegrationMode myMode = IntegrationMode.GaussLobatto;
        static double leftSpaceBorder = 0.0;
        static double rightSpaceBorder = 1.0;
        static double endTime = 1.0;
        static double timeStep = 0.001;

        public void evaluate()
        {
            try {
                Matrix errorList = computeErrorList();
                Console.WriteLine(errorList.toString(15));

                Matrix EOC = computeEOC(errorList);
                Console.WriteLine(EOC.toString(15));
            }catch(Exception err)
            {
                Console.Write(err.StackTrace);
            }
        }


        private Matrix computeErrorList()
        {
            Matrix errorList = new Matrix(elementNumber.Length, polynomOrders.Length);
            for (int i = 0; i < polynomOrders.Length; i++)
            {
                int polynomOrder = polynomOrders[i];
                for (int k = 0; k < elementNumber.Length; k++)
                {
                    DGController dgController = new DGController();
                    dgController.createDGElements(elementNumber[k], myMode, polynomOrder, leftSpaceBorder, rightSpaceBorder);
                    errorList[k, i] = dgController.computeSolution(endTime, timeStep);
                }
            }

            return errorList;
        }

        private Matrix computeEOC(Matrix errorMatrix)
        {
            Matrix EOC = new Matrix(errorMatrix.NoRows - 1, errorMatrix.NoColumns);
            for (int k = 0; k < errorMatrix.NoColumns; k++)
            {
                for (int i = 0; i < errorMatrix.NoRows - 1; i++)
                {
                    EOC[i, k] = Math.Log(errorMatrix[i + 1,k] / errorMatrix[i,k]) / Math.Log(1.0 / 2.0);
                }
            }
            return EOC;
        }
    }
}
