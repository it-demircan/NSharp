using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.DG;
using System.IO;
using Structures;
using System.Diagnostics;

namespace TaskManagement.SecondProject
{
    class TaskThree
    {
        static double leftSpaceBorder = 0.0;
        static double rightSpaceBorder = 1.0;
        static double endTime = 10.0;
        static double timeStep; //N = 7, 0.0001 nehmen.
        public void evaluate()
        {
            try {
                ComputeSolutionForUnsteady();
            }catch(Exception err)
            {
                computeDGLMatrices();
            }
        }

        private Vector GetCFL(IntegrationMode mode)
        {
            Vector cfl;
            if (mode == IntegrationMode.GaussLegendre)
            {
                cfl = new Vector(new double[] { 2, 1.4, 1.08, 0.9, 0.78, 0.68, 0.62 });
            }
            else
            {
                cfl = new Vector(new double[] { 1.5, 1.3, 1.2, 1.05, 0.95, 0.9, 0.83 });
            }
            return cfl;
        }

        public void analyseCellElements(double error)
        {
            Console.WriteLine("Start computation for No. Cell Elements...");
            IntegrationMode mode = IntegrationMode.GaussLobatto;
            int[] polynomOrders = { 1, 3, 7 };
            Vector CFLMapping = GetCFL(mode);
            double computedError = 0.0;
            int noElements ;

            for (int i = 0; i < polynomOrders.Length; i++)
            {
                noElements = 1;
                int polynomOrder = polynomOrders [i];
                do{
                    Console.Write("N = " + polynomOrder + " - N_Q = " + noElements);
                    Stopwatch sw = new Stopwatch();
                    sw.Start();
                    DGController dgController = new DGController();
                    dgController.createDGElements(noElements, mode, polynomOrder, leftSpaceBorder, rightSpaceBorder);
                    timeStep = dgController.ComputeTimeStep(CFLMapping[polynomOrders[i] - 1]);
                    computedError = dgController.computeSolution(endTime, timeStep);
                    sw.Stop();
                    Console.Write(" - L2 Error: " +);
                    Console.WriteLine(" - PassedTime: " + sw.Elapsed);

                    noElements++;
                }while(computedError > error);
            }

            Console.WriteLine("Analyse Cell Elements finished");
        }

        #region ErrorList
        private void ComputeErrorLists()
        {
            Console.Write("Gauss Lobatto:");
            Matrix errorList = computeErrorList(IntegrationMode.GaussLobatto);
            GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\GaussLobatto_ErrorList.txt", NSharp.Converter.MatLabConverter.ConvertMatrixToMatLabMatrix(errorList, "A"));
            Matrix EOC = computeEOC(errorList);
            GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\GaussLobatto_EOCList.txt", NSharp.Converter.MatLabConverter.ConvertMatrixToMatLabMatrix(EOC, "A"));
            Console.Write("Gauss:");
            errorList = computeErrorList(IntegrationMode.GaussLegendre);
            GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\GaussLegendre_Error.txt", NSharp.Converter.MatLabConverter.ConvertMatrixToMatLabMatrix(errorList, "A"));
            EOC = computeEOC(errorList);
            GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\GaussLegendre_EOC.txt", NSharp.Converter.MatLabConverter.ConvertMatrixToMatLabMatrix(EOC, "A"));
        }

        private Matrix computeErrorList(IntegrationMode mode)
        {
            Console.WriteLine("Start computation for Error List...");
            int[] elementNumber = {2,4,8,16,32,64};
            int[] polynomOrders = {1,3,7};

            Vector CFLMapping = GetCFL(mode);
            Matrix errorList = new Matrix(elementNumber.Length, polynomOrders.Length);
            for (int i = 0; i < polynomOrders.Length; i++)
            {
                int polynomOrder = polynomOrders[i];
                for (int k = 0; k < elementNumber.Length; k++)
                {
                    Console.Write("N = " + polynomOrder + " - N_Q = " + elementNumber[k]);
                    Stopwatch sw = new Stopwatch();
                    sw.Start();
                    DGController dgController = new DGController();
                    dgController.createDGElements(elementNumber[k], mode, polynomOrder, leftSpaceBorder, rightSpaceBorder);
                    timeStep = dgController.ComputeTimeStep(CFLMapping[polynomOrders[i]-1]);
                    errorList[k, i] = dgController.computeSolution(endTime, timeStep);
                    sw.Stop();
                    Console.Write(" - L2 Error: " + errorList[k, i]);
                    Console.WriteLine(" - PassedTime: " + sw.Elapsed);
                }
            }

            Console.WriteLine("Error Computation finished");

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
        #endregion 

        #region DGL Matrices
        private void computeDGLMatrices()
        {
            computeDGLMatrices(IntegrationMode.GaussLobatto);
            computeDGLMatrices(IntegrationMode.GaussLegendre);
            Console.Write("DGL Matrices generated...");
        }

        private void computeDGLMatrices(IntegrationMode mode)
        {
            int[] elementNumber = {8,16,32,64};
            int[] polynomOrders = {1,2,3,4,5,6,7};
            
            for (int i = 0; i < polynomOrders.Length; i++)
            {
                int polynomOrder = polynomOrders[i];
                for (int k = 0; k < elementNumber.Length; k++)
                {
                    DGController dgController = new DGController();
                    dgController.createDGElements(elementNumber[k], mode, polynomOrder, leftSpaceBorder, rightSpaceBorder);
                    Matrix A = dgController.ConstructDGLMatrix();
                    string matrixString = NSharp.Converter.MatLabConverter.ConvertMatrixToMatlabReadable(A);
                    GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\" + mode.ToString() + "_NQ=" + elementNumber[k] + "_N =" + polynomOrder+".txt", matrixString);
                }
            }        
        }
        #endregion

        #region Unsteady
        private void ComputeSolutionForUnsteady()
        {
            Console.WriteLine("GaussLobatto");
            ComputeSolutionForUnsteady(IntegrationMode.GaussLobatto);
            Console.WriteLine("GaussLegendre");
            ComputeSolutionForUnsteady(IntegrationMode.GaussLegendre);
        }
        private void ComputeSolutionForUnsteady(IntegrationMode mode)
        {
            Console.WriteLine("Start computation for Unsteady Solution...");
            int[] elementNumber = {50};
            int[] polynomOrders = { 1, 3, 7 };

            Vector CFLMapping = GetCFL(mode);
            Matrix errorList = new Matrix(elementNumber.Length, polynomOrders.Length);
            for (int i = 0; i < polynomOrders.Length; i++)
            {
                int polynomOrder = polynomOrders[i];
                for (int k = 0; k < elementNumber.Length; k++)
                {
                    Console.WriteLine("N = " + polynomOrder + " - N_Q = " + elementNumber[k]);
                    DGController dgController = new DGController();
                    dgController.createDGElements(elementNumber[k], mode, polynomOrder, leftSpaceBorder, rightSpaceBorder);
                    timeStep = dgController.ComputeTimeStep(CFLMapping[polynomOrders[i] - 1]);
                    errorList[k, i] = dgController.computeSolution(endTime, timeStep);

                    Vector space = dgController.getOriginSpace();
                    Vector sol = dgController.getCompleteSolution();
                    string plotString = NSharp.Converter.MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(space, sol, "X-Achse", "u approx", "Approximation mit NQ = " + elementNumber[k] + " N = " + polynomOrder);
                    GeneralHelper.WriteOutputText(Directory.GetCurrentDirectory() + "\\" + mode.ToString() + "_PLOT_N =" + polynomOrder + ".txt", plotString);
                }
            }

            Console.WriteLine("Unsteady Solution finished");

        }
        #endregion 
    }
}
