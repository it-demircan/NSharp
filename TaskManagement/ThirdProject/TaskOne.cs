using NSharp.Converter;
using NSharp.Numerics.DG;
using NSharp.Numerics.DG._1DSystem;
using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskManagement.ThirdProject
{
    class TaskOne
    {

        public void TestSystemDG()
        {
            TaskOneStone();
            //TaskOneEOC();
            //TaskOneMass();
            //TaskOneEnergy();
            //TaskOneWellBalanced();
        }

        public void TaskOneStone()
        {
            DGSystemController controller = new DGSystemController();
            controller.createDGElements(4, 3, 0.0, 20.0, 2);
            Matrix Solution = controller.GetSolution();
            Vector startV = ComputeV(Solution);
            string matlabStringStartH = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(controller.GetOriginNodes(), Solution.GetColumn(0), "Raumvariable x", "h", "Lösung h zum Zeitpunkt t = 0.0");
            string matlabStringStartV = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(controller.GetOriginNodes(), startV, "Raumvariable x", "v", "Lösung v zum Zeitpunkt t = 0.0");

            controller.ComputeSolution(2.3);
            Solution = controller.GetSolution();

            Vector endV = ComputeV(Solution);
            string matlabStringV = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(controller.GetOriginNodes(), endV, "Raumvariable x", "v", "Lösung v zum Zeitpunkt t = 2.3");

            string matlabStringH = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(controller.GetOriginNodes(), Solution.GetColumn(0), "Raumvariable x", "h", "Lösung h zum Zeitpunkt t = 2.3");

        }

        private Vector ComputeV(Matrix Solution)
        {
            Vector v = new Vector(Solution.NoRows);
            Vector h = Solution.GetColumn(0);
            Vector hv = Solution.GetColumn(1);

            for(int i = 0; i < v.Length; i++){
                v[i] = hv[i] / h[i];
            }
            return v;
        }

        public void TaskOneEnergy()
        {
            DGSystemController controller = new DGSystemController();
            controller.createDGElements(101, 3, 0.0, 20.0, 2);

            double test = controller.ComputeEnergy();
            Console.WriteLine(test);
            Dictionary<double, double> Energy = controller.ComputeSolution(1.0);


            Vector time = new Vector(Energy.Count);
            Vector energyVector = new Vector(Energy.Count);
            int counter = 0;
            foreach (KeyValuePair<double, double> entry in Energy)
            {
                time[counter] = entry.Key;
                energyVector[counter] = entry.Value;
                counter++;
                Console.WriteLine("Time:" + entry.Key + " Energy:" + entry.Value);
            }

            string matlabString = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(time, energyVector, "Zeit t", "Energie e", "Energieverlauf");

        }


        private void TaskOneWellBalanced()
        {
            DGSystemController controller = new DGSystemController();
            controller.createDGElements(100, 3, 0.0, 20.0, 2);
            Dictionary<double, double> Energy = controller.ComputeSolution(1.0);

            Console.WriteLine(controller.ComputeConstant().toString(15));
            Vector nodes = controller.GetOriginNodes();
            Vector constant = controller.ComputeConstant();

            string matlabString = MatLabConverter.ConvertToMatLabPlotStringWithAxisLabelAndTitle(nodes, constant, "Ort x", "h+b", "Well-Balanced Test");
        }
        private void TaskOneMass()
        {
            DGSystemController controller = new DGSystemController();            
            controller.createDGElements((int)Math.Pow(2.0, 5), 3, 0.0, Math.PI * 2.0, 2);
            double massBefore = controller.ComputeMass()[0];
            Console.WriteLine("Masse an t = 0: " + massBefore);
            controller.ComputeSolution(1.0);
            double massAfter = controller.ComputeMass()[0];
            Console.WriteLine("Masse an t = 1.0: " + massBefore);
        }
        private void TaskOneEOC()
        {
            Matrix error = new Matrix(9, 1);
            for (int i = 1; i < 8; i++)
            {
                DGSystemController controller = new DGSystemController();
                controller.createDGElements((int)Math.Pow(2.0, i), 3, 0.0, Math.PI * 2.0, 2);
                controller.ComputeSolution(1.0);
                Matrix approx = controller.GetSolution();
                Matrix exact = controller.ComputeExactSolution(1.0);

                Vector diff = approx.GetColumn(0) - exact.GetColumn(0);
                double err = diff.GetMaxAbsValue();
                Console.WriteLine("Error:" + err + "   i = " + i);
                error[i - 1, 0] = err;
            }

            Matrix eoc = computeEOC(error);
            Console.WriteLine("EOC:");
            Console.WriteLine(eoc.toString(15));
        }

        private Matrix computeEOC(Matrix errorMatrix)
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

        public void TestOldDG()
        {
            DGController controller = new DGController();
            controller.createDGElements(8, IntegrationMode.GaussLobatto, 3, 0.0, 1.0);
            controller.computeSolution(1.0, 0.001);
            Console.Write(controller.getCompleteSolution().toString(15));
        }
    }
}
