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
            DGSystemController controller = new DGSystemController();
            controller.createDGElements(20, 3, 0.0, Math.PI*2.0, 2);
            controller.ComputeSolution(1.0,0.001);
            Matrix approx = controller.GetSolution();
            Matrix exact = controller.ComputeExactSolution(1.0);

            Console.WriteLine(approx.toString(15));
            Console.WriteLine("Diff:");
            Vector diff = approx.GetColumn(0) - exact.GetColumn(0);
            Console.Write(diff.toString(15));
            double err = Math.Sqrt(diff * diff);
            Console.Write("Error:" + err);
            Console.Write(approx.toString(15));
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
