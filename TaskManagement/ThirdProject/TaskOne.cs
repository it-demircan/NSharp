using NSharp.Numerics.DG;
using NSharp.Numerics.DG._1DSystem;
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
            controller.createDGElements(8, 3, 0.0, 1.0, 2);
            controller.ComputeSolution(1.0, 0.001);
            Console.Write(controller.GetSolution().toString(15));
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
