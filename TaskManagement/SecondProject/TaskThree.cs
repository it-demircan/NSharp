using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSharp.Numerics.DG;

namespace TaskManagement.SecondProject
{
    class TaskThree
    {
        public void evaluate()
        {
            try {
                DGController dgController = new DGController();

                dgController.createDGElements(256, 3, 0, 1);
                dgController.computeSolution(0.5, 0.001);
            }catch(Exception err)
            {
                Console.Write(err.StackTrace);
            }
        }
    }
}
