﻿using Structures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NSharp.Converter
{
    public class MatLabConverter
    {

        public static String ConvertToMatLabPlotStringWithAxisLabelAndTitle(Vector nodes, Vector evaluation, String xAxisLabel, String yAxisName, String title)
        {
            StringBuilder sb = new StringBuilder(ConvertToMatLabPlotString(nodes, evaluation));
            sb.AppendLine();
            sb.Append("title('").Append(title).Append("');").AppendLine();
            sb.Append("xlabel('").Append(xAxisLabel).Append("');").AppendLine();
            sb.Append("ylabel('").Append(yAxisName).Append("');").AppendLine();

            return sb.ToString();
        }
        public static String ConvertToMatLabPlotString(Vector nodes, Vector evaluation)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(ConvertVectorToMatLabVector(nodes, "X")).AppendLine();
            sb.Append(ConvertVectorToMatLabVector(evaluation, "Y")).AppendLine();
            sb.Append("figure").AppendLine();
            sb.Append("plot(X,Y)");
            return sb.ToString();
        }

        private static String ConvertVectorToMatLabVector(Vector array, String arrayName)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(arrayName).Append(" = ");
            sb.Append("[");
            for(int i = 0; i < array.Length - 1; i++)
            {
                sb.Append(array[i] + ";");
            }
            sb.Append(array[array.Length-1] + "];");

            return sb.ToString().Replace(",",".") ;
        }

        public static String ConvertMatrixToMatLabMatrix(Matrix matrix, String matArray)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(matArray).Append(" = ");
            sb.Append("[");
            for (int rows = 0; rows < matrix.NoRows; rows++)
            {
                for (int columns = 0; columns < matrix.NoColumns; columns++)
                {
                    sb.Append( (""+matrix[rows, columns]).Replace(",",".")).Append(",");
                }
                sb.Append(";");
            }
            sb.Append("]");
            return sb.ToString();
        }

        public static String ConvertMatrixToMatlabReadable(Matrix matrix)
        {
            StringBuilder sb = new StringBuilder();
            for (int rows = 0; rows < matrix.NoRows; rows++)
            {
                for (int columns = 0; columns < matrix.NoColumns; columns++)
                {
                    sb.Append(("" + matrix[rows, columns]).Replace(",", ".")).Append(" ");
                }
                sb.AppendLine();
            }
            return sb.ToString();
        }
    }
}
