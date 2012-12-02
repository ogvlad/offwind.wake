using System;
using System.IO;
using System.Reflection;

namespace WakeCode
{
    static class Program
    {
        private static void Main(string[] args)
        {
            var dir = "";
            if (args.Length > 0)
            {
                dir = args[0];
            }
            var generalData = new GeneralData();
            var calcData = new CalcData();
            var dataReader = new DataReader();
            var dataWriter = new DataWriter();
            var calc = new WakeCalc();

            dataReader.Read(generalData, dir);


            calcData.x = new double[generalData.GridPointsX];
            calcData.y = new double[generalData.GridPointsY];
            calcData.vell_i = new double[generalData.GridPointsX, generalData.GridPointsY];
            calcData.R_TURB = new double[generalData.TurbinesAmount];
            calcData.WPOWER = new double[generalData.TurbinesAmount];

            calcData.xc_turb = new Int32[generalData.TurbinesAmount];
            calcData.yc_turb = new Int32[generalData.TurbinesAmount];

            calc.Run(generalData, calcData);

            dataWriter.Write(generalData, calcData, dir);
            dataWriter.WritePower(generalData, calcData, dir);
        }
    }
}
