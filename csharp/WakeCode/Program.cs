using System.IO;
using System.Reflection;

namespace WakeCode
{
    static class Program
    {
        private static void Main(string[] args)
        {
            var dir = ".\\";
            if (args.Length > 0)
            {
                dir = args[0];
                if (!Directory.Exists(dir)) Directory.CreateDirectory(dir);
            }
            var solverData = new SolverData();
            var generalData = new GeneralData();
            var dataReader = new DataReader();
            var dataWriter = new DataWriter();
            var calc = new WakeCalc();

            dataReader.Read(solverData, generalData, dir);
            calc.Run(generalData, solverData);

            dataWriter.Write(solverData, generalData, dir);
            dataWriter.WritePower(generalData, dir);
        }
    }
}
