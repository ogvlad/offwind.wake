namespace WakeCode
{
    static class Program
    {
        private static void Main(string[] args)
        {
            var solverData = new SolverData();
            var generalData = new GeneralData();
            var dataReader = new DataReader();
            var dataWriter = new DataWriter();
            var calc = new WakeCalc();

            dataReader.Read(solverData, generalData);
            calc.Run(generalData, solverData);

            dataWriter.Write(solverData, generalData);
            dataWriter.WritePower(generalData);
        }
    }
}
