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

            dataReader.READ_DATA(solverData, generalData);
            calc.Run(generalData, solverData);


            dataWriter.WRITE_DATA(solverData, generalData);
            dataWriter.WRITE_DATA_power(generalData);
        }
    }
}
