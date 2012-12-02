namespace WakeCode
{
    public class SolverData
    {
        public double Ct;
        public double Dturb;
        public double Kwake;
        public double H;
        public double Uhub;
        public double Dwake;
        public double Rho;
        public double Cp;
        public double dist;

        public float[,] V = new float[1000, 1000];
        public float[,] Darea = new float[1000, 1000];
        public float[,] Darea_D = new float[1000, 1000];
    }
}