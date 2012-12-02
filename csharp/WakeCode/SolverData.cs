namespace WakeCode
{
    public class SolverData
    {
        public double TurbineThrust;
        public double TurbineDiameter;
        public double WakeDecay;
        public double TurbineHeight;
        public double VelocityAtHub;
        public double Dwake;
        public double AirDensity;
        public double Cp;
        public double PowerDistance;

        public float[,] V = new float[1000, 1000];
        public float[,] Darea = new float[1000, 1000];
        public float[,] Darea_D = new float[1000, 1000];
    }
}