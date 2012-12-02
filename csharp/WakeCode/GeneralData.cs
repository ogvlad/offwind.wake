using System;

namespace WakeCode
{
    public class GeneralData
    {
        public System.Int32 TurbinesAmount;
        public System.Int32 GridPointsX;
        public System.Int32 GridPointsY;
        public double dx;
        public double dy;
        public double pi;
        public double xmax;
        public double ymax;
        public double ymin;
        public double xmin;
        public double RotationAngle;

        public double[] x;
        public double[] y;
        public double[,] vell_i;
        public double[] x_turb;     // location of the turbine
        public double[] y_turb;     // location of the turbine
        public double[] R_TURB;     // location of the turbine
        public double[] WPOWER;     // location of the turbine
        public Int32[] xc_turb;
        public Int32[] yc_turb;
    }
}