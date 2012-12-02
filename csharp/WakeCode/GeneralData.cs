using System;

namespace WakeCode
{
    public class GeneralData
    {
        public System.Int32 N_TURB;
        public System.Int32 IMAX;
        public System.Int32 JMAX;
        public double dx;
        public double dy;
        public double pi;
        public double xmax;
        public double ymax;
        public double ymin;
        public double xmin;
        public double ang;

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