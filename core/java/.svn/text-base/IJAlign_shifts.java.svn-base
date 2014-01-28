import java.io.*;
import java.awt.image.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.process.*;
import ij.plugin.PlugIn;
import ij.measure.*;

public class IJAlign_shifts implements PlugIn
{
        
        public IJAlign_shifts() {;}
        
        public double [][] doAlign(String cmdstr, ImagePlus source, ImagePlus target) { 
            TurboReg_AK tr = new TurboReg_AK();
            int nSlices = source.getNSlices();
            double[][] shifts = new double [nSlices][2];
            double[][] sourcePoints = new double [4][2];
            double[][] targetPoints = new double [4][2];
            for (int i = 0; i < nSlices; i++) {
                source.setSlice((i+1));
                ImagePlus sourceSlice = new ImagePlus("",source.getProcessor());
                tr.run(cmdstr, sourceSlice, target);
                sourcePoints = tr.getSourcePoints();
                targetPoints = tr.getTargetPoints();
                shifts [i][0] = (sourcePoints [0][0] - targetPoints [0][0]);
                shifts [i][1] = (sourcePoints [0][1] - targetPoints [0][1]);
			}
            source.setSlice(1);
            return shifts;
        }
             
        public void run(String arg) {;}
}