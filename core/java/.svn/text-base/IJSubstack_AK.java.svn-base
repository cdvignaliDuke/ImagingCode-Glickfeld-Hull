import java.io.*;
import java.awt.image.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.process.*;
import ij.plugin.PlugIn;
import ij.measure.*;

public class IJSubstack_AK implements PlugIn
{
        public IJSubstack_AK() {;}
        
        public ImagePlus doSubstack(ImagePlus imp, int start, int end) {
            int width = imp.getWidth();
            int height = imp.getHeight();
            ImageStack stack = imp.getStack();
            ImageStack stack2 = new ImageStack(width, height, imp.getProcessor().getColorModel());
            for (int i=start; i<=end; i++) {
                ImageProcessor ip2 = stack.getProcessor(i);
                ip2 = ip2.crop();
                stack2.addSlice(stack.getSliceLabel(i), ip2);
            }
            
            ImagePlus impSubstack = imp.createImagePlus();
            impSubstack.setStack("", stack2);
            impSubstack.setCalibration(imp.getCalibration());
            return impSubstack;
        }
        
        public void run(String arg) {;}
}        