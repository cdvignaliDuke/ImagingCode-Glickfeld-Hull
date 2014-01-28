import java.io.*;
import java.awt.image.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.process.*;
import ij.plugin.PlugIn;
import ij.measure.*;

/**
		Does various reconstruction operations on a hypervolume
		(a series of 3d volumes) in a ImageJ stack.
				choice variable:
						"shuffle  (xytz -> xyzt)", 	// 0
						"unshuffle (xyzt -> xytz)",	// 1
*/
public class Hypervolume_ShufflerAK implements PlugIn
{
		public int 			choice = 0;
		protected Calibration	c;
		public int			depth = 1;
        
        public Hypervolume_ShufflerAK() {;}
        
        public void run(String arg) {;}
		public ImagePlus shuffle(ImagePlus imp, int choice, int depthin)
		{
			ImagePlus impn = imp; 
					depth = depthin;
                    ImageStack is = imp.getStack();
					c = imp.getCalibration();
					switch (choice)
					{
						case 0: {impn=shuffle(imp); return impn;}
						case 1: {impn=shufflerev(imp); return impn;}
                        default: {return impn;}
					}
		}
		private ImagePlus shuffle(ImagePlus imp)
		/*
				Change the order from xytz to xyzt.
		*/
		{
						ImageStack is = imp.getStack();
						// calculate number of volumes in hypervolume.
						int length = is.getSize() / depth;
						ImageStack isn = new ImageStack(is.getWidth(), is.getHeight());
						for (int t = 0; t < length; t++)
								for (int z = 0; z < depth; z++)
										isn.addSlice(""+ (z * length + t + 1),
												is.getProcessor(z * length + t + 1));
						ImagePlus impn = new ImagePlus("Shuffled "+imp.getTitle(), isn);
						impn.setCalibration(c);
						return impn;
		}
		private ImagePlus shufflerev(ImagePlus imp)
		/*
				Change the order from xyzt to xytz.
		*/
		{
						ImageStack is = imp.getStack();
						// calculate number of volumes in hypervolume.
						int length = is.getSize() / depth;
						ImageStack isn = new ImageStack(is.getWidth(), is.getHeight());
						for (int z = 0; z < depth; z++)
								for (int t = 0; t < length; t++)
										isn.addSlice(""+ (t * depth + z + 1),
												is.getProcessor(t * depth + z + 1));
						ImagePlus impn = new ImagePlus("Unshuffled "+imp.getTitle(), isn);
						impn.setCalibration(c);
						return impn;
		}

}


