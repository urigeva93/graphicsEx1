
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class SeamCarving {
	public static double[][] computeEnergy(BufferedImage img) {
		int h = img.getHeight();
		int w = img.getWidth();
		double sumNeighbors = 0;
		double numOfNeighbors;
		Color selfColor;
		double[][] energyMat = new double[h][w];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {

				sumNeighbors = 0;
				Pixel p = new Pixel(i, j);
				List<Pixel> neighborsList = p.getNeighbors(w, h);
				numOfNeighbors = neighborsList.size();

				selfColor = new Color(img.getRGB(j, i));
				while (neighborsList.size() > 0) {

					Pixel neighbor = neighborsList.remove(0);
					Color colorNeighbor = new Color(img.getRGB(neighbor.getY(), neighbor.getX()));
					sumNeighbors += (Math.abs(colorNeighbor.getRed() - selfColor.getRed())
							+ Math.abs(colorNeighbor.getBlue() - selfColor.getBlue())
							+ Math.abs(colorNeighbor.getGreen() - selfColor.getGreen())) / 3;
				}
				energyMat[i][j] = sumNeighbors / numOfNeighbors; // normalizing
																	
			}
		}
		return energyMat;
	}

	public static double[][] computeEntropy(BufferedImage img) {
		int h = img.getHeight();
		int w = img.getWidth();

		double[][] entropyMat = new double[h][w];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				Pixel p = new Pixel(i, j);
				List<Pixel> pixelsList = p.getEnthropyMembers(w, h);
				int sumGrey = 0;

				Iterator<Pixel> it = pixelsList.iterator();

				while (it.hasNext()) {
					Pixel pGrey = it.next();
					Color colorNeighbor = new Color(img.getRGB(pGrey.getY(), pGrey.getX()));
					sumGrey += getGrayFromRGB(colorNeighbor);
				}

				it = pixelsList.iterator();
				double H = 0;
				while (it.hasNext()) {
					Pixel neighbor = it.next();
					Color colorNeighbor = new Color(img.getRGB(neighbor.getY(), neighbor.getX()));

					double P_mn = getGrayFromRGB(colorNeighbor) / sumGrey;

					H -= (P_mn * Math.log(P_mn));
				}
				entropyMat[i][j] = H;
			}
		}
		return entropyMat;
	}

	public static double getGrayFromRGB(Color c) {
		return (c.getBlue() + c.getGreen() + c.getRed()) / 3;
	}
}