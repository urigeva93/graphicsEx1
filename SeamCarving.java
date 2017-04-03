
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

	public static double[][] transposeMat(double[][] mat) {
		int h = mat.length;
		int w = mat[1].length;
		double[][] res = new double[w][h];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				res[j][i] = mat[i][j];
			}
		}
		return res;
	}

	public static int straightSeam(double[][] energyMat) {
		int h = energyMat.length;
		int w = energyMat[1].length;
		int indexSeam = -1; // index for the seam colum
		double minSum = Double.POSITIVE_INFINITY, tempSum = 0;
		for (int j = 0; j < w; j++) {
			for (int i = 0; i < h; i++) {
				tempSum += energyMat[i][j];
			}
			if (tempSum < minSum) {
				minSum = tempSum;
				indexSeam = j;
			}
			tempSum = 0;
		}
		return indexSeam;
	}

	public static double[][] dynamicEnergyMat(double[][] energyMat) {
		int h = energyMat.length;
		int w = energyMat[1].length;
		for (int i = 1; i < h; i++) {
			for (int j = 0; j < w; j++) {

				double leftNeighbor = (j == 0) ? Double.POSITIVE_INFINITY : energyMat[i - 1][j - 1];
				double rightNeighbor = (j == w - 1) ? Double.POSITIVE_INFINITY : energyMat[i - 1][j + 1];

				energyMat[i][j] = energyMat[i][j]
						+ Math.min(leftNeighbor, Math.min(energyMat[i - 1][j], rightNeighbor));
			}
		}
		return energyMat;

	}

	public static List<Pixel> pickPixelsSeam(double[][] dynamicMat) {
		int h = dynamicMat.length;
		int w = dynamicMat[1].length;
		int minIndex = -1;
		double minLastRow = Double.POSITIVE_INFINITY; 
		List<Pixel> pixelsList = new ArrayList<Pixel>();

		// finding the minimum in the last line
		for (int j = 0; j < w; j++) {
			if (dynamicMat[h - 1][j] < minLastRow) {
				minLastRow = dynamicMat[h - 1][j];
				minIndex = j;
			}
		}
		pixelsList.add(new Pixel(h - 1, minIndex));
		int curr = minIndex;
		for (int i = h - 1; i > 0; i--) {
			double leftNeighbor = (curr == 0) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr - 1];
			double rightNeighbor = (curr == w - 1) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr + 1];

			double nextMinPixel = Math.min(leftNeighbor, Math.min(dynamicMat[i - 1][curr], rightNeighbor));

			if (nextMinPixel == dynamicMat[i - 1][curr - 1])
				curr--;
			else if (nextMinPixel == dynamicMat[i - 1][curr + 1])
				curr++;

			pixelsList.add(new Pixel(i - 1, curr));
		}
		return pixelsList;

	}
}
