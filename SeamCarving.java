import javax.imageio.ImageIO;
import java.awt.*;
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
        double sumNeighbors;
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

    public static BufferedImage transposeImg(BufferedImage img) {
        int h = img.getHeight();
        int w = img.getWidth();
        BufferedImage transImg = new BufferedImage(h, w, img.getType());
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                transImg.setRGB(i, j, img.getRGB(j, i));
            }
        }
        return transImg;
    }

    public static int straightSeam(double[][] energyMat) {
        int h = energyMat.length;
        int w = energyMat[1].length;
        int indexSeam = -1; // index for the seam colum
        double minSum = Double.POSITIVE_INFINITY, tempSum = 0;
        for (int j = 0; j < w; j++) {
            for (int i = 0; i < h; i++)
                tempSum += energyMat[i][j];
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

    public static List<Pixel> pickNextSeam(double[][] dynamicMat) {
        int h = dynamicMat.length;
        int w = dynamicMat[1].length;
        int minIndex = -1;
        double minLastRow = Double.POSITIVE_INFINITY;
        List<Pixel> pixelSeamList = new ArrayList<>();

        // find the minimum in the last line
        for (int j = 0; j < w; j++) {
            if (dynamicMat[h - 1][j] < minLastRow) {
                minLastRow = dynamicMat[h - 1][j];
                minIndex = j;
            }
        }
        pixelSeamList.add(new Pixel(h - 1, minIndex));
        int curr = minIndex;
        for (int i = h - 1; i > 0; i--) {
            double leftNeighbor = (curr == 0) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr - 1];
            double rightNeighbor = (curr == w - 1) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr + 1];

            double nextMinPixel = Math.min(leftNeighbor, Math.min(dynamicMat[i - 1][curr], rightNeighbor));

            if (nextMinPixel == leftNeighbor)
                curr--;
            else if (nextMinPixel == rightNeighbor)
                curr++;

            pixelSeamList.add(new Pixel(i - 1, curr));
        }
        return pixelSeamList;

    }


    //TODO kSeams not good
    public static List<List<Pixel>> pickNextKSeams(double[][] dynamicMat, int k) {
        List<List<Pixel>> kSeams = new ArrayList<>();
        double[][] lastDynamicMat = new double[dynamicMat.length][dynamicMat[0].length];
        //copy dynamicMat
        for (int i = 0; i < dynamicMat.length; i++) {
            for (int j = 0; j < dynamicMat[0].length; j++) {
                lastDynamicMat[i][j] = dynamicMat[i][j];
            }
        }
        for (int i = 0; i < k; i++) {
            double[][] tempDynamicMat = new double[lastDynamicMat.length][lastDynamicMat[0].length - 1];
            List<Pixel> seamList = pickNextSeam(lastDynamicMat);
            kSeams.add(seamList);
            for (int p = 0; p < lastDynamicMat.length; p++) {
                Pixel pix = seamList.remove(seamList.size() - 1);
                for (int j = 0; j < pix.getY(); j++) {
                    tempDynamicMat[p][j] = lastDynamicMat[p][j];
                }
                for (int j = pix.getY(); j < tempDynamicMat[0].length; j++) {
                    tempDynamicMat[p][j] = lastDynamicMat[p][j + 1];
                }
            }
            lastDynamicMat = tempDynamicMat;
        }
        return kSeams;
    }

    public static double[][] computeEnergyFromImage(BufferedImage img, EnergyTypes type) {
        int w = img.getWidth();
        int h = img.getHeight();
        double[][] combinedMat = new double[h][w];
        double[][] energyMat = computeEnergy(img);

        if (type == EnergyTypes.REGULAR)
            return energyMat;
        else {
            if (type == EnergyTypes.ENTROPY) {
                double[][] entropyMat = computeEntropy(img);
                for (int i = 0; i < h; i++) {
                    for (int j = 0; j < w; j++) {
                        combinedMat[i][j] = energyMat[i][j] + 3 * entropyMat[i][j];
                    }
                }
            } else {
                // return forward;
            }
        }
        return combinedMat;

    }

    public static void saveImage(String path, BufferedImage inputImg) {
        try {
            ImageIO.write(inputImg, "jpg", new File(path));
        } catch (IOException e) {
            System.out.println("Error in saving file: " + e.getMessage());
        }
    }

    public static BufferedImage deepCopy(BufferedImage img) {
        ColorModel cm = img.getColorModel();
        boolean isAlpha = cm.isAlphaPremultiplied();
        WritableRaster raster = img.copyData(null);
        return new BufferedImage(cm, raster, isAlpha, null);

    }

    public static BufferedImage removeSeamFromImg(BufferedImage oldImg, List<Pixel> seamList) {
        int h = oldImg.getHeight();
        int w = oldImg.getWidth();
        Pixel pixel;
        BufferedImage newImg = new BufferedImage(w - 1, h, oldImg.getType());
        for (int i = 0; i < h; i++) {
            pixel = seamList.remove(seamList.size() - 1);
            for (int j = 0; j < pixel.getY(); j++) { // first half of row -
                // until removed pixel
                newImg.setRGB(j, i, oldImg.getRGB(j, i));
            }
            for (int k = pixel.getY(); k < w - 1; k++) { // second half of row -
                // from removed
                // pixel
                newImg.setRGB(k, i, oldImg.getRGB(k + 1, i));
            }
        }
        return newImg;
    }

    public static BufferedImage addSeamToImg(BufferedImage oldImg, List<Pixel> seamList, boolean interpolation) {
        int h = oldImg.getHeight();
        int w = oldImg.getWidth();
        Pixel pixel;
        BufferedImage newImg = new BufferedImage(w + 1, h, oldImg.getType());
        for (int i = 0; i < h; i++) {
            pixel = seamList.remove(seamList.size() - 1);
            for (int j = 0; j < pixel.getY(); j++) { // first half of row -
                // until removed pixel
                newImg.setRGB(j, i, oldImg.getRGB(j, i));
            }

            // two added pixel
            for (int t = 0; t <= 1; t++) {
                if (interpolation && pixel.getY() != 0 && pixel.getY() != w - 1) {
                    Color leftColor = new Color(oldImg.getRGB(pixel.getY() - 1, i));
                    Color rightColor = new Color(oldImg.getRGB(pixel.getY() + 1, i));
                    Color imgColor = new Color(oldImg.getRGB(pixel.getY(), i));

                    int red = ((leftColor.getRed() + rightColor.getRed() + imgColor.getRed()) / 3);
                    int green = ((leftColor.getGreen() + rightColor.getGreen() + imgColor.getGreen()) / 3);
                    int blue = ((leftColor.getBlue() + rightColor.getBlue() + imgColor.getBlue()) / 3);
                    Color newImgColor = new Color(red, green, blue);

                    newImg.setRGB(pixel.getY()+ t, i, newImgColor.getRGB());
                } else
                    newImg.setRGB(pixel.getY() + t, i, oldImg.getRGB(pixel.getY(), i));
            }

            for (int k = pixel.getY() + 1; k < w; k++) { // second half of row -
                // from removed
                // pixel
                newImg.setRGB(k + 1, i, oldImg.getRGB(k, i));
            }
        }
        return newImg;
    }

    public static void main(String[] args) throws IOException {
        BufferedImage inputImg = ImageIO.read(new File(args[0]));
        System.out.println("Height of img: "+inputImg.getHeight());
        System.out.println("Width of img: "+inputImg.getWidth());

        int desiredRows = Integer.parseInt(args[1]);
        int desiredCols = Integer.parseInt(args[2]);
        int energyTypeInt = Integer.parseInt(args[3]);
        String outputImagePath = args[4];
        EnergyTypes energyType = EnergyTypes.values()[energyTypeInt];

        int numOfIterCols = desiredCols - inputImg.getWidth();
        int numOfIterRows = desiredRows - inputImg.getHeight();

        boolean enlargeCols = (numOfIterCols >= 0) ? true : false;
        boolean enlargeRows = (numOfIterRows >= 0) ? true : false;

        BufferedImage outputImg = deepCopy(inputImg);

        if (enlargeCols) {
//			for(int i = 0; i < Math.abs(numOfIterCols); i++){
//				double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
//				double[][] dynamicMat = dynamicEnergyMat(energyMat);
//				List<Pixel> seamList = pickNextSeam(dynamicMat);
//				outputImg = addSeamToImg(outputImg, seamList, true);
//			}
            double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
            double[][] dynamicMat = dynamicEnergyMat(energyMat);
            List<List<Pixel>> kSeams = pickNextKSeams(dynamicMat, Math.abs(numOfIterCols));
            while (kSeams.size() > 0) {
                List<Pixel> seamList = kSeams.remove(0);
                outputImg = addSeamToImg(outputImg, seamList, false);
            }
        } else {
            for (int i = 0; i < Math.abs(numOfIterCols); i++) {
                double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
                double[][] dynamicMat = dynamicEnergyMat(energyMat);
                List<Pixel> seamList = pickNextSeam(dynamicMat);
                outputImg = removeSeamFromImg(outputImg, seamList);
            }
        }

        outputImg = transposeImg(outputImg); //working on traspose img for rows

        if (enlargeRows) {
            double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
            double[][] dynamicMat = dynamicEnergyMat(energyMat);
            List<List<Pixel>> kSeams = pickNextKSeams(dynamicMat, Math.abs(numOfIterRows));
            while (kSeams.size() > 0) {
                List<Pixel> seamList = kSeams.remove(0);
                outputImg = addSeamToImg(outputImg, seamList, false);
            }
        } else {
            for (int i = 0; i < Math.abs(numOfIterRows); i++) {
                double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
                double[][] dynamicMat = dynamicEnergyMat(energyMat);
                List<Pixel> seamList = pickNextSeam(dynamicMat);
                outputImg = removeSeamFromImg(outputImg, seamList);
            }
        }

        outputImg = transposeImg(outputImg); //trasnposed again to get original
        saveImage(outputImagePath, outputImg);
    }

    public enum EnergyTypes {
        REGULAR, ENTROPY, FORWARD
    }


}