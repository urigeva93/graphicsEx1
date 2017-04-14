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
                for (int k = 0; k < neighborsList.size(); k++) {

                    Pixel neighbor = neighborsList.get(k);
                    Color colorNeighbor = new Color(img.getRGB(neighbor.getY(), neighbor.getX()));
                    sumNeighbors += (Math.abs(colorNeighbor.getRed() - selfColor.getRed())
                            + Math.abs(colorNeighbor.getBlue() - selfColor.getBlue())
                            + Math.abs(colorNeighbor.getGreen() - selfColor.getGreen()) / 3);
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
                List<Pixel> pixelsList = p.getEntropyMembers(w, h);
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

    public static ForwardEnergyCost[][] computeForwardEnergy(BufferedImage img) {
        int w = img.getWidth();
        int h = img.getHeight();
        ForwardEnergyCost[][] forwardMat = new ForwardEnergyCost[h][w];

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {

                forwardMat[i][j] = new ForwardEnergyCost(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
                double base;
                Color c1, c2;

                //compute base
                if (j != 0 && j != w - 1) {
                    c1 = new Color(img.getRGB(j + 1, i));
                    c2 = new Color(img.getRGB(j - 1, i));

                    base = (Math.abs(c1.getBlue() - c2.getBlue()) +
                            Math.abs(c1.getGreen() - c2.getGreen()) +
                            Math.abs(c1.getRed() - c2.getRed())) / 3;
                } else
                    base = 0;

                forwardMat[i][j].cl = base;
                forwardMat[i][j].cr = base;

                //compute C_U if possible (or 0 if not)
                forwardMat[i][j].cu = base;

                // compute C_L if possible
                if (i != 0 && j != 0) {
                    c1 = new Color(img.getRGB(j, i - 1));
                    c2 = new Color(img.getRGB(j - 1, i));

                    forwardMat[i][j].cl += (Math.abs(c1.getBlue() - c2.getBlue()) +
                            Math.abs(c1.getGreen() - c2.getGreen()) +
                            Math.abs(c1.getRed() - c2.getRed())) / 3;
                }

                // compute C_R if possible
                if (i != 0 && j != w - 1) {
                    c1 = new Color(img.getRGB(j, i - 1));
                    c2 = new Color(img.getRGB(j + 1, i));

                    forwardMat[i][j].cr += (Math.abs(c1.getBlue() - c2.getBlue()) +
                            Math.abs(c1.getGreen() - c2.getGreen()) +
                            Math.abs(c1.getRed() - c2.getRed())) / 3;
                }
            }

        }
        return forwardMat;
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
        int w = energyMat[0].length;
        int indexSeam = -1; // index for the seam col
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

    public static double[][] dynamicEnergyMatWithSeamMat(double[][] energyMat, double[][] seamMat) {
        int h = energyMat.length;
        int w = energyMat[0].length;
        for (int i = 1; i < h; i++) {
            for (int j = 0; j < w; j++) {

//                double leftNeighbor = (j == 0) ? Double.POSITIVE_INFINITY : energyMat[i - 1][j - 1];
//                double rightNeighbor = (j == w - 1) ? Double.POSITIVE_INFINITY : energyMat[i - 1][j + 1];
//
//                energyMat[i][j] = energyMat[i][j]
//                        + Math.min(leftNeighbor, Math.min(energyMat[i - 1][j], rightNeighbor));
                if (j == 0)
                    energyMat[i][j] = energyMat[i][j] + Math.min(energyMat[i - 1][j], energyMat[i - 1][j + 1]);
                else if (j == w - 1)
                    energyMat[i][j] = energyMat[i][j] + Math.min(energyMat[i - 1][j], energyMat[i - 1][j - 1]);
                else if (j > 0 && j < w - 1)
                    energyMat[i][j] = energyMat[i][j] + Math.min(energyMat[i - 1][j - 1], Math.min(energyMat[i - 1][j], energyMat[i - 1][j + 1]));

                energyMat[i][j] += seamMat[i][j];
            }
        }
        return energyMat;
    }

    public static List<Pixel> pickNextSeam(double[][] dynamicMat) {
        int h = dynamicMat.length;
        int w = dynamicMat[0].length;
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

            if (curr == 0)
                curr = (dynamicMat[i - 1][curr + 1] < dynamicMat[i - 1][curr]) ? curr + 1 : curr;
            if (curr == w - 1)
                curr = (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr]) ? curr - 1 : curr;
            if (curr > 0 && curr < w - 1) {
                int tempIndex = (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr + 1]) ? curr - 1 : curr + 1;
                curr = (dynamicMat[i - 1][tempIndex] < dynamicMat[i - 1][curr]) ? tempIndex : curr;
            }


//            if (curr == 0) {
//                if (dynamicMat[i - 1][curr] > dynamicMat[i - 1][curr + 1])
//                    curr++;
//            }
//            if (curr == w - 1) {
//                if (dynamicMat[i - 1][curr] > dynamicMat[i - 1][curr - 1]) {
//                    curr--;
//                }
//
//            }
//            if (curr > 0 && curr < w - 1) {
//                if (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr]) {
//                    if (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr + 1]) {
//                        curr--;
//                    } else
//                        curr++;
//                } else if (dynamicMat[i - 1][curr] > dynamicMat[i - 1][curr + 1])
//                    curr++;
//
//            }


//            double leftNeighbor = (curr == 0) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr - 1];
//            double rightNeighbor = (curr == w - 1) ? Double.POSITIVE_INFINITY : dynamicMat[i - 1][curr + 1];
//
//            double nextMinPixel = Math.min(leftNeighbor, Math.min(dynamicMat[i - 1][curr], rightNeighbor));
//
//            if (nextMinPixel == leftNeighbor)
//                curr--;
//            else if (nextMinPixel == rightNeighbor)
//                curr++;

            pixelSeamList.add(new Pixel(i - 1, curr));
        }
        return pixelSeamList;

    }

    public static List<List<Pixel>> pickNextKSeams(double[][] energyMat, int k) {
        List<List<Pixel>> kSeams = new ArrayList<>();
        double[][] lastEnergyMat = new double[energyMat.length][energyMat[0].length];
        for (int i = 0; i < energyMat.length; i++) {
            for (int j = 0; j < energyMat[0].length; j++) {
                lastEnergyMat[i][j] = energyMat[i][j];
            }
        }
        for (int i = 0; i < k; i++) {
            double[][] tempEnergyMat = new double[lastEnergyMat.length][lastEnergyMat[0].length - 1];
            List<Pixel> listSeam = pickNextSeam(lastEnergyMat);
            Iterator<Pixel> it = listSeam.iterator();
            for (int b = 0; b < tempEnergyMat.length; b++) {
                int outputColIndex = 0;
                Pixel p = it.next();
                for (int a = 0; a < tempEnergyMat[0].length; a++) {
                    if (p.getY() != a) {
                        tempEnergyMat[b][outputColIndex] = lastEnergyMat[b][a];
                        outputColIndex++;
                    }
                }
            }
            lastEnergyMat = tempEnergyMat;
            kSeams.add(listSeam);
        }
        return kSeams;


//        List<List<Pixel>> kSeams = new ArrayList<>();
//        double[][] seamMat = new double[energyMat.length][energyMat[0].length];
//
//        double[][] tempDynamicMat;
//
//
//        for (int i = 0; i < k; i++) {
//
//            tempDynamicMat = dynamicEnergyMatWithSeamMat(energyMat, seamMat);
//            List<Pixel> seamList = pickNextSeam(tempDynamicMat);
//            kSeams.add(seamList);
//
//            for (int j = 0; j < seamList.size(); j++) {
//                Pixel p = seamList.get(j);
//                // System.out.println("<"+p.getX()+","+p.getY()+">");
//                seamMat[p.getX()][p.getY()] = Double.POSITIVE_INFINITY;
//
//            }
//
//        }
//        return kSeams;
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

    public static BufferedImage readImage(String path) {
        try {
            BufferedImage inputImg = ImageIO.read(new File(path));
            return inputImg;
        } catch (IOException e) {
            System.out.println("Error in reading file: " + e.getMessage());
        }
        return null;
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
            pixel = seamList.get(seamList.size() - 1 - i);
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
            pixel = seamList.get(seamList.size() - 1 - i);
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

                    newImg.setRGB(pixel.getY() + t, i, newImgColor.getRGB());
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

        BufferedImage inputImg = readImage(args[0]);
        System.out.println("Height of img: " + inputImg.getHeight());
        System.out.println("Width of img: " + inputImg.getWidth());
        int desiredRows = Integer.parseInt(args[1]);
        int desiredCols = Integer.parseInt(args[2]);
        int energyTypeInt = Integer.parseInt(args[3]);
        String outputImagePath = args[4];
        EnergyTypes energyType = EnergyTypes.values()[energyTypeInt];

        int numOfIterCols = desiredCols - inputImg.getWidth();
        int numOfIterRows = desiredRows - inputImg.getHeight();
        boolean enlargeCols = (numOfIterCols > 0) ? true : false;
        boolean enlargeRows = (numOfIterRows > 0) ? true : false;

        BufferedImage outputImg = deepCopy(inputImg);

        if (Math.abs(numOfIterCols) > 0) {
            if (enlargeCols) {

                double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
                //double[][] dynamicMat = dynamicEnergyMat(energyMat);
                List<List<Pixel>> kSeams = pickNextKSeams(energyMat, Math.abs(numOfIterCols));
                System.out.println("salami");
                presentSeam(inputImg, kSeams);

//                for (int i = 0; i < kSeams.size(); i++) {
//                    List<Pixel> seamList = kSeams.get(i);
//
//                    //adjust rest seams for the new sizes of image
//                    for (int j = i + 1; j < kSeams.size(); j++) {
//                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
//                    }
//                }

                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);

                    outputImg = addSeamToImg(outputImg, seamList, false);

                    //adjust rest seams for the new sizes of image
                    for (int j = i + 1; j < kSeams.size(); j++) {
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
                    }
                }


            } else {
                for (int i = 0; i < Math.abs(numOfIterCols); i++) {
                    double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
                    double[][] dynamicMat = dynamicEnergyMat(energyMat);
                    List<Pixel> seamList = pickNextSeam(dynamicMat);
                    outputImg = removeSeamFromImg(outputImg, seamList);
                }
            }
        }


        if (Math.abs(numOfIterRows) > 0) {
            outputImg = transposeImg(outputImg); //working on traspose img for rows

            if (enlargeRows) {
                double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
                List<List<Pixel>> kSeams = pickNextKSeams(energyMat, Math.abs(numOfIterRows));


                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);
                    outputImg = addSeamToImg(outputImg, seamList, false);

                    //adjust rest seams for the new sizes of image
                    for (int j = i + 1; j < kSeams.size(); j++)
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
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
        }

        saveImage(outputImagePath, outputImg);
    }

    public static List<Pixel> adjustSeamForAdd(List<Pixel> currPixelList, List<Pixel> seamList) {
        for (int k = 0; k < currPixelList.size(); k++) {
            Pixel p = currPixelList.get(k);
            p.addToCol(seamList.get(k));
        }
        return currPixelList;
    }

    public static void presentSeam(BufferedImage img, List<List<Pixel>> kSeams) {
        BufferedImage outputImg = deepCopy(img);
        Color c = new Color(255, 0, 0);

        for (int i = 0; i < kSeams.size(); i++) {
            List<Pixel> listSeam = kSeams.get(i);
            for (int j = 0; j < listSeam.size(); j++) {
                Pixel p = listSeam.get(j);
                outputImg.setRGB(p.getY(), p.getX(), c.getRGB());
            }
        }

        saveImage("//Users//uri//Documents//workspace//picEx1//src//output4.jpg", outputImg);

    }

    public enum EnergyTypes {
        REGULAR, ENTROPY, FORWARD
    }
}