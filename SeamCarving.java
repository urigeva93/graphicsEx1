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
                    sumNeighbors = sumNeighbors + (Math.abs(colorNeighbor.getRed() - selfColor.getRed())
                            + Math.abs(colorNeighbor.getBlue() - selfColor.getBlue())
                            + Math.abs(colorNeighbor.getGreen() - selfColor.getGreen()));
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

                    H -= (P_mn * (P_mn != 0 ? Math.log(P_mn) : 1));
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

    public static int pickNextStraightSeam(double[][] dynamicMat) {
        int h = dynamicMat.length;
        int w = dynamicMat[0].length;
        int indexSeam = 0; // index for the seam col
        double minSum = dynamicMat[h - 1][0], tempSum = 0;
        for (int j = 0; j < w; j++) {
            for (int i = 0; i < h; i++)
                tempSum += dynamicMat[i][j];
            if (tempSum < minSum) {
                minSum = tempSum;
                indexSeam = j;
            }
            tempSum = 0;
        }
        return indexSeam;
    }

    public static double[][] controlDynamicMat(BufferedImage img, EnergyTypes type) {

        double[][] energyMat = computeEnergyFromImage(img, type);

        if (type == EnergyTypes.REGULAR || type == EnergyTypes.ENTROPY) {
            return dynamicEnergyMat(energyMat);
        } else { // type == 'FORWARD'
            ForwardEnergyCost[][] forwardMat = computeForwardEnergy(img);
            return dynamicForwardMat(energyMat, forwardMat);
        }
    }

    public static double[][] dynamicEnergyMat(double[][] energyMat) {
        int h = energyMat.length;
        int w = energyMat[0].length;
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

    public static double[][] dynamicForwardMat(double[][] energyMat, ForwardEnergyCost[][] forwardCost) {
        int h = energyMat.length;
        int w = energyMat[0].length;
        double cu, cl, cr;

        double[][] dynamicForward = new double[h][w];

        //compute first row (don't have father values)
        for (int j = 0; j < w; j++) {
            cu = forwardCost[0][j].getCu();
            cl = forwardCost[0][j].getCl();
            cr = forwardCost[0][j].getCr();

            dynamicForward[0][j] = energyMat[0][j];

            if (j == 0)
                dynamicForward[0][j] += Math.min(cu, cr);
            else if (j == w - 1)
                dynamicForward[0][j] += Math.min(cu, cl);
            else
                dynamicForward[0][j] += Math.min(cu, Math.min(cl, cr));

        }

        for (int i = 1; i < h; i++) {
            for (int j = 0; j < w; j++) {

                cu = forwardCost[i][j].getCu();
                cl = forwardCost[i][j].getCl();
                cr = forwardCost[i][j].getCr();

                dynamicForward[i][j] = energyMat[i][j];

                if (j == 0)
                    dynamicForward[i][j] += Math.min(dynamicForward[i - 1][j] + cu, dynamicForward[i - 1][j + 1] + cr);
                else if (j == w - 1)
                    dynamicForward[i][j] += Math.min(dynamicForward[i - 1][j] + cu, dynamicForward[i - 1][j - 1] + cl);
                else
                    dynamicForward[i][j] += Math.min(dynamicForward[i - 1][j] + cu, Math.min(dynamicForward[i - 1][j - 1] + cl,
                            dynamicForward[i - 1][j + 1] + cr));

            }
        }
        return dynamicForward;

    }

    public static List<Pixel> pickNextSeam(double[][] dynamicMat) {
        int h = dynamicMat.length;
        int w = dynamicMat[0].length;
        int minIndex = 0;
        double minLastRow = dynamicMat[h - 1][0];
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
            else if (curr == w - 1)
                curr = (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr]) ? curr - 1 : curr;
            else {
                int tempIndex = (dynamicMat[i - 1][curr - 1] < dynamicMat[i - 1][curr + 1]) ? curr - 1 : curr + 1;
                curr = (dynamicMat[i - 1][tempIndex] < dynamicMat[i - 1][curr]) ? tempIndex : curr;
            }

            pixelSeamList.add(new Pixel(i - 1, curr));
        }
        return pixelSeamList;

    }

    public static List<List<Pixel>> pickNextKSeams(BufferedImage img, EnergyTypes type, int k) {


        double[][] dynamicMat;
        BufferedImage workImg = deepCopy(img);
        List<List<Pixel>> seamsList = new ArrayList<>();
        for (int i = 0; i < k; i++) {
//            energyMat = computeEnergyFromImage(workImg, type);
//            dynamicMat = dynamicEnergyMat(energyMat);
            dynamicMat = controlDynamicMat(workImg, type);
            List<Pixel> seam = pickNextSeam(dynamicMat);
            seamsList.add(seam);
            workImg = removeSeamFromImg(workImg, seam);

        }
        return seamsList;
    }

    public static double[][] computeEnergyFromImage(BufferedImage img, EnergyTypes type) {
        double alpha = 0.25;
        int w = img.getWidth();
        int h = img.getHeight();
        double[][] compMat = new double[h][w];

        if (type == EnergyTypes.REGULAR || type == EnergyTypes.FORWARD)
            compMat = computeEnergy(img);
        if (type == EnergyTypes.ENTROPY) {
            double[][] entropyMat = computeEntropy(img);
            double[][] energyMat = computeEnergy(img);
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++)
                    compMat[i][j] = alpha * energyMat[i][j] + (1 - alpha) * entropyMat[i][j];
            }
        }
        return compMat;
    }

    public static void saveImage(String path, BufferedImage inputImg) {
        try {
            ImageIO.write(inputImg, "jpg", new File(path));
        } catch (IOException e) {
            System.out.println("Error in saving image: " + e.getMessage());
        }
    }

    public static BufferedImage readImage(String path) {
        try {
            BufferedImage inputImg = ImageIO.read(new File(path));
            return inputImg;
        } catch (IOException e) {
            System.out.println("Error in reading image: " + e.getMessage());
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
            if (enlargeCols) { //add seams

                List<List<Pixel>> kSeams = pickNextKSeams(outputImg, energyType, Math.abs(numOfIterCols));
                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);

                    //adjust rest seams for the original size
                    for (int j = i + 1; j < kSeams.size(); j++) {
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
                    }
                }

                presentSeam(inputImg, kSeams);

                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);

                    outputImg = addSeamToImg(outputImg, seamList, false);
                    //adjust rest seams for the new sizes of image
                    for (int j = i + 1; j < kSeams.size(); j++) {
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
                    }
                }


            } else { //remove seams
                for (int i = 0; i < Math.abs(numOfIterCols); i++) {
//                    double[][] energyMat = computeEnergyFromImage(outputImg, energyType);
//                    double[][] dynamicMat = dynamicEnergyMat(energyMat);
                    double[][] dynamicMat = controlDynamicMat(outputImg, energyType);
                    List<Pixel> seamList = pickNextSeam(dynamicMat);
                    outputImg = removeSeamFromImg(outputImg, seamList);
                }
            }
        }


        if (Math.abs(numOfIterRows) > 0) {
            outputImg = transposeImg(outputImg); //working on transposed img for rows

            if (enlargeRows) {
                List<List<Pixel>> kSeams = pickNextKSeams(outputImg, energyType, Math.abs(numOfIterRows));

                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);

                    //adjust rest seams for the original size
                    for (int j = i + 1; j < kSeams.size(); j++) {
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
                    }
                }

                for (int i = 0; i < kSeams.size(); i++) {
                    List<Pixel> seamList = kSeams.get(i);

                    outputImg = addSeamToImg(outputImg, seamList, false);

                    //adjust rest seams for the new sizes of image
                    for (int j = i + 1; j < kSeams.size(); j++) {
                        kSeams.set(j, adjustSeamForAdd(kSeams.get(j), seamList));
                    }
                }

            } else {
                for (int i = 0; i < Math.abs(numOfIterRows); i++) {
                    double[][] dynamicMat = controlDynamicMat(outputImg, energyType);
                    List<Pixel> seamList = pickNextSeam(dynamicMat);
                    outputImg = removeSeamFromImg(outputImg, seamList);
                }
            }

            outputImg = transposeImg(outputImg); //transposed again to get original
        }
        //save the output image
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