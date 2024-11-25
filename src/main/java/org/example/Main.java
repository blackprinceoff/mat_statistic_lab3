package org.example;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.Arrays;
import java.util.Random;

public class Main {
    private static final int NUM_BINS = 10;
    private static final int SAMPLE_SIZE = 350;
    private static final double CONFIDENCE_LEVEL = 0.95;

    public static void main(String[] args) {
        double[] data = generateRandomSample();

        System.out.println("Згенеровані значення:");
        for (double value : data) {
            System.out.printf("%.4f ", value);
        }
        System.out.println("\n");

        double mean = mean(data);
        double stddev = standardDeviation(data);
        System.out.printf("Математичне сподівання: %.4f%n", mean);
        System.out.printf("Стандартне відхилення: %.4f%n", stddev);

        createCumulativeHistogram(data, NUM_BINS);
        createHistogramWithNormalApproximation(data);
        calculateConfidenceIntervals(data);
        performChiSquareTest(data);
    }


    private static void createCumulativeHistogram(double[] data, int numBins) {
        double min = Arrays.stream(data).min().orElse(0.0);
        double max = Arrays.stream(data).max().orElse(0.0);
        double binWidth = (max - min) / numBins;

        double[] cumulativeFrequencies = calculateCumulativeRelativeFrequencies(data, numBins);

        XYIntervalSeries seriesCumulativeFrequencies = new XYIntervalSeries("Накопичені відносні частоти");
        for (int i = 0; i < numBins; i++) {
            double binStart = min + i * binWidth;
            double binEnd = binStart + binWidth;
            double binCenter = (binStart + binEnd) / 2;
            seriesCumulativeFrequencies.add(binCenter, binStart, binEnd, cumulativeFrequencies[i], cumulativeFrequencies[i], cumulativeFrequencies[i]);
        }

        XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
        dataset.addSeries(seriesCumulativeFrequencies);

        JFreeChart chart = ChartFactory.createXYBarChart(
                "Гістограма накопичених відносних частот",
                "Значення",
                false,
                "Накопичена частота",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = chart.getXYPlot();
        XYBarRenderer renderer = new XYBarRenderer();
        renderer.setMargin(0.0); // Set margin to zero to connect the rectangles
        renderer.setDrawBarOutline(false);

        // Add tooltips
        renderer.setDefaultToolTipGenerator((dataset1, series, item) -> {
            double binStart = min + item * binWidth;
            double binEnd = binStart + binWidth;
            double frequency = cumulativeFrequencies[item];
            return String.format("Інтервал: [%.2f, %.2f], Частота: %.4f", binStart, binEnd, frequency);
        });

        plot.setRenderer(renderer);

        JFrame frame = new JFrame("Гістограма накопичених відносних частот");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }



    private static void createHistogramWithNormalApproximation(double[] data) {
        double min = Arrays.stream(data).min().orElse(0.0);
        double max = Arrays.stream(data).max().orElse(0.0);
        double binWidth = (max - min) / NUM_BINS;

        double[] cumulativeFrequencies = calculateCumulativeRelativeFrequencies(data, NUM_BINS);

        // Create dataset for the histogram
        XYIntervalSeries seriesCumulativeFrequencies = new XYIntervalSeries("Накопичені відносні частоти");
        for (int i = 0; i < NUM_BINS; i++) {
            double binStart = min + i * binWidth;
            double binEnd = binStart + binWidth;
            double binCenter = (binStart + binEnd) / 2;
            seriesCumulativeFrequencies.add(binCenter, binStart, binEnd, cumulativeFrequencies[i], cumulativeFrequencies[i], cumulativeFrequencies[i]);
        }

        XYIntervalSeriesCollection histogramDataset = new XYIntervalSeriesCollection();
        histogramDataset.addSeries(seriesCumulativeFrequencies);

        // Create dataset for the normal approximation line
        XYSeries seriesNormalApproximation = new XYSeries("Апроксимація нормальним розподілом");
        NormalDistribution normalDistribution = new NormalDistribution(mean(data), standardDeviation(data));
        for (double x = min; x <= max; x += (max - min) / 200) {
            seriesNormalApproximation.add(x, normalDistribution.cumulativeProbability(x));
        }

        XYSeriesCollection lineDataset = new XYSeriesCollection();
        lineDataset.addSeries(seriesNormalApproximation);

        // Create the chart
        JFreeChart chart = ChartFactory.createXYBarChart(
                "Гістограма та апроксимація нормальним розподілом",
                "Значення",
                false,
                "Накопичена частота",
                histogramDataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = chart.getXYPlot();

        // Configure the histogram renderer
        XYBarRenderer barRenderer = new XYBarRenderer();
        barRenderer.setMargin(0.0); // Remove margin to connect bars
        barRenderer.setBarAlignmentFactor(1.0); // Fully fill the space
        barRenderer.setDrawBarOutline(false); // Remove bar outlines
        plot.setDataset(0, histogramDataset);
        plot.setRenderer(0, barRenderer);

        // Configure the line renderer
        XYLineAndShapeRenderer lineRenderer = new XYLineAndShapeRenderer(true, false);
        plot.setDataset(1, lineDataset);
        plot.setRenderer(1, lineRenderer);

        // Ensure the line is drawn above the bars
        plot.setDatasetRenderingOrder(org.jfree.chart.plot.DatasetRenderingOrder.FORWARD);

        // Show the chart in a frame
        JFrame frame = new JFrame("Гістограма з апроксимацією нормальним розподілом");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }


    private static double[] generateRandomSample() {
        int groupNumber = 15;
        double mean = groupNumber - 10; // Математичне сподівання
        double stddev = 3 + (double) groupNumber / 10; // Середнє-квадратичне відхилення

        Random random = new Random();
        double[] sample = new double[SAMPLE_SIZE];
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            sample[i] = mean + stddev * random.nextGaussian();
        }
        return sample;
    }


    private static double[] calculateCumulativeRelativeFrequencies(double[] data, int numBins) {
        double min = Arrays.stream(data).min().orElse(0.0);
        double max = Arrays.stream(data).max().orElse(0.0);
        double binWidth = (max - min) / numBins;
        int[] binCounts = new int[numBins];
        double[] cumulativeFrequencies = new double[numBins];

        for (double value : data) {
            int binIndex = (int) ((value - min) / binWidth);
            if (binIndex >= numBins) binIndex = numBins - 1;
            binCounts[binIndex]++;
        }

        int cumulativeCount = 0;
        for (int i = 0; i < numBins; i++) {
            cumulativeCount += binCounts[i];
            cumulativeFrequencies[i] = (double) cumulativeCount / SAMPLE_SIZE;
        }

        return cumulativeFrequencies;
    }

    private static double mean(double[] data) {
        return Arrays.stream(data).average().orElse(0.0);
    }

    private static double standardDeviation(double[] data) {
        double mean = mean(data);
        return Math.sqrt(Arrays.stream(data).map(x -> Math.pow(x - mean, 2)).average().orElse(0.0));
    }

    private static void calculateConfidenceIntervals(double[] data) {
        double mean = mean(data);
        double stddev = standardDeviation(data);
        int n = data.length;

        TDistribution tDist = new TDistribution(n - 1);
        double tCritical = tDist.inverseCumulativeProbability(1 - (1 - CONFIDENCE_LEVEL) / 2);
        double marginOfError = tCritical * stddev / Math.sqrt(n);
        double meanLowerBound = mean - marginOfError;
        double meanUpperBound = mean + marginOfError;
        System.out.println("Довірчий інтервал для математичного сподівання: [" + meanLowerBound + ", " + meanUpperBound + "]");

        ChiSquaredDistribution chi2DistLower = new ChiSquaredDistribution(n - 1);
        ChiSquaredDistribution chi2DistUpper = new ChiSquaredDistribution(n - 1);
        double chi2CriticalLower = chi2DistLower.inverseCumulativeProbability((1 - CONFIDENCE_LEVEL) / 2);
        double chi2CriticalUpper = chi2DistUpper.inverseCumulativeProbability(1 - (1 - CONFIDENCE_LEVEL) / 2);
        double variance = stddev * stddev;
        double varianceLowerBound = (n - 1) * variance / chi2CriticalUpper;
        double varianceUpperBound = (n - 1) * variance / chi2CriticalLower;
        System.out.println("Довірчий інтервал для дисперсії: [" + varianceLowerBound + ", " + varianceUpperBound + "]");
    }

    private static void performChiSquareTest(double[] data) {
        int numBins = NUM_BINS;

        double min = Arrays.stream(data).min().orElse(0.0);
        double max = Arrays.stream(data).max().orElse(0.0);
        double binWidth = (max - min) / numBins;

        long[] observedFrequencies = new long[numBins];
        for (double value : data) {
            int binIndex = (int) ((value - min) / binWidth);
            if (binIndex >= numBins) binIndex = numBins - 1;
            observedFrequencies[binIndex]++;
        }

        double[] expectedFrequencies = new double[numBins];
        NormalDistribution normalDist = new NormalDistribution(mean(data), standardDeviation(data));
        for (int i = 0; i < numBins; i++) {
            double binStart = min + i * binWidth;
            double binEnd = binStart + binWidth;
            expectedFrequencies[i] = SAMPLE_SIZE * (normalDist.cumulativeProbability(binEnd) - normalDist.cumulativeProbability(binStart));
        }

        ChiSquareTest chiSquareTest = new ChiSquareTest();
        double chiSquareStat = chiSquareTest.chiSquare(expectedFrequencies, observedFrequencies);
        double pValue = chiSquareTest.chiSquareTest(expectedFrequencies, observedFrequencies);

        // Calculate critical Z value
        double alpha = 0.05; // Рівень значущості
        ChiSquaredDistribution chi2Dist = new ChiSquaredDistribution(numBins - 1); // ступені свободи
        double zCritical = chi2Dist.inverseCumulativeProbability(1 - alpha);

        // Display results
        System.out.println("Хі-квадрат статистика: " + chiSquareStat);
        System.out.println("p-значення: " + pValue);
        System.out.println("Критична точка Zкр: " + zCritical);

        if (chiSquareStat >= zCritical) {
            System.out.println("Відхиляємо нульову гіпотезу, дані не відповідають нормальному розподілу");
        } else {
            System.out.println("Не відхиляємо нульову гіпотезу, дані відповідають нормальному розподілу");
        }
    }

}
