package meanPartition;

import java.util.Arrays;
import java.nio.file.Paths;

public class IrisData_test {
	
	private static void print2DimArray(int[][] arr){	
		for(int i =0;i<arr.length;i++)
			System.out.println(Arrays.toString(arr[i]));	
	}
	
	public static String arrayToRVector(String arrayName, double[] arr, int numberLimit) {
		String retStr = arrayName + " = c(";
		int maxNumbersPerLine = 200, numbersCounter = 0;
		for(int i = 0;i < arr.length;i++) {
			numbersCounter++;
			retStr += arr[i];
			if(i == numberLimit)
				break;
			if(i < arr.length - 1) {
				retStr += ",";
				if(numbersCounter == maxNumbersPerLine) {
					retStr += "\n";
					numbersCounter = 0;
				}
			}
		}
		retStr += ")";
		
		return(retStr);
	}
	
	public static String observationsToRVectors(double[][] observations) {
		String str = "";
		for(int i = 0; i < observations[0].length;i++) {
			double[] v = new double[observations.length];
			for(int j = 0; j < v.length;j++) {
				v[j] = observations[j][i];
			}
			str += arrayToRVector("x" + i,v,4000) + "\n";			
		}
		return str;
	}
	
	public static double getMean(double[] values) {
		double sumOf = 0;
		for(int i = 0;i < values.length;i++)
			sumOf += values[i];
		return sumOf / values.length;
	}
	
	public static double getBesselCorrectedVar(double[] values) {
		double mean = getMean(values);
		double sumOfSquaredDevs = 0;
		for(int i = 0;i < values.length;i++)
			sumOfSquaredDevs += (values[i] - mean) * (values[i] - mean) ;
		return sumOfSquaredDevs / (values.length - 1);
	}
	
	public static double getBesselCorrectedStd(double[] values) {
		return Math.sqrt(getBesselCorrectedVar(values));
	}
	
	//Calculates the length of a half bar of the 95% confidence interval according to:
	//https://de.wikipedia.org/wiki/Konfidenzintervall#.C3.9Cbersicht_f.C3.BCr_stetige_Verteilungen
	public static double getHalfBarLengthOf95CI(double[] values) {
		double ZValueFor5PercentAlpha = 1.959964;
		
		//calculate standard error
		double stdErr = getBesselCorrectedStd(values) / Math.sqrt(values.length);
		
		//calculate half bar length
		double halfBarLength = ZValueFor5PercentAlpha * stdErr;
		
		return halfBarLength;		
	}
	
	public static String generatePlotCodeForClustering(String datafile, int[] partition) {
		String str = "";
		
		//build absolute path
		datafile = Paths.get(".").toAbsolutePath().normalize().toString() + "\\" + datafile;
		datafile = datafile.replace("\\", "\\\\");
		
		str += "data = read.csv('" + datafile + "', header = TRUE, sep = ',')\n";
		str += "x0 = data[,1]\n";
		str += "x1 = data[,2]\n";
		str += "x2 = data[,3]\n";
		str += "x3 = data[,4]\n";		
		str += arrayToRVector("z",Arrays.stream(partition).asDoubleStream().toArray(), 4000) + "\n";
		
		str += "par(mfrow=c(4,4))\n";
		str += "plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab='', xlab='', xaxt='n',yaxt='n')\n";
		str += "legend('center', c('Sepal.length'), bty = 'n')\n";
		str += "plot(x1, x0, col = z, ylab='', xlab='')\n";
		str += "plot(x2, x0, col = z, ylab='', xlab='')\n";
		str += "plot(x3, x0, col = z, ylab='', xlab='')\n";
		str += "plot(x0, x1, col = z, ylab='', xlab='')\n";
		str += "plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab='', xlab='', xaxt='n',yaxt='n')\n";
		str += "legend('center', c('Sepal.Width'), bty = 'n')\n";
		str += "plot(x2, x1, col = z, ylab='', xlab='')\n";
		str += "plot(x3, x1, col = z, ylab='', xlab='')\n";
		str += "plot(x0, x2, col = z, ylab='', xlab='')\n";
		str += "plot(x1, x2, col = z, ylab='', xlab='')\n";
		str += "plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab='', xlab='', xaxt='n',yaxt='n')\n";
		str += "legend('center', c('Petal.Length'), bty = 'n')\n";
		str += "plot(x3, x2, col = z, ylab='', xlab='')\n";
		str += "plot(x0, x3, col = z, ylab='', xlab='')\n";
		str += "plot(x1, x3, col = z, ylab='', xlab='')\n";
		str += "plot(x2, x3, col = z, ylab='', xlab='')\n";
		str += "plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab='', xlab='', xaxt='n',yaxt='n')\n";
		str += "legend('center', c('Petal.Width'), bty = 'n')";

		return str;	
	}
	
	public static int getNumOfClusters(int[] partition) {
		boolean[] doesExist  = new boolean[partition.length];
		int numOfClusters = 0;
		for(int i = 0; i < partition.length;i++)
			if(!doesExist[partition[i]]) {
				doesExist[partition[i]] = true;
				numOfClusters++;
			}
		return numOfClusters;
	}
	
	public static void main(String[] args) {
		
		//First, the parameters for the base distribution ...
		int numOfDims = 4;
		double[] priorMu = 		new double[numOfDims];
		double priorKappa = 	numOfDims;		
		double priorNu = 		numOfDims;
		double[][] priorPsi = new double[numOfDims][numOfDims];
		for(int i = 0; i< numOfDims;i++)
			priorPsi[i][i] = 1;
	
		//specifiy here the datafile to use
		String[] fileToSelect = new String[6];
		fileToSelect[0] = "Datasets\\irisdata_50.csv";
		fileToSelect[1] = "Datasets\\irisdata_250.csv";
		fileToSelect[2] = "Datasets\\irisdata_450.csv";
		fileToSelect[3] = "Datasets\\irisdata_650.csv";
		fileToSelect[4] = "Datasets\\irisdata_850.csv";
		fileToSelect[5] = "Datasets\\irisdata_1050.csv";
		String datafile = fileToSelect[0];
		
		double[][] observations = CSV.read(datafile);
		
		//...then the concentration parameter for the DP
		double alpha = 1;

		//sampling settings
		int numOfSamples = 100;
		int thinning = 1;
		
		DPMixture dpm = new DPMixture(new DependentNormalsBase(priorMu, priorPsi, priorKappa, priorNu, observations), alpha);
		
		int[][] samples = dpm.simulate(0, numOfSamples, thinning);
		//print2DimArray(samples);
		
		//calculate mean number of clusters in all partitions
		double meanNumOfClusters = 0;
		for(int i = 0; i < samples.length; i++)
			meanNumOfClusters += getNumOfClusters(samples[i]);
		meanNumOfClusters /= samples.length;
	
		int[] correspondingMeanPartition = null;
		int numOfMeasurements = 100;
		double[] elapsedTime1 = new double[numOfMeasurements], elapsedTime2 =  new double[numOfMeasurements];
		int sumOfMPDistances = 0;
		for(int measurement = 0; measurement < numOfMeasurements; measurement++) {
			System.out.println("start measurement: " + (measurement + 1));
		
			long startTime, curElapsedTime1, curElapsedTime2;
			PartitionDistribution PD1 = new PartitionDistribution(samples);
			PartitionDistribution PD2 = new PartitionDistribution(samples);
			PD1.prepareArrayToKeepIterations(1);
			PD2.prepareArrayToKeepIterations(1);
			
			startTime = System.nanoTime();
			
			int[] meanPartition1 = PD1.getMeanPartition1();
			
			curElapsedTime1 = System.nanoTime() - startTime;
			startTime = System.nanoTime();
			
			int[] meanPartition2 = PD2.getMeanPartition2();
			
			curElapsedTime2 = System.nanoTime() - startTime;
		
			sumOfMPDistances += PD1.getPartitionDistance(meanPartition1, meanPartition2);			
			if(sumOfMPDistances != 0) {		
				System.out.print("\nmean partition 1 is:\n");
				System.out.println(Arrays.toString(meanPartition1) + "\n");
				
				System.out.print("\nmean partition 2 is:\n");
				System.out.println(Arrays.toString(meanPartition2) + "\n");
			}else
				correspondingMeanPartition = meanPartition2;
			elapsedTime1[measurement] += curElapsedTime1 / 1000000;
			elapsedTime2[measurement] += curElapsedTime2 / 1000000;			
		}
		double meanElapsedTime1 = getMean(elapsedTime1);
		double halfBar95CIForElapsedTime1 = getHalfBarLengthOf95CI(elapsedTime1);
		double meanElapsedTime2 = getMean(elapsedTime2);
		double halfBar95CIForElapsedTime2 = getHalfBarLengthOf95CI(elapsedTime2);
		System.out.println("\nmean number of clusters in all sample partitions: " + meanNumOfClusters);
		System.out.println("partition distance between both mean partitions: " + sumOfMPDistances);
		System.out.println("Mean elapsed time (in ms) for original algorithm: " + meanElapsedTime1 + " (+/- " + halfBar95CIForElapsedTime1 + ")");
		System.out.println("Mean elapsed time (in ms) for improved algorithm: " + meanElapsedTime2 + " (+/- " + halfBar95CIForElapsedTime2 + ")");
		System.out.println("Compared to the normal version, the 1. improved version is " + meanElapsedTime1 / meanElapsedTime2 + " times faster\n");
				
		System.out.println(generatePlotCodeForClustering(datafile, correspondingMeanPartition));
		System.out.println("");	
	}
}
