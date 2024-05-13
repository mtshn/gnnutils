package ru.ac.phyche.gnnri;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.tuple.Pair;

public class GraphInputForRI {

	private static int[] ints(int n) {
		int[] ints = new int[n];
		for (int i = 0; i < ints.length; i++) {
			ints[i] = i;
		}
		Random rnd = new Random();
		for (int c = 0; c < 5; c++) {
			for (int i = 0; i < ints.length; i++) {
				int k = ints[i];
				int j = rnd.nextInt(ints.length);
				int m = ints[j];
				ints[j] = k;
				ints[i] = m;
			}
		}
		return ints;
	}

	public static Pair<GraphC[], String[]> createGraphsForMultipleMols(String[] smiles, String[] lines,
			boolean nonBondedEdges, float nonBondedCutOff) {
		GraphC[] result = new GraphC[smiles.length];
		String[] errors = new String[smiles.length];

		AtomicInteger ai = new AtomicInteger(0);
		Arrays.stream(ints(smiles.length)).parallel().forEach(i -> {
			int j = ai.incrementAndGet();
			if (j % 1000 == 0) {
				System.out.println("Computing features for GNN: " + j + " out of " + smiles.length);
			}
			try {
				result[i] = CreateGNNFeatures.molToMolGraph(smiles[i], lines[i] + " ", nonBondedEdges, nonBondedCutOff);
			} catch (Throwable e) {
				e.printStackTrace();
				errors[i] = smiles[i] + " " + e.getMessage();
			}
		});
		ArrayList<GraphC> resultList = new ArrayList<GraphC>();
		ArrayList<String> errorsList = new ArrayList<String>();
		for (int i = 0; i < smiles.length; i++) {
			if (result[i] != null) {
				resultList.add(result[i]);
			}
			if (errors[i] != null) {
				errorsList.add(errors[i]);
			}
		}
		GraphC[] result1 = resultList.toArray(new GraphC[resultList.size()]);
		String[] errors1 = errorsList.toArray(new String[errorsList.size()]);
		return Pair.of(result1, errors1);
	}

	public static void convertRIToGraphs(String outDir, String inputFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(inputFile));

		ArrayList<String> lib = new ArrayList<String>();
		String lineml = br.readLine();
		while (lineml != null) {
			if (!lineml.trim().equals("")) {
				lib.add(lineml.trim());
			}
			lineml = br.readLine();
		}
		br.close();
		int n = lib.size();
//		n = 20000;
		String[] smiles = new String[n];
		String[] lines = new String[n];
		for (int i = 0; i < n; i++) {
			smiles[i] = lib.get(i).split("\\s+")[0].trim();
			lines[i] = lib.get(i);
		}
		Files.createDirectories(Paths.get(outDir));
		Pair<GraphC[], String[]> p = createGraphsForMultipleMols(smiles, lines, true, 5);
		FileWriter fw = new FileWriter(new File(outDir, "graphs1nb.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToString() + "\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "graphs2nb.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToStringMultiline() + "\n\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "graphs3nb.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToStringMultiline(true) + "\n\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "graphs4nb.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToStringMultilineOnlyEdges() + "\n\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "errorsnb.txt"));
		for (int i = 0; i < p.getRight().length; i++) {
			fw.write(p.getRight()[i] + "\n");
		}
		fw.close();

		for (int i = 0; i < p.getLeft().length; i++) {
			p.getLeft()[i].removeNonBondedEdgesAndFeatures(7, false,
					new int[] { 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 });
		}

		fw = new FileWriter(new File(outDir, "graphs1.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToString() + "\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "graphs2.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToStringMultiline() + "\n\n");
		}
		fw.close();
		fw = new FileWriter(new File(outDir, "graphs3.txt"));
		for (int i = 0; i < p.getLeft().length; i++) {
			fw.write(p.getLeft()[i].graphToStringMultiline(true) + "\n\n");
		}
		fw.close();

		fw = new FileWriter(new File(outDir, "errors.txt"));
		for (int i = 0; i < p.getRight().length; i++) {
			fw.write(p.getRight()[i] + "\n");
		}

		fw.close();
		System.out.println("Complete!");

	}

	public static void main(String[] args) throws Exception {
		convertRIToGraphs("./GMD", "./dataSet/gmd1.ri");
		//convertRIToGraphs("./VAL", "./dataSet/val.ri");
		// convertRIToGraphs("/home/mat/julia/TEST", "/home/mat/julia/test.ri");
		// convertRIToGraphs("/home/mat/julia/TRAIN", "/home/mat/julia/train.ri");
	}
}
