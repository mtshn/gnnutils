package ru.ac.phyche.gnnri;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class GraphC {

	private static final int DIGITS = 4;
	public static final DecimalFormat FORMAT = getDecimalFormat();

	public static class EdgeC {
		public int node1;
		public int node2;
		public float[] features;
	}

	public static class MSOfFragment {
		public int[] mzs;
		public float[] intens;

		public static MSOfFragment fromMap(HashMap<Integer, Float> m) {
			MSOfFragment result = new MSOfFragment();
			Integer[] a1 = m.keySet().toArray(new Integer[m.keySet().size()]);
			result.mzs = new int[a1.length];
			result.intens = new float[a1.length];
			for (int i = 0; i < a1.length; i++) {
				result.mzs[i] = a1[i];
			}
			Arrays.sort(result.mzs);
			for (int i = 0; i < a1.length; i++) {
				result.intens[i] = m.get(result.mzs[i]);
			}
			return result;
		}

		public String asString(float threshold) {
			if (mzs.length != intens.length) {
				throw (new RuntimeException("Wrong array length"));
			}
			String result = "";
			for (int i = 0; i < mzs.length; i++) {
				if (intens[i] > threshold) {
					result = result + mzs[i] + " " + FORMAT.format(intens[i]) + " ";
				}
			}
			return result.trim();
		}

		public HashMap<Integer, Float> asHashMap() {
			HashMap<Integer, Float> result = new HashMap<Integer, Float>();
			if (mzs.length != intens.length) {
				throw (new RuntimeException("Wrong array length"));
			}
			for (int i = 0; i < mzs.length; i++) {
				result.put(mzs[i], intens[i]);
			}
			return result;
		}

	}

	private int nEdgeFeatures = 0;
	private EdgeC[] edges;
	private float[][] nodeFeatures;
	private String comment;
	private int[] elementAtomicNumber;
	private int[] nHydrogens;

	public int getNNodes() {
		return nodeFeatures.length;
	}

	private static DecimalFormat getDecimalFormat() {
		DecimalFormat f = new DecimalFormat();
		f.setMaximumFractionDigits(DIGITS);
		DecimalFormatSymbols s = DecimalFormatSymbols.getInstance();
		s.setDecimalSeparator('.');
		f.setDecimalFormatSymbols(s);
		f.setGroupingUsed(false);
		return f;
	}

	public GraphC(int nNodes, int nNodeFeatures, int nEdges, int nEdgeFeatures) {
		this.edges = new EdgeC[nEdges];
		this.nEdgeFeatures = nEdgeFeatures;
		nodeFeatures = new float[nNodes][nNodeFeatures];
		elementAtomicNumber = new int[nNodes];
		nHydrogens = new int[nNodes];
	}

	public GraphC(int nNodes, int nNodeFeatures, int nEdges, int nEdgeFeatures, String comment) {
		this.edges = new EdgeC[nEdges];
		this.nEdgeFeatures = nEdgeFeatures;
		nodeFeatures = new float[nNodes][nNodeFeatures];
		this.comment = comment;
		elementAtomicNumber = new int[nNodes];
		nHydrogens = new int[nNodes];
	}

	public void setFeaturesOfNode(int n, int elementAtomicNumber, int nHydrogens, float[] features) {
		nodeFeatures[n] = features;
		this.elementAtomicNumber[n] = elementAtomicNumber;
		this.nHydrogens[n] = nHydrogens;
	}

	public void addEdge(int node1, int node2, float[] features) {
		if (features.length != nEdgeFeatures) {
			throw new RuntimeException();
		}
		int i = 0;
		while (edges[i] != null) {
			i++;
		}
		edges[i] = new EdgeC();
		edges[i].node1 = node1;
		edges[i].node2 = node2;
		edges[i].features = features;
	}

	public String graphToString() {
		String result = "";
		if (comment != null) {
			result = result + comment;
		}
		result = result + nodeFeatures.length + " " + nodeFeatures[0].length + " " + edges.length + " " + nEdgeFeatures;
		for (int i = 0; i < elementAtomicNumber.length; i++) {
			result = result + " " + elementAtomicNumber[i];
		}
		for (int i = 0; i < elementAtomicNumber.length; i++) {
			result = result + " " + nHydrogens[i];
		}
		for (int i = 0; i < edges.length; i++) {
			result = result + " " + edges[i].node1 + " " + edges[i].node2;
		}
		for (int i = 0; i < edges.length; i++) {
			for (int j = 0; j < nEdgeFeatures; j++) {
				result = result + " " + edges[i].features[j];
			}
		}
		for (int i = 0; i < nodeFeatures.length; i++) {
			for (int j = 0; j < nodeFeatures[0].length; j++) {
				result = result + " " + nodeFeatures[i][j];
			}
		}
		return result;
	}

	public String graphToStringMultiline() {
		return graphToStringMultiline(false);
	}

	public String graphToStringMultiline(boolean onlyAtoms) {
		String result = "";
		if (!onlyAtoms) {
			if (comment != null) {
				result = result + comment + "\n";
			}
			result = result + nodeFeatures.length + " " + nodeFeatures[0].length + " " + edges.length + " "
					+ nEdgeFeatures + "\n";
			for (int i = 0; i < edges.length; i++) {
				result = result + "edge " + i + " " + edges[i].node1 + " " + edges[i].node2 + "\n";
			}
			for (int i = 0; i < edges.length; i++) {
				result = result + "edge " + i;
				for (int j = 0; j < nEdgeFeatures; j++) {
					result = result + " " + edges[i].features[j];
				}
				result = result + "\n";
			}
		}
		for (int i = 0; i < nodeFeatures.length; i++) {
			result = result + "node " + i + " " + elementAtomicNumber[i] + " " + nHydrogens[i];
			for (int j = 0; j < nodeFeatures[0].length; j++) {
				result = result + " " + nodeFeatures[i][j];
			}
			result = result + "\n";
		}
		return result;
	}

	public String graphToStringMultilineOnlyEdges() {
		String result = "";
		if (comment != null) {
			result = result + comment + "\n";
		}
		for (int i = 0; i < edges.length; i++) {
			result = result + "edge " + i;
			for (int j = 0; j < nEdgeFeatures; j++) {
				result = result + " " + edges[i].features[j];
			}
			result = result + "\n";
		}
		result = result + "\n";
		return result;
	}

	public void removeNonBondedEdgesAndFeatures(int nOfFeatureThatMarksNonbondedEdge, boolean nonbondedMarkIsZero,
			int[] nonbondedFeatures) {
		ArrayList<EdgeC> result = new ArrayList<EdgeC>();
		float v = nonbondedMarkIsZero ? 0.0f : 1.0f;
		int nFeatures = edges[0].features.length;
		for (int i = 0; i < edges.length; i++) {
			EdgeC e = edges[i];
			if (e.features[nOfFeatureThatMarksNonbondedEdge] != v) {
				float[] newFeatures = new float[nFeatures - nonbondedFeatures.length];
				int j = 0;
				for (int k = 0; k < nFeatures; k++) {
					boolean found = false;
					for (int m = 0; m < nonbondedFeatures.length; m++) {
						if (k == nonbondedFeatures[m]) {
							found = true;
						}
					}
					if (!found) {
						newFeatures[j] = e.features[k];
						j++;
					}
				}
				e.features = newFeatures;
				result.add(e);
			}
		}
		EdgeC[] newEdges = result.toArray(new EdgeC[result.size()]);
		this.edges = newEdges;
		this.nEdgeFeatures=newEdges[0].features.length;
	}

}
