package ru.ac.phyche.gnnri;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import org.openscience.cdk.charges.GasteigerMarsiliPartialCharges;
import org.openscience.cdk.charges.GasteigerPEPEPartialCharges;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.forcefield.mmff.Mmff;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.LonePairElectronChecker;

public class CreateGNNFeatures {

	// Only types presents for more then 0.02% atoms
	public static final String[] mmffAtomTypes = new String[] { "NC=N", "NC=O", "STHI", "NAZT", "NC=S", "PO2", "PO4",
			"PO3", "OM2", "NC%N", "C=C", "NSP", "C=SN", "C5A", "C5B", "CR3R", "CONN", "O=C", "C=O", "C=N", "S=C", "C=S",
			"NC=C", "OPO3", "OPO2", "O=N", "S=O", "O=S", "F", "I", "O2CM", "BR", "P", "CR4R", "C5", "S", "COOO", "NN=N",
			"SI", "SO3", "S-P", "SO2", "=C=", "SM", "OC=C", "CB", "CGD", "COO", "CL", "NPYD", "N3OX", "NSO2", "CR",
			"NPYL", "N=C", "N5B", "N5A", "N=N", "N=O", "OSO2", "OSO3", "OPO", "NO2", "OFUR", "PTET", "O=CR", "O=CN",
			"O=CO", "SO2N", "CE4R", "OC=N", "C=OS", "OC=O", "C=OR", "C=ON", "-O-", "=N=", "OC=S", "OS=O", "NR", "ONO2",
			"-OP", "O2NO", "OP", "CSP", "CSS", "OR", "NR+", "O2N", "O2S", "PO" };

	// Most common atom types by CDK
	public static final String[] mostCommonAtomTypes = new String[] { "C.sp3", "C.sp2", "O.sp3", "O.sp2", "N.sp2",
			"N.amide", "N.sp3", "C.sp", "S.3", "O.planar3", "N.nitro", "N.planar3", "S.planar3" };

	// Atom types by CDK (uncommon types - only element letter is retained)
	public static final String[] atomTypes = new String[] { "C", "F", "Cl", "I", "S.planar3", "C.sp2", "C.sp3", "N",
			"O", "P", "N.nitro", "S", "X", "S.3", "O.planar3", "N.planar3", "Br", "N.amide", "N.sp2", "Si", "N.sp3",
			"C.sp", "O.sp3", "O.sp2" };

	public static final String[] hybridizations = new String[] { "PLANAR3", "SP3D2", "SP2", "SP1", "SP3D1", "SP3" };

	public static float[] bondFeatures(IBond b, boolean fullAtomPairsFeatures) {
		float[] result = new float[fullAtomPairsFeatures ? 22 : 7];
		result[0] = MolGraphUtils.length(b);
		result[1] = b.isAromatic() ? 1.0f : 0.0f;
		result[2] = b.isInRing() ? 1.0f : 0.0f;
		result[3] = b.getOrder().numeric();
		result[4] = result[3] == 1 ? 1.0f : 0.0f;
		result[5] = result[3] == 2 ? 1.0f : 0.0f;
		result[6] = result[3] == 3 ? 1.0f : 0.0f;
		return result;
	}

	public static float[] atompairFeatures(IAtom a1, IAtom a2) {
		float[] result = new float[18];
		result[7] = 1;
		result[8] = MolGraphUtils.distance(a1, a2);
		result[9] = result[8] <= 3 ? 1.0f : 0.0f;
		result[10] = result[8] <= 3.5 ? 1.0f : 0.0f;
		result[11] = result[8] <= 4 ? 1.0f : 0.0f;
		result[12] = result[8] <= 4.5 ? 1.0f : 0.0f;
		result[13] = MolGraphUtils.atoms13(a1, a2) ? 1.0f : 0.0f;
		result[14] = MolGraphUtils.atoms14(a1, a2) ? 1.0f : 0.0f;
		result[15] = MolGraphUtils.atoms15(a1, a2) ? 1.0f : 0.0f;
		result[16] = MolGraphUtils.atoms16(a1, a2) ? 1.0f : 0.0f;
		result[17] = MolGraphUtils.graphDistance(a1, a2);
		result = mergeArrays(result, MolGraphUtils.inSame3456Ring(a1, a2));
		return result;
	}

	public static GraphC molToMolGraph(String smiles, String comment, boolean nonBondedEdges, float nonBondedCutOff)
			throws CDKException, IllegalArgumentException, IOException {
		MolCheckUpUtils.checkElements(smiles);
		MolCheckUpUtils.checkConnectivityCycles(smiles);
		MolCheckUpUtils.checkMW(smiles);
		String smilesCan = MolCheckUpUtils.canonicalizeMol(MolCheckUpUtils.canonicalizeMol(smiles));
		MolCheckUpUtils.smilesCanCheck(smilesCan);
		MolCheckUpUtils.checkConnectivityCycles(smilesCan);
		MolCheckUpUtils.checkMW(smilesCan);

		IAtomContainer mol = OBabel.sdfStringToMol(OBabel.smilesTo3DSDFString(smilesCan));
		mol = MolCheckUpUtils.preprocessMol(mol);
		MolCheckUpUtils.checkMolIdenlity(mol, MolCheckUpUtils.smilesToMol(smilesCan));
		MolCheckUpUtils.checkMolIdenlity(mol, MolCheckUpUtils.smilesToMol(smiles));

		int heavyAtoms = 0;

		HashMap<Integer, float[]> atomTypes1 = atomTypes(mol);
		HashMap<Integer, float[]> hybridization1 = hybridization(mol);
		HashMap<Integer, float[]> simpleAtomProperties1 = AtomicDescriptors.simpleAtomProperties(mol);
		HashMap<Integer, float[]> gasteigerCharges = sigmaPiGasteigerCharges(mol);
		HashMap<Integer, float[]> descriptors = AtomicDescriptors.selectedAtomicDescriptors(mol);
		HashMap<Integer, float[]> atomTypesMMFF = atomTypesMMFF94(mol);

		HashMap<Integer, Integer> graphNumberByCDKNumber = new HashMap<Integer, Integer>();
		ArrayList<float[]> atomsFeatures = new ArrayList<float[]>();
		ArrayList<Integer> atomsElementAtomicNumber = new ArrayList<Integer>();
		ArrayList<Integer> atomsNHydrogensNeigh = new ArrayList<Integer>();
		ArrayList<float[]> bondsFeatures = new ArrayList<float[]>();
		ArrayList<Integer> firstNodesOfBonds = new ArrayList<Integer>();
		ArrayList<Integer> secondNodesOfBonds = new ArrayList<Integer>();

		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				heavyAtoms++;
				int n = a.getIndex();
				graphNumberByCDKNumber.put(n, atomsFeatures.size());
				float[] features = mergeArrays(
						mergeArrays(atomTypes1.get(n), hybridization1.get(n), simpleAtomProperties1.get(n)),
						mergeArrays(gasteigerCharges.get(n), descriptors.get(n), atomTypesMMFF.get(n)));
				for (int i = 0; i < features.length; i++) {
					if (Float.isNaN(features[i]) || Float.isInfinite(features[i])) {
						features[i] = 0;
					}
				}
				atomsFeatures.add(features);
				atomsElementAtomicNumber.add(a.getAtomicNumber());
				atomsNHydrogensNeigh.add(AtomicDescriptors.nHydrogensNeigh(a));
			}
		}
		if (heavyAtoms != atomsFeatures.size()) {
			throw new RuntimeException("Unknown error");
		}

		for (IBond b : mol.bonds()) {
			boolean nonH = true;
			for (IAtom a : b.atoms()) {
				if (a.getAtomicNumber().equals(1)) {
					nonH = false;
				}
			}
			if (nonH) {
				bondsFeatures.add(bondFeatures(b, nonBondedEdges));
				int a1 = graphNumberByCDKNumber.get(b.getBegin().getIndex());
				int a2 = graphNumberByCDKNumber.get(b.getEnd().getIndex());
				int a11 = a1 < a2 ? a1 : a2;
				int a21 = a1 < a2 ? a2 : a1;
				if (a1 == a2) {
					throw new RuntimeException("Unknown error");
				}
				firstNodesOfBonds.add(a11);
				secondNodesOfBonds.add(a21);
			}
		}

		if (nonBondedEdges) {
			for (IAtom at1 : mol.atoms()) {
				if (!at1.getAtomicNumber().equals(1)) {
					for (IAtom at2 : mol.atoms()) {
						if (!at2.getAtomicNumber().equals(1)) {
							int a1 = at1.getIndex();
							int a2 = at2.getIndex();
							if (a2 > a1) {
								if (!at1.equals(at2)) {
									if (MolGraphUtils.distance(at1, at2) < nonBondedCutOff) {
										if (!MolGraphUtils.atomsLinked(at1, at2)) {
											a1 = graphNumberByCDKNumber.get(a1);
											a2 = graphNumberByCDKNumber.get(a2);
											int a11 = a1 < a2 ? a1 : a2;
											int a21 = a1 < a2 ? a2 : a1;
											if (a1 == a2) {
												throw new RuntimeException("Unknown error");
											}
											bondsFeatures.add(atompairFeatures(at1, at2));
											firstNodesOfBonds.add(a11);
											secondNodesOfBonds.add(a21);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		int nAtomFeatures = atomsFeatures.get(0).length;
		int nBondFeatures = bondsFeatures.get(0).length;
		GraphC result = new GraphC(heavyAtoms, nAtomFeatures, bondsFeatures.size(), nBondFeatures, comment);
		for (int i = 0; i < atomsFeatures.size(); i++) {
			result.setFeaturesOfNode(i, atomsElementAtomicNumber.get(i), atomsNHydrogensNeigh.get(i),
					atomsFeatures.get(i));
		}
		for (int i = 0; i < bondsFeatures.size(); i++) {
			result.addEdge(firstNodesOfBonds.get(i), secondNodesOfBonds.get(i), bondsFeatures.get(i));
		}

		return result;
	}

	public static float[] mergeArrays(float a[], float b[]) {
		float[] result = new float[a.length + b.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (i < a.length) ? a[i] : b[i - a.length];
		}
		return result;
	}

	public static float[] mergeArrays(float a[], float b[], float c[]) {
		return mergeArrays(a, mergeArrays(b, c));
	}

	private static HashMap<Integer, float[]> sigmaPiGasteigerCharges(IAtomContainer mol) throws CDKException {

		HashMap<Integer, Float> piGastCharges = new HashMap<Integer, Float>();
		HashMap<Integer, Float> sigmaGastCharges = new HashMap<Integer, Float>();

		LonePairElectronChecker lpcheck = new LonePairElectronChecker();
		for (IAtom a : mol.atoms()) {
			a.setCharge(0.0);
		}

		GasteigerPEPEPartialCharges pepe = new GasteigerPEPEPartialCharges();
		try {
			pepe.assignGasteigerPiPartialCharges(mol, true);
		} catch (Exception e) {
			throw (new CDKException("Gasteiger pi charges failed: " + e.getMessage()));
		}
		lpcheck.saturate(mol);

		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				piGastCharges.put((Integer) a.getIndex(), (Float) (float) (double) a.getCharge());
			}
		}

		for (IAtom a : mol.atoms()) {
			a.setCharge(0.0);
		}

		GasteigerMarsiliPartialCharges peoe = new GasteigerMarsiliPartialCharges();
		try {
			peoe.assignGasteigerMarsiliSigmaPartialCharges(mol, true);
		} catch (Exception e) {
			throw (new CDKException("Gasteiger charges failed: " + e.getMessage()));
		}

		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				sigmaGastCharges.put((Integer) a.getIndex(), (Float) (float) (double) a.getCharge());
			}
		}

		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				int i = a.getIndex();
				float sigma = sigmaGastCharges.get(i);
				float pi = piGastCharges.get(i);
				result.put(i, new float[] { sigma, pi });
			}
		}
		return result;
	}

	private static String atomTypeToString(IAtom at) {
		HashSet<String> mostCommonAtomTypesSet = new HashSet<String>();
		mostCommonAtomTypesSet.addAll(Arrays.asList(mostCommonAtomTypes));
		String atomType = at.getAtomTypeName();
		if (mostCommonAtomTypesSet.contains(atomType.trim())) {
			return atomType.trim();
		} else {
			return atomType.split("\\.")[0];
		}
	}

	private static float[] atomType(IAtom at) throws CDKException {
		int atType = -1;
		String atomType = atomTypeToString(at);
		for (int i = 0; i < atomTypes.length; i++) {
			if (atomTypes[i].equals(atomType)) {
				atType = i;
			}
		}
		if (atType == -1) {
			throw new CDKException("Unexpected atom type " + at.getAtomTypeName());
		}
		float[] result = new float[atomTypes.length];
		result[atType] = 1.0f;
		return result;
	}

	// Call before MMFF94! MMFF94 changes mol irreversible
	private static HashMap<Integer, float[]> atomTypes(IAtomContainer mol) throws CDKException {
		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				result.put(a.getIndex(), atomType(a));
			}
		}
		return result;
	}

	// Call after atomTypes! MMFF94 changes mol irreversible
	private static HashMap<Integer, float[]> atomTypesMMFF94(IAtomContainer mol) throws CDKException {
		Mmff n = new Mmff();
		n.assignAtomTypes(mol);
		String[] t = mmffAtomTypes;
		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				float[] f = new float[t.length];
				String ty = a.getAtomTypeName();
				for (int i = 0; i < t.length; i++) {
					if (ty.trim().equals(t[i].trim())) {
						f[i] = 1.0f;
					}
				}
				result.put(a.getIndex(), f);
			}
		}
		return result;
	}

	private static HashMap<Integer, float[]> hybridization(IAtomContainer mol) throws CDKException {
		String[] t = hybridizations;
		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				float[] f = new float[t.length];
				String ty = a.getHybridization().toString();
				for (int i = 0; i < t.length; i++) {
					if (ty.trim().equals(t[i].trim())) {
						f[i] = 1.0f;
					}
				}
				result.put(a.getIndex(), f);
			}
		}
		return result;
	}

}
