package ru.ac.phyche.gnnri;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.openscience.cdk.IImplementationSpecification;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;

public class AtomicDescriptors {
	
	private static final List<IDescriptor> atomicDescriptorList = getDescriptorList();

	// All atomic descriptors (this array isn't actually used)
	public static final String[] allAtomicDescriptorsNames = new String[] { "pepeT", "pepe", "g3r_1", "g3r_2", "g3r_3",
			"g3r_4", "g3r_5", "g3r_6", "g3r_7", "g3r_8", "g3r_9", "g3r_10", "g3r_11", "g3r_12", "g3r_13", "aNeg",
			"indAtomSoftness", "hybr", "periodicTablePosition", "gHrTop_1", "gHrTop_2", "gHrTop_3", "gHrTop_4",
			"gHrTop_5", "gHrTop_6", "gHrTop_7", "gHrTop_8", "gHrTop_9", "gHrTop_10", "gHrTop_11", "gHrTop_12",
			"gHrTop_13", "gHrTop_14", "gHrTop_15", "elecSigmA", "protonTotalPartialCharge1",
			"protonTotalPartialCharge2", "protonTotalPartialCharge3", "protonTotalPartialCharge4",
			"protonTotalPartialCharge5", "protonAffiHOSE", "stabilPlusC", "gSr_1", "gSr_2", "gSr_3", "gSr_4", "gSr_5",
			"gSr_6", "gSr_7", "elecPiA", "bondsToAtom", "indAtomHardnesss", "vdwRadius", "distanceToAtom",
			"protonInArmaticSystem", "protonInConjSystem", "aHyb", "effAtomPol", "partialSigmaCharge", "val",
			"RDF_GHR_0", "RDF_GHR_1", "RDF_GHR_2", "RDF_GHR_3", "RDF_GHR_4", "RDF_GHR_5", "RDF_GHR_6", "RDF_GHR_7",
			"RDF_GHR_8", "RDF_GHR_9", "RDF_GHR_10", "RDF_GHR_11", "RDF_GHR_12", "RDF_GHR_13", "RDF_GHR_14", "gDr_1",
			"gDr_2", "gDr_3", "gDr_4", "gDr_5", "gDr_6", "gDr_7", "covalentRadius", "ipAtomicHOSE", "partialTCMMFF94",
			"ipAtomicLearning" };

	// Actually used atomic desriptors
	public static final String[] selectedAtomicDescriptors = new String[] { "aNeg", "indAtomSoftness", "hybr",
			"periodicTablePosition", "elecSigmA", "protonTotalPartialCharge1", "protonTotalPartialCharge2",
			"protonTotalPartialCharge3", "protonTotalPartialCharge4", "protonAffiHOSE", "stabilPlusC",
			"indAtomHardnesss", "vdwRadius", "aHyb", "effAtomPol", "val", "covalentRadius", "ipAtomicHOSE",
			"partialTCMMFF94" };
	
	private static List<IDescriptor> getDescriptorList() {
		List<String> classes = DescriptorEngine
				.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.atomic", null);
		DescriptorEngine engine = new DescriptorEngine(classes, null);
		List<IDescriptor> inst = engine.instantiateDescriptors(classes);
		List<IImplementationSpecification> specs = engine.initializeSpecifications(inst);
		engine.setDescriptorInstances(inst);
		engine.setDescriptorSpecifications(specs);
		return engine.getDescriptorInstances();
	}
	
	private static int nHydrogensNeigh13(IAtom a) {
		int count = 0;
		for (IBond b : a.bonds()) {
			if (!b.getOther(a).getAtomicNumber().equals(1)) {
				count = count + nHydrogensNeigh(b.getOther(a));
			}
		}
		return count;
	}

	private static float[] centerMassMol(IAtomContainer mol) {
		float[] result = new float[3];
		float m = 0;
		for (IAtom a : mol.atoms()) {
			result[0] = result[0] + (float) (a.getExactMass() * a.getPoint3d().x);
			result[1] = result[1] + (float) (a.getExactMass() * a.getPoint3d().y);
			result[2] = result[2] + (float) (a.getExactMass() * a.getPoint3d().z);
			m = m + (float) (double) a.getExactMass();
		}
		result[0] = result[0] / m;
		result[1] = result[1] / m;
		result[2] = result[2] / m;
		return result;
	}

	private static float[] centerMol(IAtomContainer mol) {
		float[] result = new float[3];
		float m = 0;
		for (IAtom a : mol.atoms()) {
			result[0] = result[0] + (float) (a.getPoint3d().x);
			result[1] = result[1] + (float) (a.getPoint3d().y);
			result[2] = result[2] + (float) (a.getPoint3d().z);
			m = m + 1;
		}
		result[0] = result[0] / m;
		result[1] = result[1] / m;
		result[2] = result[2] / m;
		return result;
	}

	private static float distanceToCenter(IAtom a, float[] center) {
		float dx = (float) (a.getPoint3d().x) - center[0];
		float dy = (float) (a.getPoint3d().y) - center[1];
		float dz = (float) (a.getPoint3d().z) - center[2];
		return (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

	public static int nHydrogensNeigh(IAtom a) {
		int count = 0;
		for (IBond b : a.bonds()) {
			if (b.getOther(a).getAtomicNumber().equals(1)) {
				count = count + 1;
			}
		}
		if (a.getImplicitHydrogenCount() != 0) {
			throw (new RuntimeException("Implicit hydrogens are not allowed here"));
		}
		return count;
	}	
	
	public static HashMap<Integer, float[]> simpleAtomProperties(IAtomContainer mol) throws CDKException {
		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		float[] center1 = centerMol(mol);
		float[] center2 = centerMassMol(mol);
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				float atomicNumber = a.getAtomicNumber();
				float isC = atomicNumber == 6 ? 1.0f : 0.0f;
				float isN = atomicNumber == 7 ? 1.0f : 0.0f;
				float isO = atomicNumber == 8 ? 1.0f : 0.0f;
				float isF = atomicNumber == 9 ? 1.0f : 0.0f;
				float isSi = atomicNumber == 14 ? 1.0f : 0.0f;
				float isP = atomicNumber == 15 ? 1.0f : 0.0f;
				float isS = atomicNumber == 16 ? 1.0f : 0.0f;
				float isCl = atomicNumber == 17 ? 1.0f : 0.0f;
				float isBr = atomicNumber == 35 ? 1.0f : 0.0f;
				float isI = atomicNumber == 53 ? 1.0f : 0.0f;
				float isHalogen = ((isF == 1.0f) || (isCl == 1.0f) || (isBr == 1.0f) || (isI == 1.0f)) ? 1.0f : 0.0f;
				float bondCount = a.getBondCount();
				float bondOrderSum = (float) (double) a.getBondOrderSum();
				float charge = a.getFormalCharge();
				float valency = a.getValency();
				float exactMass = (float) (double) a.getExactMass();
				float hydrogens = nHydrogensNeigh(a);
				float oneBond = (bondCount == 1.0f) ? 1.0f : 0.0f;
				float twoBonds = (bondCount == 2.0f) ? 1.0f : 0.0f;
				float threeBonds = (bondCount == 3.0f) ? 1.0f : 0.0f;
				float fourBonds = (bondCount == 4.0f) ? 1.0f : 0.0f;
				float noH = (hydrogens == 0.0f) ? 1.0f : 0.0f;
				float oneH = (hydrogens == 1.0f) ? 1.0f : 0.0f;
				float twoH = (hydrogens == 2.0f) ? 1.0f : 0.0f;
				float threeH = (hydrogens == 3.0f) ? 1.0f : 0.0f;
				float hydrogens13 = nHydrogensNeigh13(a);
				float distance1 = distanceToCenter(a, center1);
				float distance2 = distanceToCenter(a, center2);
				float aromatic = a.isAromatic() ? 1.0f : 0.0f;
				float inRing = a.isInRing() ? 1.0f : 0.0f;
				float[] f = new float[] { atomicNumber, isC, isN, isO, isF, isSi, isP, isS, isCl, isBr, isI, isHalogen,
						bondCount, bondOrderSum, charge, valency, exactMass, hydrogens, oneBond, twoBonds, threeBonds,
						fourBonds, noH, oneH, twoH, threeH, hydrogens13, distance1, distance2, aromatic, inRing };
				result.put(a.getIndex(), f);
			}
		}
		return result;
	}

	public static float[] atomicDescriptors(IAtom at, IAtomContainer mol, String[] descriptorNames)
			throws CDKException {
		HashSet<String> descriptorNamesSet = new HashSet<String>(Arrays.asList(descriptorNames));
		HashMap<String, Float> resultMap = new HashMap<String, Float>();
		List<IDescriptor> descriptors = atomicDescriptorList;

		for (IDescriptor d : descriptors) {
			IAtomicDescriptor md = (IAtomicDescriptor) d;
			HashSet<String> thisDescriptorNamesSet = new HashSet<String>(Arrays.asList(md.getDescriptorNames()));
			thisDescriptorNamesSet.retainAll(descriptorNamesSet);

			if (thisDescriptorNamesSet.size() > 0) {
				boolean failedToCalculate = false;
				DescriptorValue dv = null;
				try {
					dv = computeDescriptor(at, mol, md);
				} catch (Throwable e) {
					failedToCalculate = true;
					if ((md == null) || (mol == null) || (at == null)) {
						throw e;
					}
				}
				if (dv == null) {
					failedToCalculate = true;
				}
				String[] values = null;
				if (!failedToCalculate) {
					IDescriptorResult dr = dv.getValue();
					if (dr != null) {
						values = dr.toString().split(",");
					} else {
						failedToCalculate = true;
					}
				}
				String[] names = md.getDescriptorNames();

				for (int j = 0; j < names.length; j++) {
					float value = Float.NaN;
					if (!failedToCalculate) {
						if (values[j].trim().equals("true")) {
							value = 1;
						} else {
							if (values[j].trim().equals("false")) {
								value = 1;
							} else {
								value = Float.parseFloat(values[j]);
							}
						}
					}
					resultMap.put(names[j], value);
				}
			}
		}
		float[] result = new float[descriptorNames.length];

		for (int k = 0; k < descriptorNames.length; k++) {
			if ((resultMap.get(descriptorNames[k]) == null)) {
				result[k] = Float.NaN;
			} else {
				float res = resultMap.get(descriptorNames[k]);
				result[k] = Float.isInfinite(res) ? Float.NaN : res;
			}
		}
		return result;
	}

	public static HashMap<Integer, float[]> selectedAtomicDescriptors(IAtomContainer mol) throws CDKException {
		HashMap<Integer, float[]> result = new HashMap<Integer, float[]>();
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				result.put(a.getIndex(), atomicDescriptors(a, mol, selectedAtomicDescriptors));
			}
		}
		return result;
	}

	private static DescriptorValue computeDescriptor(IAtom at, IAtomContainer mol, IAtomicDescriptor md) {
		try {
			DescriptorValue rslt = md.calculate(at, mol);
			return rslt;
		} catch (Throwable e) {
			return null;
		}
	}
}
