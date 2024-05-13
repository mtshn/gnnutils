package ru.ac.phyche.gnnri;

import java.io.IOException;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import io.github.dan2097.jnainchi.InchiStatus;

public class MolCheckUpUtils {
	public static final int MAXMW = 1000;
	public static final int SMILES_MAX_LENGTH = 300;
	
	public static String canonicalizeMol(String smiles) throws CDKException, IOException {
		IAtomContainer mol = smilesToMol(smiles);
		String inchikey14First = molToInchiKey(mol).split("\\-")[0];
		mol = preprocessMol(inchiToMol(molToInchi(preprocessMol(mol))));
		String inchikey14Second = molToInchiKey(mol).split("\\-")[0];
		String result = molToSmiles(mol);
		String inchikey14Third = molToInchiKey(smilesToMol(result)).split("\\-")[0];
		if (!inchikey14First.equals(inchikey14Second)) {
			throw (new CDKException("inchikey changed during conversion mol-inchi-mol using CDK"));
		}
		if (!inchikey14Third.equals(inchikey14Second)) {
			throw (new CDKException("inchikey changed during conversion mol-smiles-mol using CDK"));
		}
		return result;
	}

	public static String molToSmiles(IAtomContainer mol) throws CDKException {
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		AtomContainerManipulator.suppressHydrogens(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		for (IAtom a : mol.atoms()) {
			a.setMassNumber(null);
		}
		SmilesGenerator sg;
		sg = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols ^ SmiFlavor.AtomicMass);
		String smiles = sg.create(mol);
		return smiles;
	}

	public static IAtomContainer inchiToMol(String inchi) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIToStructure intostruct = factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
		InchiStatus ret = intostruct.getStatus();
		if (ret.equals(InchiStatus.ERROR)) {
			throw (new CDKException("Inchi status failed!"));
		}
		IAtomContainer mol = intostruct.getAtomContainer();
		return mol;
	}

	public static IAtomContainer smilesToMol(String smiles) throws CDKException {
		SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer mol = parser.parseSmiles(smiles.trim());
		return mol;
	}

	public static IAtomContainer preprocessMol(IAtomContainer mol) throws CDKException, IOException {
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		// AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
		for (IAtom a : mol.atoms()) {
			a.setMassNumber(Isotopes.getInstance().getMajorIsotope(a.getAtomicNumber()).getMassNumber());
		}
		Isotopes.getInstance().configureAtoms(mol);
		// AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		return mol;
	}


	public static IAtomContainer checkElements(String smiles) throws CDKException {
		IAtomContainer mol1 = smilesToMol(smiles);
		for (IAtom a : mol1.atoms()) {
			int n = a.getAtomicNumber();
			if (!((n == 1) || (n == 6) || (n == 7) || (n == 8) || (n == 9) || (n == 14) || (n == 15) || (n == 16)
					|| (n == 17) || (n == 35) || (n == 53))) {
				throw (new CDKException("Element type isn't supported. Only H,C,N,O,Si,S,P,F,Cl,Br,I! + Element: "
						+ a.getAtomicNumber() + a.getSymbol()));
			}
		}
		return mol1;
	}

	public static void checkConnectivityCycles(String smiles) throws CDKException {
		if (smiles.contains("0")) {
			throw (new CDKException("Not more than 9 cycles are allowed " + smiles));
		}
		if (smiles.contains(".")) {
			throw (new CDKException("Only one connectivity component is allowed " + smiles));
		}
	}

	public static float weight(IAtomContainer mol) throws CDKException {
		WeightDescriptor molarMass = new WeightDescriptor();
		return ((float) ((DoubleResult) molarMass.calculate(mol).getValue()).doubleValue());
	}

	public static void checkMW(String smiles) throws CDKException {
		IAtomContainer mol1 = smilesToMol(smiles);
		if (weight(mol1) > MAXMW) {
			throw (new CDKException("Too large MW " + weight(mol1) + "  " + smiles));
		}
	}

	public static void smilesCanCheck(String smilesCan) throws CDKException, IOException {
		if (!canonicalizeMol(smilesCan).equals(smilesCan)) {
			throw (new CDKException("Smiles changes after second canonicalization   " + smilesCan));
		}
		if (smilesCan.trim().length() > SMILES_MAX_LENGTH) {
			throw (new CDKException("Smiles too long   " + smilesCan.length() + " " + smilesCan));
		}
	}
	
	public static String molToInchiKey(IAtomContainer mol) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchiKey();
		return inchi;
	}

	public static String molToInchi(IAtomContainer mol) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchi();
		return inchi;
	}
	
	public static void checkMolIdenlity(IAtomContainer mol1, IAtomContainer mol2) throws CDKException {
		String inchikey14First = molToInchiKey(mol1).split("\\-")[0];
		String inchikey14Second = molToInchiKey(mol2).split("\\-")[0];
		if (!inchikey14First.equals(inchikey14Second)) {
			throw (new CDKException("inchikey changed during 3D coordinate generation using OB"));
		}
	}
	
}
