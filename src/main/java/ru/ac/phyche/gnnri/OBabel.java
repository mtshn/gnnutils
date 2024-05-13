package ru.ac.phyche.gnnri;

import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Scanner;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

public class OBabel {

	public static String smilesTo3DSDFString(String smiles) throws IOException, CDKException {
		ProcessBuilder b = new ProcessBuilder("sh", "./ob.sh");
		b.directory(new File("./"));
		Process p = b.start();

		BufferedWriter w = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));

		w.write(smiles + "\n\n");
		w.flush();
		w.close();

		Scanner sc = new Scanner(p.getInputStream());
		ArrayList<String> l = new ArrayList<String>();
		boolean t = false;
		while (sc.hasNextLine()) {
			String s = sc.nextLine();
			String[] sp = s.split("\\s+");
			if (sp.length > 4) {
				if (sp[sp.length - 1].trim().equals("V2000")) {
					t = true;
				}
			}
			if (t) {
				l.add(s);
			}
			if (s.trim().equals("$$$$")) {
				t = false;
			}
		}
		sc.close();
		if (t) {
			throw (new CDKException("SDF mol created by openbabel is incorrect"));
		}
		String mol = "\nheader\n\n";
		for (String s : l) {
			mol = mol + s + "\n";
		}
		return mol;
	}

	public static IAtomContainer sdfStringToMol(String sdfString)
			throws IllegalArgumentException, CDKException, IOException {
		MDLV2000Reader r = new MDLV2000Reader(new ByteArrayInputStream(sdfString.getBytes(StandardCharsets.UTF_8)));
		IAtomContainer m = r.read(SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
		r.close();
		return m;
	}
                                                                             
}
