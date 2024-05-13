package ru.ac.phyche.gnnri;

import java.util.ArrayList;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

public class MolGraphUtils {
	public static float length(IBond b) {
		double dx = b.getBegin().getPoint3d().x - b.getEnd().getPoint3d().x;
		double dy = b.getBegin().getPoint3d().y - b.getEnd().getPoint3d().y;
		double dz = b.getBegin().getPoint3d().z - b.getEnd().getPoint3d().z;
		return (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

	public static float distance(IAtom a1, IAtom a2) {
		double dx = a1.getPoint3d().x - a2.getPoint3d().x;
		double dy = a1.getPoint3d().y - a2.getPoint3d().y;
		double dz = a1.getPoint3d().z - a2.getPoint3d().z;
		return (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

	public static boolean atomsLinked(IAtom a1, IAtom a2) {
		for (IBond b : a1.bonds()) {
			if (b.getOther(a1).equals(a2)) {
				return true;
			}
		}
		return false;
	}

	public static boolean atoms13(IAtom a1, IAtom a2) {
		for (IBond b1 : a1.bonds()) {
			for (IBond b2 : a2.bonds()) {
				if (!a1.equals(a2)) {
					if (b1.getOther(a1).equals(b2.getOther(a2))) {
						return true;
					}
				}
			}
		}
		return false;
	}

	public static ArrayList<IAtom[]> atoms14Atoms(IAtom a1, IAtom a2) {
		ArrayList<IAtom[]> result = new ArrayList<IAtom[]>();
		if (!a1.equals(a2)) {
			for (IBond b1 : a1.bonds()) {
				for (IBond b2 : a2.bonds()) {
					IAtom a3 = b1.getOther(a1);
					IAtom a4 = b2.getOther(a2);
					if ((!a3.equals(a4)) && (!a3.equals(a2)) && (!a4.equals(a1))) {
						if (atomsLinked(a3, a4)) {
							result.add(new IAtom[] { a1, a3, a4, a2 });
						}
					}
				}
			}
		}
		return result;
	}

	public static ArrayList<IAtom[]>  atoms15Atoms(IAtom a1, IAtom a2) {
		ArrayList<IAtom[]> result = new ArrayList<IAtom[]>();
		if (!a1.equals(a2)) {
			for (IBond b1 : a1.bonds()) {
				IAtom a3 = b1.getOther(a1);
				ArrayList<IAtom[]> links14 = atoms14Atoms(a3, a2);
				for (IAtom[] link14 : links14) {
					boolean x = true;
					for (int i = 0; i < link14.length; i++) {
						if (link14[i].equals(a1)) {
							x = false;
						}
					}
					if (x) {
						result.add(mergeArrays(new IAtom[] { a1 }, link14));
					}
				}
			}
		}
		return result;
	}

	public static ArrayList<IAtom[]>  atoms16Atoms(IAtom a1, IAtom a2) {
		ArrayList<IAtom[]> result = new ArrayList<IAtom[]>();
		if (!a1.equals(a2)) {
			for (IBond b1 : a1.bonds()) {
				IAtom a3 = b1.getOther(a1);
				ArrayList<IAtom[]> links15 = atoms15Atoms(a3, a2);
				for (IAtom[] link15 : links15) {
					boolean x = true;
					for (int i = 0; i < link15.length; i++) {
						if (link15[i].equals(a1)) {
							x = false;
						}
					}
					if (x) {
						result.add(mergeArrays(new IAtom[] { a1 }, link15));
					}
				}
			}
		}
		return result;
	}
	
	public static boolean atoms16(IAtom a1, IAtom a2) {
		return (atoms16Atoms(a1, a2).size() != 0);
	}

	public static boolean atoms14(IAtom a1, IAtom a2) {
		return (atoms14Atoms(a1, a2).size() != 0);
	}

	public static boolean atoms15(IAtom a1, IAtom a2) {
		return (atoms15Atoms(a1, a2).size() != 0);
	}

	public static int graphDistance(IAtom a1, IAtom a2) {
		if (atomsLinked(a1, a2)) {
			return 1;
		}
		if (atoms13(a1, a2)) {
			return 2;
		}
		if (atoms14(a1, a2)) {
			return 3;
		}
		if (atoms15(a1, a2)) {
			return 4;
		}
		if (atoms16(a1, a2)) {
			return 5;
		}
		return 0;
	}

	public static float[] inSame3456Ring(IAtom a1, IAtom a2) {
		float[] result = new float[4];
		boolean x12 = atomsLinked(a1, a2);
		boolean x13 = atoms13(a1, a2);
		boolean x14 = atoms14(a1, a2);
		boolean x15 = atoms15(a1, a2);
		boolean x16 = atoms16(a1, a2);

		if (x12 && x13) {
			result[0] = 1.0f;
		}
		if (x12 && x14) {
			result[1] = 1.0f;
		}
		if (x12 && x15) {
			result[2] = 1.0f;
		}
		if (x12 && x16) {
			result[3] = 1.0f;
		}
		if (x13 && x14) {
			result[2] = 1.0f;
		}
		if (x13 && x15) {
			result[3] = 1.0f;
		}
		if (x14) {
			if (atoms14Atoms(a1,a2).size() > 1) {
				result[3] = 1.0f;
			}
		}
		if (x13) {
			int count = 0;
			for (IBond b1 : a1.bonds()) {
				for (IBond b2 : a2.bonds()) {
					if (b1.getOther(a1).equals(b2.getOther(a2))) {
						count = count + 1;
					}
				}
			}
			if (count > 1) {
				result[1] = 1.0f;
			}
		}

		return result;
	}
	
	public static IAtom[] mergeArrays(IAtom a[], IAtom b[]) {
		IAtom[] result = new IAtom[a.length + b.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (i < a.length) ? a[i] : b[i - a.length];
		}
		return result;
	}

	
}
