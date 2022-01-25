package fields.vectors;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.vectors.GenericPIDModule.Mod;

class GenericPIDModuleTest {

	@Test
	void testZ4() {
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<IntE> syzygy1 = new ArrayList<>();
		syzygy1.add(z.getInteger(1));
		syzygy1.add(z.getInteger(-2));
		List<IntE> syzygy2 = new ArrayList<>();
		syzygy2.add(z.getInteger(2));
		syzygy2.add(z.getInteger(0));
		List<IntE> syzygy3 = new ArrayList<>();
		syzygy3.add(z.getInteger(0));
		syzygy3.add(z.getInteger(4));
		syzygies.add(syzygy1);
		syzygies.add(syzygy2);
		syzygies.add(syzygy3);
		GenericPIDModule<IntE, Vector<IntE>> module = GenericPIDModule.fromSyzygies(z, 2, syzygies);
		List<IntE> baseVector1 = new ArrayList<>();
		baseVector1.add(z.one());
		baseVector1.add(z.zero());
		List<IntE> baseVector2 = new ArrayList<>();
		baseVector2.add(z.zero());
		baseVector2.add(z.one());
		Vector<IntE> base1 = new Vector<>(baseVector1);
		Vector<IntE> base2 = new Vector<>(baseVector2);
		System.out.println(module.reduce(base1));
		System.out.println(module.reduce(base2));
	}

	@Test
	void testZ6() {
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<IntE> syzygy1 = new ArrayList<>();
		syzygy1.add(z.getInteger(3));
		syzygy1.add(z.getInteger(0));
		List<IntE> syzygy2 = new ArrayList<>();
		syzygy2.add(z.getInteger(0));
		syzygy2.add(z.getInteger(2));
		syzygies.add(syzygy1);
		syzygies.add(syzygy2);
		GenericPIDModule<IntE, Vector<IntE>> module = GenericPIDModule.fromSyzygies(z, 2, syzygies);
		for (Mod<Vector<IntE>> basis : module.getModuleGenerators()) {
			System.out.println(basis);
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 2; j++) {
				List<IntE> vector = new ArrayList<>();
				vector.add(z.getInteger(i));
				vector.add(z.getInteger(j));
				Vector<IntE> v = new Vector<>(vector);
				System.out.println("Not reduced: " + v);
				System.out.println("Reduced: " + module.reduce(v));
			}
		}
	}

	@Test
	void testZ2Z4() {
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<IntE> syzygy1 = new ArrayList<>();
		syzygy1.add(z.getInteger(4));
		syzygy1.add(z.getInteger(0));
		List<IntE> syzygy2 = new ArrayList<>();
		syzygy2.add(z.getInteger(0));
		syzygy2.add(z.getInteger(4));
		List<IntE> syzygy3 = new ArrayList<>();
		syzygy3.add(z.getInteger(2));
		syzygy3.add(z.getInteger(-2));
		syzygies.add(syzygy1);
		syzygies.add(syzygy2);
		syzygies.add(syzygy3);
		GenericPIDModule<IntE, Vector<IntE>> module = GenericPIDModule.fromSyzygies(z, 2, syzygies);
		for (Mod<Vector<IntE>> basis : module.getModuleGenerators()) {
			System.out.println(basis);
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				List<IntE> vector = new ArrayList<>();
				vector.add(z.getInteger(i));
				vector.add(z.getInteger(j));
				Vector<IntE> v = new Vector<>(vector);
				System.out.println("Not reduced: " + v);
				System.out.println("Reduced: " + module.reduce(v));
			}
		}
		System.out.println();
		for (Mod<Vector<IntE>> vector : module) {
			System.out.println(vector);
		}
		System.out.println();
			}
	
	@Test
	void testZ2Z() {
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<IntE> syzygy1 = new ArrayList<>();
		syzygy1.add(z.getInteger(2));
		syzygy1.add(z.getInteger(0));
		syzygies.add(syzygy1);
		GenericPIDModule<IntE, Vector<IntE>> module = GenericPIDModule.fromSyzygies(z, 2, syzygies);
		for (Mod<Vector<IntE>> basis : module.getModuleGenerators()) {
			System.out.println(basis);
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				List<IntE> vector = new ArrayList<>();
				vector.add(z.getInteger(i));
				vector.add(z.getInteger(j));
				Vector<IntE> v = new Vector<>(vector);
				System.out.println("Not reduced: " + v);
				System.out.println("Reduced: " + module.reduce(v));
			}
		}
	}
}
