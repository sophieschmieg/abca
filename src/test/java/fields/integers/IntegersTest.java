package fields.integers;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.integers.Integers.IntE;
import fields.integers.Integers.SmallestIntegerSolutionPreparation;
import fields.vectors.Vector;

class IntegersTest {

	@Test
	void testFactors() {
		Integers z = Integers.z();
		SortedSet<IntE> factors = new TreeSet<>();
		factors.addAll(z.factors(z.getInteger(60)));
		assertEquals(12, factors.size());
		assertTrue(factors.contains(z.getInteger(1)));
		assertTrue(factors.contains(z.getInteger(2)));
		assertTrue(factors.contains(z.getInteger(3)));
		assertTrue(factors.contains(z.getInteger(4)));
		assertTrue(factors.contains(z.getInteger(5)));
		assertTrue(factors.contains(z.getInteger(6)));
		assertTrue(factors.contains(z.getInteger(10)));
		assertTrue(factors.contains(z.getInteger(12)));
		assertTrue(factors.contains(z.getInteger(15)));
		assertTrue(factors.contains(z.getInteger(20)));
		assertTrue(factors.contains(z.getInteger(30)));
		assertTrue(factors.contains(z.getInteger(60)));
	}

	@Test
	void testSmallestIntegerSolutionOneDimension() {
		Integers z = Integers.z();
		Vector<IntE> two = new Vector<>(z.getInteger(2));
		Vector<IntE> three = new Vector<>(z.getInteger(3));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(two);
		generators.add(three);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators);
		for (int i = -9; i <= 12; i++) {
			Vector<IntE> solution = z.smallestIntegerSolution(new Vector<>(z.getInteger(i)), preparation);
			System.out.println(i + " = " + solution.get(1) + "*2 + " + solution.get(2) + "*3");
		}
	}

	@Test
	void testSmallestIntegerSolutionTwoDimensionTrivial() {
		Integers z = Integers.z();
		Vector<IntE> a = new Vector<>(z.getInteger(2), z.getInteger(3));
		Vector<IntE> b = new Vector<>(z.getInteger(3), z.getInteger(5));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(a);
		generators.add(b);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators);
		for (int i = -9; i <= 12; i++) {
			for (int j = -9; j <= 12; j++) {
				Vector<IntE> target = new Vector<>(z.getInteger(i), z.getInteger(j));
				Vector<IntE> solution = z.smallestIntegerSolution(target, preparation);
				System.out.println(target + " = " + solution.get(1) + "*" + a + " + " + solution.get(2) + "*" + b);
			}
		}
	}

	@Test
	void testSmallestIntegerSolutionTwoDimensionNonTrivial() {
		Integers z = Integers.z();
		Vector<IntE> a = new Vector<>(z.getInteger(2), z.getInteger(3));
		Vector<IntE> b = new Vector<>(z.getInteger(3), z.getInteger(5));
		Vector<IntE> c = new Vector<>(z.getInteger(7), z.getInteger(9));
		Vector<IntE> d = new Vector<>(z.getInteger(5), z.getInteger(-3));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(a);
		generators.add(b);
		generators.add(c);
		generators.add(d);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators);
		for (int i = -9; i <= 12; i++) {
			for (int j = -9; j <= 12; j++) {
				Vector<IntE> target = new Vector<>(z.getInteger(i), z.getInteger(j));
				Vector<IntE> solution = z.smallestIntegerSolution(target, preparation);
				System.out.println(target + " = " + solution.get(1) + "*" + a + " + " + solution.get(2) + "*" + b
						+ " + " + solution.get(3) + "*" + c + " + " + solution.get(4) + "*" + d);
			}
		}
	}

	@Test
	void testSmallestModSolutionOneDimension() {
		Integers z = Integers.z();
		Vector<IntE> two = new Vector<>(z.getInteger(2));
		Vector<IntE> three = new Vector<>(z.getInteger(3));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(two);
		generators.add(three);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators, z.getInteger(5));
		for (int i = -9; i <= 12; i++) {
			Vector<IntE> solution = z.smallestIntegerSolution(new Vector<>(z.getInteger(i)), preparation);
			System.out.println(i + " = " + solution.get(1) + "*2 + " + solution.get(2) + "*3");
		}
	}

	@Test
	void testSmallestModSolutionTwoDimensionTrivial() {
		Integers z = Integers.z();
		Vector<IntE> a = new Vector<>(z.getInteger(2), z.getInteger(3));
		Vector<IntE> b = new Vector<>(z.getInteger(3), z.getInteger(5));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(a);
		generators.add(b);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators, z.getInteger(5));
		for (int i = -9; i <= 12; i++) {
			for (int j = -9; j <= 12; j++) {
				Vector<IntE> target = new Vector<>(z.getInteger(i), z.getInteger(j));
				Vector<IntE> solution = z.smallestIntegerSolution(target, preparation);
				System.out.println(target + " = " + solution.get(1) + "*" + a + " + " + solution.get(2) + "*" + b);
			}
		}
	}

	@Test
	void testSmallestModSolutionTwoDimensionNonTrivial() {
		Integers z = Integers.z();
		Vector<IntE> a = new Vector<>(z.getInteger(2), z.getInteger(3));
		Vector<IntE> b = new Vector<>(z.getInteger(3), z.getInteger(5));
		Vector<IntE> c = new Vector<>(z.getInteger(7), z.getInteger(9));
		Vector<IntE> d = new Vector<>(z.getInteger(5), z.getInteger(-3));
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.add(a);
		generators.add(b);
		generators.add(c);
		generators.add(d);
		SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(generators, z.getInteger(5));
		for (int i = -9; i <= 12; i++) {
			for (int j = -9; j <= 12; j++) {
				Vector<IntE> target = new Vector<>(z.getInteger(i), z.getInteger(j));
				Vector<IntE> solution = z.smallestIntegerSolution(target, preparation);
				System.out.println(target + " = " + solution.get(1) + "*" + a + " + " + solution.get(2) + "*" + b
						+ " + " + solution.get(3) + "*" + c + " + " + solution.get(4) + "*" + d);
			}
		}
	}
}
