package fields.finitefields;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Set;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField.FFE;

class FiniteFieldTest {

	@Test
	void listAllElementsTest() {
		FiniteField f = FiniteField.getFiniteField(25);
		Set<FFE> elements = new TreeSet<>();
		for (FFE p : f) {
			elements.add(p);
		}
		assertEquals(25, elements.size());
	}

//	@Test
//	void discreteLogarithmTest() {
//		FiniteField field = FiniteField.getFiniteField(25);
//		field.discreteLogarithm(field.primitiveRoot(), field.getInteger(3));
//		Integers z = Integers.z();
//		List<IntE> smoothnessBase = new ArrayList<>();
//		Iterator<IntE> it = z.primes();
//		for (int i = 0; i < 10; i++) {
//			smoothnessBase.add(it.next());
//		}
//		SortedSet<IntE> smoothBaseSet = new TreeSet<>();
//		smoothBaseSet.addAll(smoothnessBase);
//		UnivariatePolynomialRing<IntE> polynomials = z.getUnivariatePolynomialRing();
//		UnivariatePolynomial<IntE> f = polynomials.getPolynomial(z.getInteger(-30), z.one(), z.one());
//		List<IntE> smooth = z.polynomialSieveSingleVariable(f, z.zero(), z.getInteger(500), smoothnessBase);
//		Set<IntE> smoothSet = new TreeSet<>();
//		smoothSet.addAll(smooth);
//		for (int i = 0; i < 500; i++) {
//			IntE eval = z.getUnivariatePolynomialRing().evaluate(f, z.getInteger(i));
//			if (!eval.equals(z.zero()) && smoothBaseSet.containsAll(z.uniqueFactorization(eval).primeFactors())) {
//				System.out.println(i + ": " + eval);
//				assertTrue(smoothSet.contains(z.getInteger(i)));
//			} else {
//				assertFalse(smoothSet.contains(z.getInteger(i)));
//			}
//		}
//		PolynomialRing<IntE> int2 = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
//		Polynomial<IntE> g = polynomials.homogenize(f);
//		List<Vector<IntE>> smoothPairs = z.polynomialSieve(g, z.getInteger(-500), z.getInteger(500), smoothnessBase);
//		Set<Vector<IntE>> smoothPairSet = new TreeSet<>();
//		smoothPairSet.addAll(smoothPairs);
//		for (int i = -500; i < 500; i++) {
//			for (int j = -500; j < 500; j++) {
//				IntE eval = int2.evaluate(g, z.getInteger(i), z.getInteger(j));
//				if (!eval.equals(z.zero()) && smoothBaseSet.containsAll(z.uniqueFactorization(eval).primeFactors())) {
//					System.out.println("(" + i + ", " + j + "): " + eval);
//					assertTrue(smoothPairSet.contains(new Vector<>(z.getInteger(i), z.getInteger(j))));
//				} else {
//					assertFalse(smoothPairSet.contains(new Vector<>(z.getInteger(i), z.getInteger(j))));
//				}
//			}
//		}
//	}
}
