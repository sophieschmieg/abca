package util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.vectors.Vector;

class PellsEquationTest {

	private int pell(Vector<IntE> solution, int d) {
		Integers z = Integers.z();
		IntE x = solution.get(1);
		IntE y = solution.get(2);
		int result = z.subtract(z.multiply(x, x), z.multiply(d, y, y)).intValueExact();
		System.out.println(x + "^2 - " + d + "*" + y + "^2 = " + result);
		return result;
	}

	@Test
	void testPellsEquation() {
		Integers z = Integers.z();
		List<Vector<IntE>> result = MiscAlgorithms.pellsEquation(z.getInteger(5), z.getInteger(-1), false);
		pell(result.get(0), 5);
		assertEquals(z.getInteger(2), result.get(0).get(1));
		result = MiscAlgorithms.pellsEquation(z.getInteger(7), z.getInteger(1), false);
		pell(result.get(0), 7);
		result = MiscAlgorithms.pellsEquation(z.getInteger(5), z.getInteger(-4), false);
		pell(result.get(0), 5);
		result = MiscAlgorithms.pellsEquation(z.getInteger(5), z.getInteger(-8), true);
		assertTrue(result.isEmpty());
		result = MiscAlgorithms.pellsEquation(z.getInteger(5), z.getInteger(-16), true);
		assertEquals(-16, pell(result.get(0), 5));
		result = MiscAlgorithms.pellsEquation(z.getInteger(257), z.getInteger(-4), false);
		if (!result.isEmpty()) {
			pell(result.get(0), 257);
		}
		result = MiscAlgorithms.pellsEquation(z.getInteger(257), z.getInteger(4), false);
		if (!result.isEmpty()) {
			pell(result.get(0), 257);
		}
		// assertEquals(-4, result.get());
		result = MiscAlgorithms.pellsEquation(z.getInteger(257), z.getInteger(32), true);
		pell(result.get(0), 257);
		// assertEquals(-32, result.get());
	}

	@Test
	void testNumberField() {
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> p = z.getUnivariatePolynomialRing();
		List<Vector<IntE>> result = MiscAlgorithms.pellsEquation(z.getInteger(7), z.getInteger(1), false);
		NumberField nf = NumberField.getNumberFieldFromIntegerPolynomial(p.subtract(p.getVarPower(2), p.getInteger(7)));
		System.out.println(nf.classNumber());
	}
}
