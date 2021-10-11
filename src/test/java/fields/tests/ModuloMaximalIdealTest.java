package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.ModuloMaximalIdealResult;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;

class ModuloMaximalIdealTest {

	@Test
	void integersTest() {
		Integers z = Integers.z();
		PolynomialRing<IntE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
		List<Polynomial<IntE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(2)));
		generators.add(
				polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)),
						polynomialRing.getInteger(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2),
				polynomialRing.multiply(2, polynomialRing.getVar(1))));
		generators.add(polynomialRing.getInteger(13));
		PolynomialIdeal<IntE> ideal = polynomialRing.getIdeal(generators);
		@SuppressWarnings("unchecked")
		ModuloMaximalIdealResult<Polynomial<IntE>, FFE> mod = (ModuloMaximalIdealResult<Polynomial<IntE>, FFE>) polynomialRing
				.moduloMaximalIdeal(ideal);
		FiniteField fq = (FiniteField) mod.getField();
		assertEquals(BigInteger.valueOf(13), fq.characteristic());
		assertEquals(BigInteger.valueOf(13 * 13 * 13), fq.getNumberOfElements());
		assertEquals(fq.getInteger(2), fq.power(mod.getReduction().evaluate(polynomialRing.getVar(1)), 3));
		assertEquals(fq.getInteger(4), fq.power(mod.getReduction().evaluate(polynomialRing.getVar(2)), 3));
		assertTrue(ideal.contains(polynomialRing.subtract(
				polynomialRing.power(mod.getLift().evaluate(mod.getReduction().evaluate(polynomialRing.getVar(1))), 3),
				polynomialRing.getInteger(2))));
	}
}
