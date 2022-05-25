package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;

class CoordinateRingTest {

	@Test
	void testIdeal() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.REVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(
				polynomialRing.subtract(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1)),
				polynomialRing.getVarPower(2, 2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> cr = ideal.divideOut();
		CoordinateIdeal<PFE> crIdeal = cr.getIdeal(Collections.singletonList(cr.multiply(cr.getVar(1), cr.getVar(2))));
		System.out.println(crIdeal);
		List<CoordinateRingElement<PFE>> g2 = new ArrayList<>();
		g2.add(cr.getVar(1));
		g2.add(cr.getVar(2));
		CoordinateIdeal<PFE> crIdeal2 = cr.getIdeal(g2);
		System.out.println(crIdeal2);
		System.out.println(cr.multiply(crIdeal2, crIdeal2));
		System.out.println(cr.divideChecked(cr.subtract(cr.power(cr.getVar(1), 4), cr.power(cr.getVar(1), 2)), cr.multiply(cr.getVar(1), cr.getVar(2))));
	}

	@Test
	void testInverseCbrt2() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(2)));
		generators.add(
				polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)),
						polynomialRing.getInteger(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2),
				polynomialRing.multiply(2, polynomialRing.getVar(1))));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = ideal.divideOut();
		assertEquals(3, coordinateRing.degree());
		assertEquals(
				coordinateRing.getEmbedding(polynomialRing.divideScalar(polynomialRing.getVar(2), fp.getElement(2))),
				coordinateRing.inverse(coordinateRing.getVar(1)));
		assertEquals(coordinateRing.getVar(1),
				coordinateRing.divideChecked(coordinateRing.getVar(2), coordinateRing.getVar(1)));
	}

	@Test
	void testInverseImaginaryTwice() {
		PrimeField fp = PrimeField.getPrimeField(11);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 2), polynomialRing.one()));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(2, 2), polynomialRing.one()));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = ideal.divideOut();
		assertEquals(4, coordinateRing.degree());
		CoordinateRingElement<PFE> x = coordinateRing.getVar(1);
		CoordinateRingElement<PFE> y = coordinateRing.getVar(2);
		assertTrue(coordinateRing.isUnit(x));
		assertTrue(coordinateRing.isUnit(y));
		assertEquals(coordinateRing.negative(x), coordinateRing.inverse(x));
		assertEquals(coordinateRing.negative(y), coordinateRing.inverse(y));
		assertEquals(x, coordinateRing.inverse(coordinateRing.negative(x)));
		assertEquals(y, coordinateRing.inverse(coordinateRing.negative(y)));
		CoordinateRingElement<PFE> xpy = coordinateRing.add(x, y);
		CoordinateRingElement<PFE> xmy = coordinateRing.subtract(x, y);
		assertEquals(coordinateRing.zero(), coordinateRing.multiply(xpy, xmy));
		assertFalse(coordinateRing.isUnit(xpy));
		assertFalse(coordinateRing.isUnit(xmy));
		CoordinateRingElement<PFE> xpy2 = coordinateRing.multiply(xpy, xpy);
		assertTrue(coordinateRing.isDivisible(xpy2, xpy));
		CoordinateRingElement<PFE> xpy2byXpy = coordinateRing.divideChecked(xpy2, xpy);
		assertEquals(xpy2, coordinateRing.multiply(xpy2byXpy, xpy));
		CoordinateRingElement<PFE> xmy2 = coordinateRing.multiply(xmy, xmy);
		assertTrue(coordinateRing.isDivisible(xmy2, xmy));
		CoordinateRingElement<PFE> xmy2byXmy = coordinateRing.divideChecked(xmy2, xmy);
		assertEquals(xmy2, coordinateRing.multiply(xmy2byXmy, xmy));
	}

	@Test
	void testDegree3() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(1)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2), polynomialRing.getVar(2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		assertEquals(3, coordinateRing.degree());
		System.out.println(coordinateRing.hilbertPolynomial());
	}

	@Test
	void testDegreeEllipticIntersection() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(
				polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1), polynomialRing.one()),
				polynomialRing.getVarPower(2, 2)));
		generators.add(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		assertEquals(5, coordinateRing.degree());
		System.out.println(coordinateRing.hilbertPolynomial());
	}

	@Test
	void testDegreeElliptic() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(
				polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1), polynomialRing.one()),
				polynomialRing.getVarPower(2, 2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		assertEquals(3, coordinateRing.degree());
		assertEquals(1, coordinateRing.genus());
	}

	@Test
	void testDegreeCommon() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 4, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(
				polynomialRing.subtract(polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1)),
						polynomialRing.getVarPower(2, 2)));
		generators.add(
				polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(4)),
						polynomialRing.multiply(polynomialRing.getVar(2), polynomialRing.getVar(3))));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 4), polynomialRing.getVarPower(2, 4),
				polynomialRing.getVarPower(3, 2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		int degree = coordinateRing.degree();
		generators.add(polynomialRing.getVar(4));
		PolynomialIdeal<PFE> ideal2 = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing2 = new CoordinateRing<>(polynomialRing, ideal2);
		assertEquals(degree, coordinateRing2.degree());
		System.out.println(coordinateRing.hilbertPolynomial());
	}

	@Test
	void testDegreeDimension3Common() {
		PrimeField fp = PrimeField.getPrimeField(13);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 10, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(
				polynomialRing.subtract(polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1)),
						polynomialRing.getVarPower(2, 2)));
		generators.add(
				polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(4)),
						polynomialRing.multiply(polynomialRing.getVar(2), polynomialRing.getVar(3))));
		generators.add(polynomialRing.add(polynomialRing.getVarPower(1, 4), polynomialRing.getVarPower(2, 4),
				polynomialRing.getVarPower(3, 2)));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		int degree = coordinateRing.degree();
		generators.add(polynomialRing.getVar(4));
		PolynomialIdeal<PFE> ideal2 = polynomialRing.getIdeal(generators);
		CoordinateRing<PFE> coordinateRing2 = new CoordinateRing<>(polynomialRing, ideal2);
		assertEquals(degree, coordinateRing2.degree());
		System.out.println(coordinateRing.hilbertPolynomial());
	}

}
