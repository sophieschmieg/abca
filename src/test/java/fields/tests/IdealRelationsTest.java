package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Vector;

class IdealRelationsTest {

	private <T extends Element<T>> void testIdeal(Ring<T> ring, Ideal<T> ideal) {
		List<T> generators = ideal.generators();
		List<Vector<T>> relations = ideal.getSyzygies();
		for (Vector<T> relation : relations) {
			assertEquals(generators.size(), relation.dimension());
			T result = ring.zero();
			for (int i = 0; i < generators.size(); i++) {
				result = ring.add(ring.multiply(relation.get(i + 1), generators.get(i)), result);
			}
			assertEquals(ring.zero(), result);
		}
	}

	@Test
	void polynomialIdealRelationsTest() {
		PrimeField fp = PrimeField.getPrimeField(65537);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getVar(2)));
		generators.add(
				polynomialRing.subtract(polynomialRing.multiply(polynomialRing.getVar(1), polynomialRing.getVar(2)),
						polynomialRing.getInteger(2)));
		generators.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 2),
				polynomialRing.multiply(2, polynomialRing.getVar(1))));
		PolynomialIdeal<PFE> ideal = polynomialRing.getIdeal(generators);
		testIdeal(polynomialRing, ideal);
	}

	@Test
	void coordinateIdealRelationsTest() {
		PrimeField fp = PrimeField.getPrimeField(65537);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		Polynomial<PFE> curve = polynomialRing.subtract(
				polynomialRing.add(polynomialRing.getVarPower(1, 3), polynomialRing.getVar(1), polynomialRing.one()),
				polynomialRing.getVarPower(2, 2));
		CoordinateRing<PFE> coordinateRing = polynomialRing.getIdeal(Collections.singletonList(curve)).divideOut();
		List<CoordinateRingElement<PFE>> generators = new ArrayList<>();
		generators.add(coordinateRing.getVar(1));
		generators.add(coordinateRing.subtract(coordinateRing.getVar(2), coordinateRing.one()));
		CoordinateIdeal<PFE> ideal = coordinateRing.getIdeal(generators);
		testIdeal(coordinateRing, ideal);
	}

	@Test
	void numberFieldIdealSqrtMinus5RelationsTest() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(5), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NFE> generators = new ArrayList<>();
		generators.add(order.getInteger(2));
		generators.add(order.add(order.one(), nf.alpha()));
		NumberFieldIdeal ideal = order.getIdeal(generators);
		testIdeal(order, ideal);
	}

	@Test
	void numberFieldIdealDegree3RelationsTest() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				rationalPolynomialRing.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NFE> generators = new ArrayList<>();
		generators.add(order.getInteger(2));
		generators.add(order.add(order.one(), nf.alpha()));
		NumberFieldIdeal ideal = order.getIdeal(generators);
		testIdeal(order, ideal);
	}

	@Test
	void numberFieldIdealSqrt5RelationsTest() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(-5), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NFE> generators = new ArrayList<>();
		generators.add(order.getInteger(2));
		generators.add(order.add(order.one(), nf.alpha()));
		NumberFieldIdeal ideal = order.getIdeal(generators);
		testIdeal(order, ideal);
	}

	@Test
	void numberFieldIdealGaussianIntegersRelationsTest() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NFE> generators = new ArrayList<>();
		generators.add(order.getInteger(2));
		generators.add(order.add(order.one(), nf.alpha()));
		NumberFieldIdeal ideal = order.getIdeal(generators);
		testIdeal(order, ideal);
	}

}
