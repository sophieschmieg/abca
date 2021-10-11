package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.local.Value;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;

class CMLiftTest {

	@Test
	void testSaturation() {
		FiniteField fq = FiniteField.getFiniteField(25);
		PolynomialRing<FFE> p = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.LEX);
		Polynomial<FFE> mod = p.subtract(p.getVarPower(2, 2), p.add(p.getVarPower(1, 3), p.getVar(1), p.one()));
		PolynomialIdeal<FFE> ideal = p.getIdeal(Collections.singletonList(mod));
		System.out.println("Ring: " + p);
		System.out.println("Ideal: " + ideal);
		CoordinateRing<FFE> ring = ideal.divideOut();
		List<CoordinateRingElement<FFE>> list = new ArrayList<>();
		list.add(ring.getEmbedding(p.getVar(1)));
		list.add(ring.getEmbedding(p.subtract(p.getVar(2), p.one())));
		LocalizedCoordinateRing<FFE> lcr = new LocalizedCoordinateRing<>(fq, ring, ring.getIdeal(list));
		LocalizedElement<FFE> le1 = lcr.getEmbedding(p.subtract(p.getVar(2), p.one()));

		assertEquals(Value.ONE, lcr.valuation(le1));
		Polynomial<FFE> mod2 = p.subtract(
				p.add(p.getVarPower(1, 3), p.getVarPower(2, 3), p.multiply(p.getVar(1), p.getVarPower(2, 2))),
				p.getVar(2));
		CoordinateRing<FFE> ring2 = p.getIdeal(Collections.singletonList(mod2)).divideOut();
		List<CoordinateRingElement<FFE>> list2 = new ArrayList<>();
		list2.add(ring.getEmbedding(p.getVar(1)));
		list2.add(ring.getEmbedding(p.getVar(2)));
		LocalizedCoordinateRing<FFE> lcr2 = new LocalizedCoordinateRing<>(fq, ring2, ring2.getIdeal(list2));
		assertEquals(new Value(3), lcr2.valuation(lcr2.getEmbedding(p.getVar(2))));
		assertEquals(new Value(12), lcr2.valuation(lcr2.getEmbedding(p.getVarPower(2, 4))));
		LocalizedElement<FFE> testElement = lcr2.divide(lcr2.getEmbedding(p.getVarPower(2, 4)),
				lcr2.getEmbedding(p.getVarPower(1, 12)));
		assertEquals(Value.ZERO, lcr2.valuation(testElement));
		assertEquals(fq.one(), lcr2.reduce(testElement));
	}

//	@Test
//	void testMultiplicationMorphism() {
//		FiniteField fq = FiniteField.getFiniteField(25);
//		EllipticCurve<FFE> c = new EllipticCurve<>(fq, fq.one(), fq.one());
//		for (int n = -5; n <= 5; n++) {
//			System.out.println(n);
//			Morphism<FFE> m = c.multiplicationMorphism(n);
//			for (ProjectivePoint<FFE> p : c) {
//				System.out.println(p);
//				assertEquals(c.multiply(n, p), m.evaluate(p));
//			}
//		}
//	}

//	@Test
//	void testKernelPointMorphism() {
//		FiniteField fq = FiniteField.getFiniteField(25);
//		EllipticCurve<FFE> c = new EllipticCurve<>(fq, fq.one(), fq.one());
//		for (ProjectivePoint<FFE> point : c) {
//			System.out.println(point);
//			Isogeny<FFE> m = new KernelPointIsogeny<>(c, point, c.getOrder(point).intValueExact());
//			for (ProjectivePoint<FFE> p : c) {
//				System.out.println(p);
//				assertEquals(m.evaluate(p), m.asMorphism().evaluate(p));
//			}
//		}
//	}

//	@Test
//	void testTranslationMorphism() {
//		FiniteField fq = FiniteField.getFiniteField(25);
//		EllipticCurve<FFE> c = new EllipticCurve<>(fq, fq.one(), fq.one());
//		for (ProjectivePoint<FFE> point : c) {
//			if (point.equals(c.neutral())) {
//				continue;
//			}
//			System.out.println(point);
//			Morphism<FFE> m = c.translationMorphism(point);
//			for (ProjectivePoint<FFE> p : c) {
//				System.out.println(p);
//				assertEquals(c.add(p, point), m.evaluate(p));
//			}
//		}
//	}
//
//	@Test
//	void test5() {
//		// PAdicField qp = new PAdicField(BigInteger.valueOf(23), 5);
//		FiniteField fq = FiniteField.getFiniteField(625);
//		UnivariatePolynomialRing<FFE> polynomials = fq.getUnivariatePolynomialRing();
//		UnivariatePolynomial<FFE> rhs = polynomials.getPolynomial(fq.zero(), fq.one(), fq.zero(), fq.one());
//		System.out.println(rhs);
//		UnivariatePolynomial<FFE> dx = polynomials.getPolynomial(fq.one(), fq.zero(), fq.getInteger(3));
//		System.out.println(dx);
//		ExtendedEuclideanResult<Polynomial<FFE>> ee = polynomials.extendedEuclidean(polynomials.multiply(2, rhs), dx);
//		System.out.println("2*Y*Y*(" + ee.getCoeff2() + ") + (" + dx + ")*(" + ee.getCoeff1() + ")=" + ee.getGcd());
//	}
//
//	@Test
//	void test25() {
//		// PAdicField qp = new PAdicField(BigInteger.valueOf(23), 5);
//		ModularIntegerRing fq = new ModularIntegerRing(25);
//		UnivariatePolynomialRing<ModularIntegerRingElement> polynomials = fq.getUnivariatePolynomialRing();
//		UnivariatePolynomial<ModularIntegerRingElement> rhs = polynomials.getPolynomial(fq.zero(), fq.one(), fq.zero(),
//				fq.one());
//		UnivariatePolynomial<ModularIntegerRingElement> dx = polynomials.getPolynomial(fq.one(), fq.zero(),
//				fq.getInteger(3));
//		System.out.println("Y^2 = " + rhs);
//		System.out.println("3X^2+a = " + dx);
//		System.out
//				.println(
//						polynomials.divideScalar(
//								polynomials.subtract(polynomials.power(rhs, 5),
//										polynomials.substitute(rhs,
//												Collections.singletonList(polynomials.getVarPower(5)))),
//								fq.getInteger(5)));
//		// rhs = polynomials.toUnivariate(polynomials.power(rhs, 5));
//		System.out.println("Y^10 = " + rhs);
//		// dx = polynomials.substitute(dx,
//		// Collections.singletonList(polynomials.getVarPower(5)));
//		// System.out.println("3(X^5)^2+a = " + dx);
//		ExtendedResultantResult<ModularIntegerRingElement> ee = polynomials
//				.extendedResultant(polynomials.multiply(2, rhs), dx);
//		ModularIntegerRingElement inverse = fq.inverse(ee.getResultant());
//		System.out.println("2*Y*Y*(" + polynomials.multiply(inverse, ee.getCoeff2()) + ") + (" + dx + ")*("
//				+ polynomials.multiply(inverse, ee.getCoeff1()) + ")=" + polynomials.multiply(inverse, ee.getGcd()));
//		PolynomialRing<ModularIntegerRingElement> p2 = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.REVLEX);
//		Polynomial<ModularIntegerRingElement> mod = p2.subtract(p2.getVarPower(2, 2),
//				p2.add(p2.getVarPower(1, 3), p2.getVar(1)));
//		System.out.println(mod);
//		CoordinateRing<ModularIntegerRingElement> cr = new CoordinateRing<>(p2,
//				p2.getIdeal(Collections.singletonList(mod)));
//		System.out.println(cr);
//		CoordinateRingElement<ModularIntegerRingElement> phiX = cr.getEmbedding(p2.getVarPower(1, 5));
//		CoordinateRingElement<ModularIntegerRingElement> phiY = cr.getEmbedding(p2.getVarPower(2, 5));
//		System.out.println("X^5=" + phiX);
//		System.out.println("Y^5=" + phiY);
//		System.out.println(cr.subtract(cr.power(phiY, 2), cr.add(cr.power(phiX, 3), phiX)));
//		System.out.println();
//		System.out.println();
//	}

}
