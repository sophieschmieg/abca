package varieties.curves;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField.FFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.EllipticCurve.TatesAlgorithmResult;
import varieties.curves.elliptic.EllipticCurve.WithComplexMultiplication;

class CMLiftTest {

	@Test
	void testCMCurvesSqrt23() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(23), q.zero(), q.one()));
		WithComplexMultiplication withComplexMultiplication = EllipticCurve.withComplexMultiplication(nf,
				nf.maximalOrder().getModuleGenerators().get(1));
		System.out.println(withComplexMultiplication.getCurve());
		NumberFieldIntegers order = ((NumberField) withComplexMultiplication.getCurve().getField()).maximalOrder();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(z.getInteger(40)) > 0) {
				break;
			}
			for (NumberFieldIdeal primeIdeal : order.idealsOver(prime)) {
				TatesAlgorithmResult<NFE, FFE> reduction = withComplexMultiplication.getCurve()
						.tatesAlgorithm(order.localizeAndQuotient(primeIdeal));
				System.out.println(reduction.getReducedScheme());
				System.out.println(reduction.hasGoodReduction());
			}
		}
	}

	@Test
	void testCMCurvesSqrt2() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(2), q.zero(), q.one()));
		WithComplexMultiplication withComplexMultiplication = EllipticCurve.withComplexMultiplication(nf,
				nf.maximalOrder().getModuleGenerators().get(1));
		NumberFieldIntegers order = ((NumberField) withComplexMultiplication.getCurve().getField()).maximalOrder();
		System.out.println(order.idealFactorization(order.getIdeal(withComplexMultiplication.getCurve().jInvariant())));
		TatesAlgorithmResult<NFE, FFE> reduction = withComplexMultiplication.getCurve().tatesAlgorithm(order.localizeAndQuotient(order.idealsOver(11).get(0)));
		System.out.println(reduction.getReducedCurve());
		System.out.println(withComplexMultiplication.getCurve());
	}

	@Test
	void testCMCurvesGauss() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));

		WithComplexMultiplication withComplexMultiplication = EllipticCurve.withComplexMultiplication(nf, nf.alpha());
		System.out.println(withComplexMultiplication.getCurve());
	}

	@Test
	void testCMCurves() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(5), q.zero(), q.one()));
		WithComplexMultiplication withComplexMultiplication = EllipticCurve.withComplexMultiplication(nf, nf.alpha());
		System.out.println(withComplexMultiplication.getCurve());
	}

	// @Test
	void testCMCurvesRay() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(36), q.zero(), q.one()));
		WithComplexMultiplication withComplexMultiplication = EllipticCurve.withComplexMultiplication(nf, nf.alpha());
		System.out.println(withComplexMultiplication.getCurve());
		System.out.println(withComplexMultiplication.getRayClassField());
		System.out.println(withComplexMultiplication.getCurve().getField());
		NumberField field = (NumberField) withComplexMultiplication.getCurve().getField();
		System.out.println(field.maximalOrder().idealsOver(5));
		TatesAlgorithmResult<NFE, FFE> tate = withComplexMultiplication.getCurve().tatesAlgorithm(
				field.maximalOrder().localizeAndQuotient(field.maximalOrder().idealsOver(5).get(0))
				);
		System.out.println(tate.getReducedCurve());
		System.out.println(tate.getReducedCurve().getField());
		System.out.println(tate.getReducedCurve().jInvariant());
	}
//		FiniteField fq = FiniteField.getFiniteField(83*83);
//		
//		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(8000));
//		System.out.println(curve);
//		System.out.println(curve.isSupersingular());
//		System.out.println(curve.getNumberOfElements());
//		System.out.println(curve.trace());
//		System.out.println(curve.getTorsionPoints(2));
//		Graph<EllipticCurve<FFE>, Isogeny<FFE>> isogenyGraph = curve.isogenyGraph(3);
//		for (EllipticCurve<FFE> vertex : isogenyGraph.vertexSet()) {
//			for (Edge<EllipticCurve<FFE>, Isogeny<FFE>> edge : isogenyGraph.edges(vertex)) {
//			System.out.println(edge.firstVertex().jInvariant()+"-->"+edge.secondVertex().jInvariant());
//			}
//		}

//		Complex c = Complex.c(40);
//		Map<ComplexNumber, Integer> sqrts2 = c.roots(c.getUnivariatePolynomialRing().getPolynomial(c.getInteger(5), c.zero(), c.one()));
//		ComplexNumber sqrt2 = null;
//		for (ComplexNumber t : sqrts2.keySet()) {
//			if (t.complexPart().compareTo(c.getField().zero()) > 0) {
//				sqrt2 = t;
//				break;
//			}
//		} // sqrt(-2) => 8000
//		//	 sqrt(-5) => 308736+.7i
//		System.out.println(EllipticCurve.fromTau(c, sqrt2).getCurve().jInvariant());
//	}

//	@Test
//	void testSaturation() {
//		FiniteField fq = FiniteField.getFiniteField(25);
//		PolynomialRing<FFE> p = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.LEX);
//		Polynomial<FFE> mod = p.subtract(p.getVarPower(2, 2), p.add(p.getVarPower(1, 3), p.getVar(1), p.one()));
//		PolynomialIdeal<FFE> ideal = p.getIdeal(Collections.singletonList(mod));
//		System.out.println("Ring: " + p);
//		System.out.println("Ideal: " + ideal);
//		CoordinateRing<FFE> ring = ideal.divideOut();
//		List<CoordinateRingElement<FFE>> list = new ArrayList<>();
//		list.add(ring.getEmbedding(p.getVar(1)));
//		list.add(ring.getEmbedding(p.subtract(p.getVar(2), p.one())));
//		LocalizedCoordinateRing<FFE> lcr = new LocalizedCoordinateRing<>(fq, ring, ring.getIdeal(list));
//		LocalizedElement<FFE> le1 = lcr.getEmbedding(p.subtract(p.getVar(2), p.one()));
//
//		assertEquals(Value.ONE, lcr.valuation(le1));
//		Polynomial<FFE> mod2 = p.subtract(
//				p.add(p.getVarPower(1, 3), p.getVarPower(2, 3), p.multiply(p.getVar(1), p.getVarPower(2, 2))),
//				p.getVar(2));
//		CoordinateRing<FFE> ring2 = p.getIdeal(Collections.singletonList(mod2)).divideOut();
//		List<CoordinateRingElement<FFE>> list2 = new ArrayList<>();
//		list2.add(ring.getEmbedding(p.getVar(1)));
//		list2.add(ring.getEmbedding(p.getVar(2)));
//		LocalizedCoordinateRing<FFE> lcr2 = new LocalizedCoordinateRing<>(fq, ring2, ring2.getIdeal(list2));
//		assertEquals(new Value(3), lcr2.valuation(lcr2.getEmbedding(p.getVar(2))));
//		assertEquals(new Value(12), lcr2.valuation(lcr2.getEmbedding(p.getVarPower(2, 4))));
//		LocalizedElement<FFE> testElement = lcr2.divide(lcr2.getEmbedding(p.getVarPower(2, 4)),
//				lcr2.getEmbedding(p.getVarPower(1, 12)));
//		assertEquals(Value.ZERO, lcr2.valuation(testElement));
//		assertEquals(fq.one(), lcr2.reduce(testElement));
//	}

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
