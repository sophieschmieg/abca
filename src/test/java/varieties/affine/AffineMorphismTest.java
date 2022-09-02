package varieties.affine;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;

class AffineMorphismTest {

	@Test
	void testPreimage() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> one = q.getUnivariatePolynomialRing();
		PolynomialRing<Fraction> two = AbstractPolynomialRing.getPolynomialRing(q, 2, Monomial.GREVLEX);
		AffineScheme<Fraction> line = new AffineScheme<>(q, one.getZeroIdeal().divideOut());
		AffineScheme<Fraction> plane = new AffineScheme<>(q, two.getZeroIdeal().divideOut());
		List<Polynomial<Fraction>> nodalMapList = new ArrayList<>();
		nodalMapList.add(one.subtract(one.getVarPower(2), one.one()));
		nodalMapList.add(one.multiply(nodalMapList.get(0), one.getVar()));
		AffineMorphism<Fraction> morphism = AffineMorphism.fromPolynomials(line, plane, nodalMapList);
		System.out.println(morphism);
		System.out.println(morphism.getDomain());
		System.out.println(morphism.getRange());
		assertEquals(new AffinePoint<>(q, q.zero(), q.zero()), morphism.evaluate(new AffinePoint<>(q, q.one())));
		System.out.println(morphism.evaluate(new AffinePoint<>(q, q.getInteger(2))));
		System.out.println(morphism.preimage(new AffinePoint<>(q, q.zero(), q.zero())).getDomain());
		System.out.println(morphism.preimage(new AffinePoint<>(q, q.zero(), q.zero())).closure().getDomain());
		System.out.println(morphism.preimage(new AffinePoint<>(q, q.getInteger(3), q.getInteger(6))).getDomain());
		List<CoordinateRingElement<Fraction>> pointList = new ArrayList<>();
		pointList.add(plane.getCoordinateRing().getVar(1));
		pointList.add(plane.getCoordinateRing().getVar(2));
		AffineMorphism<Fraction> embedding = AffineMorphism.getClosedImmersion(plane,
				plane.getCoordinateRing().getIdeal(pointList));
		AffineMorphism<Fraction> pre = morphism.preimage(embedding);
		System.out.println(pre.getDomain());
		System.out.println(pre.closure().getDomain());
		System.out.println(morphism.closure().getDomain());
		List<Polynomial<Fraction>> lineList = new ArrayList<>();
		lineList.add(one.getVar());
		lineList.add(one.multiply(2, one.getVar()));// subtract(one.getVar(), one.one()));
		AffineMorphism<Fraction> lineMorphism = AffineMorphism.fromPolynomials(line, plane, lineList);
		System.out.println(lineMorphism.closure().getDomain());
		pre = morphism.preimage(lineMorphism);
		System.out.println(pre.getDomain());
		System.out.println(pre.closure().getDomain());
	}

	@Test
	void immersionTest() throws IOException {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> one = q.getUnivariatePolynomialRing();
		PolynomialRing<Fraction> two = AbstractPolynomialRing.getPolynomialRing(q, 2, Monomial.GREVLEX);
		AffineScheme<Fraction> nodal = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.parse("-1*Y^2 + X^3 + X^2"))).divideOut());
		AffineScheme<Fraction> line = new AffineScheme<>(q, one.getZeroIdeal().divideOut());
		List<Polynomial<Fraction>> normalization = new ArrayList<>();
		normalization.add(one.getPolynomial(q.getInteger(-1), q.zero(), q.one()));
		normalization.add(one.getPolynomial(q.zero(), q.getInteger(-1), q.zero(), q.one()));
		AffineMorphism<Fraction> resolve = AffineMorphism.fromPolynomials(line, nodal, normalization);
		resolve.closure();
		AffineScheme<Fraction> cross = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.parse("X*Y"))).divideOut());
		List<Polynomial<Fraction>> ltc = new ArrayList<>();
		ltc.add(one.getVar());
		ltc.add(one.zero());
		AffineMorphism<Fraction> lineIntoCross = AffineMorphism.fromPolynomials(line, cross, ltc);
		assertTrue(lineIntoCross.isClosedImmersion());
		assertFalse(lineIntoCross.isOpenImmersion());
		assertTrue(lineIntoCross.isBirational());
		lineIntoCross.birationalInverse();
		AffineMorphism<Fraction> crossWithoutLine = AffineMorphism.getOpenImmersion(cross, two.getVar(2))
				.getOpenImmersion();
		assertTrue(crossWithoutLine.isOpenImmersion());
		assertTrue(crossWithoutLine.isBirational());
		crossWithoutLine.birationalInverse();
		AffineScheme<Fraction> plane = new AffineScheme<>(q, two.getZeroIdeal().divideOut());
		List<Polynomial<Fraction>> nodalMapList = new ArrayList<>();
		nodalMapList.add(two.getVar(1));
		nodalMapList.add(two.getVar(2));
		AffineMorphism<Fraction> morphism = AffineMorphism.fromPolynomials(nodal, plane, nodalMapList);
		morphism.graph().getProduct().simplify();
		assertTrue(morphism.isSeparated());
		assertTrue(morphism.isImmersion());
		assertTrue(morphism.isClosedImmersion());
		assertTrue(morphism.isFinite());
		AffineScheme<Fraction> onePointMissing = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.parse("X*Y + 5*Y + -1"))).divideOut());
		List<Polynomial<Fraction>> embedding = new ArrayList<>();
		embedding.add(two.getVar(1));
		AffineMorphism<Fraction> openImmersion = AffineMorphism.fromPolynomials(onePointMissing, line, embedding);
		// AffineScheme.restrictAwayFrom(openImmersion);
		assertTrue(openImmersion.isImmersion());
		assertFalse(openImmersion.isClosedImmersion());
		assertTrue(openImmersion.isOpenImmersion());
		assertFalse(openImmersion.isFinite());
		List<Polynomial<Fraction>> messedUp = new ArrayList<>();
		messedUp.add(two.subtract(two.getVar(1), two.getInteger(2)));
		messedUp.add(two.one());
//		AffineMorphism<Fraction> messedUpImmersion = AffineMorphism.fromPolynomials(onePointMissing, plane, messedUp);
		// AffineScheme.restrictAwayFrom(messedUpImmersion);
		AffineMorphism<Fraction> twoCover = AffineMorphism.fromPolynomials(nodal, line,
				Collections.singletonList(two.getVar(1)));
		assertFalse(twoCover.isImmersion());
		assertFalse(twoCover.isClosedImmersion());
		assertFalse(twoCover.isOpenImmersion());
		assertTrue(twoCover.isFinite());
		AffineScheme<Fraction> twoPointsMissing = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.parse("X^2*Y + -1*Y + -1"))).divideOut());
		List<Polynomial<Fraction>> map = new ArrayList<>();
		map.add(two.subtract(two.getVar(1), two.getInteger(4)));
		map.add(two.multiply(two.getVar(2), two.subtract(two.getVar(1), two.getInteger(1))));
		AffineMorphism<Fraction> twoToOne = AffineMorphism.fromPolynomials(twoPointsMissing, onePointMissing, map);
		assertTrue(twoToOne.isOpenImmersion());// AffineScheme.restrictAwayFrom(twoToOne);
		one.setVariableName("X");
		AffineMorphism<Fraction> translation = AffineMorphism.fromPolynomials(line, line,
				Collections.singletonList(one.parse("X + 1")));
		assertTrue(translation.isImmersion());
		assertTrue(translation.isClosedImmersion());
		assertTrue(translation.isOpenImmersion());
		assertTrue(translation.isFinite());
		System.out.println(translation.inverse());
		AffineScheme<Fraction> lineSuperfluousX = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.getVar(1))).divideOut());
		lineSuperfluousX.simplify();
		AffineScheme<Fraction> lineSuperfluousY = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.getVar(2))).divideOut());
		lineSuperfluousY.simplify();
		ArrayList<Polynomial<Fraction>> yList = new ArrayList<>();
		yList.add(one.zero());
		yList.add(one.getVar());
		ArrayList<Polynomial<Fraction>> xList = new ArrayList<>();
		xList.add(one.getVar());
		xList.add(one.zero());
		AffineMorphism<Fraction> lineToLineSupX = AffineMorphism.fromPolynomials(line, lineSuperfluousX, yList);
		AffineMorphism<Fraction> lineToLineSupY = AffineMorphism.fromPolynomials(line, lineSuperfluousY, xList);
		AffineMorphism<Fraction> lineSupXToLine = AffineMorphism.fromPolynomials(lineSuperfluousX, line,
				Collections.singletonList(two.getVar(2)));
		AffineMorphism<Fraction> lineSupYToLine = AffineMorphism.fromPolynomials(lineSuperfluousY, line,
				Collections.singletonList(two.getVar(1)));
		assertTrue(lineToLineSupX.isImmersion());
		assertTrue(lineToLineSupX.isClosedImmersion());
		assertTrue(lineToLineSupX.isFinite());
		assertTrue(lineToLineSupY.isImmersion());
		assertTrue(lineToLineSupY.isClosedImmersion());
		assertTrue(lineToLineSupY.isFinite());
		assertTrue(lineSupXToLine.isImmersion());
		assertTrue(lineSupXToLine.isClosedImmersion());
		assertTrue(lineSupXToLine.isFinite());
		assertTrue(lineSupYToLine.isImmersion());
		assertTrue(lineSupYToLine.isClosedImmersion());
		assertTrue(lineSupYToLine.isFinite());
	}

}
