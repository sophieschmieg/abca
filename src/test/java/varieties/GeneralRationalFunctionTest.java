package varieties;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.affine.AffineMorphism;
import varieties.affine.AffineScheme;
import varieties.curves.ProjectiveLine;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

class GeneralRationalFunctionTest {

	@Test
	void evaluateLineTest() throws IOException {
		Rationals q = Rationals.q();
		ProjectiveLine<Fraction> line = new ProjectiveLine<>(q);
		AffineScheme<Fraction> affineLine = line.getAffineCover().get(0);
		PolynomialRing<Fraction> two = AbstractPolynomialRing.getPolynomialRing(q, 2, Monomial.GREVLEX);
		AffineScheme<Fraction> onePointMissing = new AffineScheme<>(q,
				two.getIdeal(Collections.singletonList(two.parse("X*Y + -1"))).divideOut());
		List<Polynomial<Fraction>> embedding = new ArrayList<>();
		embedding.add(two.getVar(1));
		AffineMorphism<Fraction> openImmersion = AffineMorphism.fromPolynomials(onePointMissing, affineLine, embedding);
		GeneralRationalFunction<Fraction, ProjectivePoint<Fraction>, ProjectivePoint<Fraction>> inverse = new GeneralRationalFunction<>(
				line, line, openImmersion, 0, openImmersion, 1);
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.one(), q.one())),
				inverse.evaluate(new ProjectivePoint<>(q, q.one(), q.one())));
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.one(), q.getInteger(5))),
				inverse.evaluate(new ProjectivePoint<>(q, q.getInteger(5), q.one())));
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.one(), q.zero())),
				inverse.evaluate(new ProjectivePoint<>(q, q.zero(), q.one())));
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.zero(), q.one())),
				inverse.evaluate(new ProjectivePoint<>(q, q.one(), q.zero())));
	}

	@Test
	void evaluateCurveTest() throws IOException {
		Rationals q = Rationals.q();
		ProjectiveLine<Fraction> line = new ProjectiveLine<>(q);
		AffineScheme<Fraction> affineLine = line.getAffineCover().get(1);
		PolynomialRing<Fraction> two = AbstractPolynomialRing.getPolynomialRing(q, 2, Monomial.GREVLEX);
		EllipticCurve<Fraction> curve = new EllipticCurve<>(q, q.getInteger(-1), q.zero());
		AffineScheme<Fraction> affineCurve = curve.getAffineCover().get(2);
		List<Polynomial<Fraction>> cover = new ArrayList<>();
		cover.add(two.getVar(1));
		AffineMorphism<Fraction> twoCover = AffineMorphism.fromPolynomials(affineCurve, affineLine, cover);
		GeneralRationalFunction<Fraction, ProjectivePoint<Fraction>, ProjectivePoint<Fraction>> rationalTwoCover = new GeneralRationalFunction<>(
				curve, line, twoCover, 2, affineCurve.identityMorphism(), 1);
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.one(), q.one())),
				rationalTwoCover.evaluate(new ProjectivePoint<>(q, q.one(), q.zero(), q.one())));
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.zero(), q.one())),
				rationalTwoCover.evaluate(new ProjectivePoint<>(q, q.zero(), q.zero(), q.one())));
		assertEquals(Optional.of(new ProjectivePoint<>(q, q.one(), q.zero())),
				rationalTwoCover.evaluate(new ProjectivePoint<>(q, q.zero(), q.one(), q.zero())));
		List<Polynomial<Fraction>> projPolynomials = new ArrayList<>();
		projPolynomials.add(curve.asGenericProjectiveScheme().homogenousPolynomialRing().getVar(1));
		projPolynomials.add(curve.asGenericProjectiveScheme().homogenousPolynomialRing().getVar(3));
		ProjectiveMorphism<Fraction> projTwoCover = new ProjectiveMorphism<>(curve.asGenericProjectiveScheme(),
				line.asGenericProjectiveScheme(), projPolynomials);
		assertEquals(new ProjectivePoint<>(q, q.one(), q.one()),
				projTwoCover.evaluate(new ProjectivePoint<>(q, q.one(), q.zero(), q.one())));
		assertEquals(new ProjectivePoint<>(q, q.zero(), q.one()),
				projTwoCover.evaluate(new ProjectivePoint<>(q, q.zero(), q.zero(), q.one())));
		assertEquals(new ProjectivePoint<>(q, q.one(), q.zero()),
				projTwoCover.evaluate(new ProjectivePoint<>(q, q.zero(), q.one(), q.zero())));
	}

}
