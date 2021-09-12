package varieties.affine;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

class AffineMorphismTest {

	@Test
	void testPreimage() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> one = q.getUnivariatePolynomialRing();
		PolynomialRing<Fraction> two = AbstractPolynomialRing.getPolynomialRing(q, 2, Monomial.GREVLEX);
		AffineScheme<Fraction> line = new AffineScheme<>(q, new CoordinateRing<>(one, one.getZeroIdeal()));
		AffineScheme<Fraction> plane = new AffineScheme<>(q, new CoordinateRing<>(two, two.getZeroIdeal()));
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
		System.out.println(morphism.preimage(new AffinePoint<>(q, q.zero(), q.zero())).image().getDomain());
		System.out.println(morphism.preimage(new AffinePoint<>(q, q.getInteger(3), q.getInteger(6))).getDomain());
		List<CoordinateRingElement<Fraction>> pointList = new ArrayList<>();
		pointList.add(plane.getCoordinateRing().getVar(1));
		pointList.add(plane.getCoordinateRing().getVar(2));
		AffineMorphism<Fraction> embedding = AffineScheme.restrictAwayFrom(plane, pointList);
		AffineMorphism<Fraction> pre = morphism.preimage(embedding);
		System.out.println(pre.getDomain());
		System.out.println(pre.image().getDomain());
		System.out.println(morphism.image().getDomain());
		List<Polynomial<Fraction>> lineList = new ArrayList<>();
		lineList.add(one.getVar());
		lineList.add(one.multiply(2, one.getVar()));//subtract(one.getVar(), one.one()));
		AffineMorphism<Fraction> lineMorphism = AffineMorphism.fromPolynomials(line, plane, lineList);
		System.out.println(lineMorphism.image().getDomain());
		pre = morphism.preimage(lineMorphism);
		System.out.println(pre.getDomain());
		System.out.println(pre.image().getDomain());
	}

}
