package varieties;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.affine.AffineScheme;
import varieties.curves.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class SchemesTest {

	@Test
	void linesTest() {
		FiniteField fq = FiniteField.getFiniteField(25);
		PolynomialRing<FFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.GREVLEX);
		Polynomial<FFE> xy = polynomialRing.getEmbedding(fq.one(), new int[] { 1, 1 });
		AffineScheme<FFE> twoLines = new AffineScheme<>(fq,
				polynomialRing.getIdeal(Collections.singletonList(xy)).divideOut());
		for (AffineScheme<FFE> line : twoLines.irreducibleComponents()) {
			System.out.println(line);
		}
	}

	@Test
	void pointTest() {
		FiniteField fq = FiniteField.getFiniteField(25);
		PolynomialRing<FFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.GREVLEX);
		Polynomial<FFE> xy = polynomialRing.getEmbedding(fq.one(), new int[] { 1, 1 });
		Polynomial<FFE> xpy = polynomialRing.add(polynomialRing.getVar(1), polynomialRing.getVar(2));
		List<Polynomial<FFE>> list = new ArrayList<>();
		list.add(xy);
		list.add(xpy);
		AffineScheme<FFE> doublePoint = new AffineScheme<>(fq, polynomialRing.getIdeal(list).divideOut());
		for (AffineScheme<FFE> point : doublePoint.irreducibleComponents()) {
			System.out.println(point);
		}
	}

	@Test
	void pointsTest() {
		FiniteField fq = FiniteField.getFiniteField(25);
		PolynomialRing<FFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.GREVLEX);
		List<Polynomial<FFE>> list1 = new ArrayList<>();
		list1.add(polynomialRing.getVar(1));
		list1.add(polynomialRing.getVar(2));
		AffineScheme<FFE> point1 = new AffineScheme<>(fq, polynomialRing.getIdeal(list1).divideOut());
		List<Polynomial<FFE>> list2 = new ArrayList<>();
		list2.add(polynomialRing.subtract(polynomialRing.getVar(1), polynomialRing.one()));
		list2.add(polynomialRing.subtract(polynomialRing.getVar(2), polynomialRing.one()));
		AffineScheme<FFE> point2 = new AffineScheme<>(fq, polynomialRing.getIdeal(list2).divideOut());
		List<Polynomial<FFE>> list3 = new ArrayList<>();
		list3.add(polynomialRing.getVar(1));
		list3.add(polynomialRing.subtract(polynomialRing.getVar(2), polynomialRing.one()));
		AffineScheme<FFE> point3 = new AffineScheme<>(fq, polynomialRing.getIdeal(list3).divideOut());
		List<Polynomial<FFE>> list4 = new ArrayList<>();
		list4.add(polynomialRing.subtract(polynomialRing.getVar(1), polynomialRing.one()));
		list4.add(polynomialRing.getVar(2));
		AffineScheme<FFE> point4 = new AffineScheme<>(fq, polynomialRing.getIdeal(list4).divideOut());
		AffineScheme<FFE> union1 = AffineScheme.union(point1, point2).getUnion();
		AffineScheme<FFE> union2 = AffineScheme.union(point3, point4).getUnion();
		AffineScheme<FFE> union = AffineScheme.union(union1, union2).getUnion();
		for (AffineScheme<FFE> point : union.irreducibleComponents()) {
			System.out.println(point);
		}
	}

	@Test
	void intersectionTest() {
		FiniteField fq = FiniteField.getFiniteField(25);
		PolynomialRing<FFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fq, 2, Monomial.GREVLEX);
		Polynomial<FFE> polynomial = polynomialRing.getVarPower(2, 2);
		polynomial = polynomialRing.subtract(polynomial, polynomialRing.getVarPower(1, 3));
		polynomial = polynomialRing.subtract(polynomial, polynomialRing.getVar(1));
		polynomial = polynomialRing.subtract(polynomial, polynomialRing.one());
		EllipticCurve<FFE> ec = new EllipticCurve<>(fq, fq.one(), fq.one());
		for (ProjectivePoint<FFE> point : ec) {
			System.out.println(point);
		}
		Polynomial<FFE> line = polynomialRing.add(
				polynomialRing.multiply(fq.subtract(fq.multiply(2, fq.alpha()), fq.one()), polynomialRing.getVar(1)),
				polynomialRing.multiply(fq.negative(fq.add(fq.one(), fq.alpha())), polynomialRing.getVar(2)),
				polynomialRing.getEmbedding(fq.add(fq.one(), fq.alpha())));
		AffineScheme<FFE> curve = new AffineScheme<>(fq,
				polynomialRing.getIdeal(Collections.singletonList(polynomial)).divideOut());
		AffineScheme<FFE> l = new AffineScheme<>(fq,
				polynomialRing.getIdeal(Collections.singletonList(line)).divideOut());
		AffineScheme<FFE> intersection = AffineScheme.intersect(curve, l).getIntersection();
		for (AffineScheme<FFE> point : intersection.irreducibleComponents()) {
			System.out.println(point);
		}
	}
}
