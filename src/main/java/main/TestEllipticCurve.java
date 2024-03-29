package main;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomialRing;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class TestEllipticCurve<T extends Element<T>> {
	public TestEllipticCurve(Field<T> field, T a, T b) {
		EllipticCurve<T> curve = new EllipticCurve<T>(field, a, a);
		UnivariatePolynomialRing<T> univar = field.getUnivariatePolynomialRing();
		System.out.println("Curve: " + curve);
		for (int i = 1; i < 10; i++) {
			Polynomial<T> div = curve.getDivisionPolynomial(i);
			System.out.println(i + " Division Polynomial: " + div);
			if (i % 2 == 1) {
				Map<T, Integer> roots = field.roots(univar.getEmbedding(div, new int[] {0}));
				FactorizationResult<Polynomial<T>, T> factors = field.factorization(univar.getEmbedding(div, new int[] {0}));
				System.out.println("Roots: " + roots);
				System.out.println("Factors: " + factors);
				Polynomial<T> product = univar.one();
				for (Polynomial<T> factor : factors.primeFactors()) {
					product = univar.multiply(product, factor);
				}
				if (!product.equals(univar.getEmbedding(div, new int[] {0}))) {
					throw new RuntimeException("Factorization wrong");
				}
			}
			for (ProjectivePoint<T> point : curve.getTorsionPoints(i)) {
				System.out.println("[" + i + "]" + point.toString() + " = " + curve.neutral());
			}
		}
		int num = curve.getNumberOfElements().intValue();
		System.out.println("Number of elements: " + num);
		System.out.println("Field: " + field);
		System.out.println("Field Characteristic: " + field.characteristic());
		System.out.println("Number of Field Elements: " + field.getNumberOfElements());
		List<Integer> divisor = new ArrayList<>();
		for (int i = 1; i <= num; i++) {
			if (num % i == 0) {
				divisor.add(i);
			}
		}
		int counter = 1;
		for (ProjectivePoint<T> p : curve) {
			int i = -1;
			for (int d : divisor) {
				if (curve.multiply(d, p).equals(curve.neutral())) {
					i = d;
					break;
				}
			}
			if (i > 2 && i < 40) {
				Polynomial<T> divPoly = curve.getDivisionPolynomial(i);
				System.out.println(i + " division polynomial evaluated at " + p + ": " + divPoly.getPolynomialRing().evaluate(divPoly, p.getDehomogenous(3).getCoords()));
			}
			//			T t = tau(p);
			System.out.println(counter + " Point: " + p + " Order: " + i); //+ " Tau: " + t/ + " Root: " + this.field.hasRoot(t, i));
			counter++;
			if (counter > 200) {
				break;
			}
		}
	}
}
