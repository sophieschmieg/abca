package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class KernelPointIsogeny<T extends Element<T>> implements Isogeny<T> {
	private Field<T> field;
	private EllipticCurve<T> domain;
	private EllipticCurve<T> range;

	private ProjectivePoint<T> kernelGenerator;
	private int kernelDegree;

	private ProjectiveMorphism<T> asMorphism;

	public KernelPointIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernelGenerator, int kernelDegree) {
		this.field = domain.getField();
		this.domain = domain;
		if (!domain.multiply(kernelDegree, kernelGenerator).equals(domain.neutral())) {
			throw new ArithmeticException("kernel degree wrong!");
		}
		this.kernelGenerator = kernelGenerator;
		this.kernelDegree = kernelDegree;

		T a = domain.getA();
		T b = domain.getB();
		T v = field.zero();
		T w = field.zero();
		for (int i = 1; i < (kernelDegree + 1) / 2; i++) {
			ProjectivePoint<T> kernelPoint = domain.multiply(i, kernelGenerator);
			if (kernelPoint.equals(domain.neutral())) {
				throw new ArithmeticException("kernel degree wrong!");
			}
			T pointX = kernelPoint.getDehomogenisedCoord(1, 3);
			T pointY = kernelPoint.getDehomogenisedCoord(2, 3);
			T gxp = field.add(field.multiply(3, pointX, pointX), a);
			T vp = field.multiply(2, gxp);
			T gyp = field.multiply(-2, pointY);
			T up = field.multiply(gyp, gyp);
			v = field.add(v, vp);
			w = field.add(w, up, field.multiply(pointX, vp));
		}
		if (kernelDegree % 2 == 0) {
			ProjectivePoint<T> kernelPoint = domain.multiply(kernelDegree / 2, kernelGenerator);
			if (kernelPoint.equals(domain.neutral())) {
				throw new ArithmeticException("kernel degree wrong!");
			}
			T pointX = kernelPoint.getDehomogenisedCoord(1, 3);
			T gxp = field.add(field.multiply(3, pointX, pointX), a);
			v = field.add(v, gxp);
			w = field.add(w, field.multiply(pointX, gxp));

		}
		T rangeA = field.subtract(a, field.multiply(5, v));
		T rangeB = field.subtract(b, field.multiply(7, w));
		this.range = new EllipticCurve<T>(field, rangeA, rangeB);
	}

	public EllipticCurve<T> getDomain() {
		return domain;
	}

	public EllipticCurve<T> getRange() {
		return range;
	}

	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		if (!domain.hasRationalPoint(point)) {
			throw new ArithmeticException("Point not on curve");
		}
		if (point.equals(domain.neutral())) {
			return range.neutral();
		}
		T x = point.getDehomogenisedCoord(1, 3);
		T y = point.getDehomogenisedCoord(2, 3);
		for (int i = 1; i < kernelDegree; i++) {
			ProjectivePoint<T> kernelPoint = domain.multiply(i, kernelGenerator);
			if (kernelPoint.equals(point)) {
				return range.neutral();
			}
			ProjectivePoint<T> added = domain.add(point, kernelPoint);
			if (added.equals(domain.neutral())) {
				return range.neutral();
			}
			x = field.add(x, added.getDehomogenisedCoord(1, 3));
			x = field.subtract(x, kernelPoint.getDehomogenisedCoord(1, 3));
			y = field.add(y, added.getDehomogenisedCoord(2, 3));
			y = field.subtract(y, kernelPoint.getDehomogenisedCoord(2, 3));
		}
		ProjectivePoint<T> result = new ProjectivePoint<T>(field, x, y, field.one());
		return result;
	}

	public BigInteger getDegree() {
		return BigInteger.valueOf(kernelDegree);
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		if (asMorphism == null) {
			PolynomialRing<T> polynomials = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
			Polynomial<T> x = polynomials.getVar(1);
			Polynomial<T> y = polynomials.getVar(2);
			Polynomial<T> z = polynomials.getVar(3);
			for (int i = 1; i < kernelDegree; i++) {
				ProjectivePoint<T> kernelPoint = domain.multiply(i, kernelGenerator);
				ProjectiveMorphism<T> added = domain.translationMorphism(kernelPoint);
				Polynomial<T> newDenom = added.asPolynomials().get(2);
				x = polynomials.add(polynomials.multiply(x, newDenom),
						polynomials.multiply(added.asPolynomials().get(0), z));
				y = polynomials.add(polynomials.multiply(y, newDenom),
						polynomials.multiply(added.asPolynomials().get(1), z));
				z = polynomials.multiply(z, newDenom);
				x = polynomials.subtract(x, polynomials.multiply(kernelPoint.getDehomogenisedCoord(1, 3), z));
				y = polynomials.subtract(y, polynomials.multiply(kernelPoint.getDehomogenisedCoord(2, 3), z));
			}
			List<Polynomial<T>> result = new ArrayList<>();
			result.add(x);
			result.add(y);
			result.add(z);
			asMorphism = new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), range.asGenericProjectiveScheme(), result);
		}
		return asMorphism;
	}

	public Isogeny<T> getDual() {
		ProjectivePoint<T> dualKernelPoint = null;
		torsionPointsLoop: for (ProjectivePoint<T> torsionPoint : domain.getTorsionPoints(kernelDegree)) {
			for (int i = 1; i < kernelDegree; i++) {
				if (domain.multiply(i, torsionPoint).equals(domain.neutral())) {
					continue torsionPointsLoop;
				}
			}
			if (!evaluate(torsionPoint).equals(range.neutral())) {
				dualKernelPoint = torsionPoint;
				break;
			}
		}
		if (dualKernelPoint == null) {
			throw new ArithmeticException("Dual isogeny not definied over this field. Please extend the field.");
		}
		return new KernelPointIsogeny<T>(range, dualKernelPoint, kernelDegree);
	}
}
