package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Field.Extension;
import varieties.curves.elliptic.EllipticCurve.EllipticCurveExtensionResult;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class KernelPointIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private Field<T> field;
	private EllipticCurve<T> domain;
	private EllipticCurve<T> range;

	private ProjectivePoint<T> kernelGenerator;
	private int kernelDegree;

	private ProjectiveMorphism<T> asMorphism;

	public KernelPointIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernelGenerator, int kernelDegree) {
		this.field = domain.getField();
		this.domain = domain;
		if (!Integers.z().isPrime(Integers.z().getInteger(kernelDegree))) {
			throw new ArithmeticException("Not of prime degree!");
		}
		if (kernelGenerator.equals(domain.neutral())) {
			throw new ArithmeticException("Kernel point is not generator!");
		}
		if (!domain.multiply(kernelDegree, kernelGenerator).equals(domain.neutral())) {
			throw new ArithmeticException("kernel degree wrong! [" + kernelDegree + "]" + kernelGenerator + " = "
					+ domain.multiply(kernelDegree, kernelGenerator));
		}
		this.kernelGenerator = kernelGenerator;
		this.kernelDegree = kernelDegree;

		T a1 = domain.getA1();
		T a2 = domain.getA2();
		T a3 = domain.getA3();
		T a4 = domain.getA4();
		T a6 = domain.getA6();
		T t = field.zero();
		T u = field.zero();
		for (int i = 1; i < (kernelDegree + 1) / 2; i++) {
			ProjectivePoint<T> kernelPoint = domain.multiply(i, kernelGenerator);
			if (kernelPoint.equals(domain.neutral())) {
				throw new ArithmeticException("kernel degree wrong!");
			}
			T pointX = kernelPoint.getDehomogenisedCoord(1, 3);
			T pointY = kernelPoint.getDehomogenisedCoord(2, 3);
			T gxp = field.add(field.add(field.multiply(3, pointX, pointX), field.multiply(2, a2, pointX), a4),
					field.multiply(-1, a1, pointY));
			T gyp = field.add(field.multiply(-2, pointY), field.multiply(-1, a1, pointX), field.multiply(-1, a3));
			T tp = field.add(field.multiply(2, gxp), field.multiply(-1, a1, gyp));
			T up = field.multiply(gyp, gyp);
			t = field.add(t, tp);
			u = field.add(u, up, field.multiply(pointX, tp));
		}
		if (kernelDegree % 2 == 0) {
			ProjectivePoint<T> kernelPoint = domain.multiply(kernelDegree / 2, kernelGenerator);
			if (kernelPoint.equals(domain.neutral())) {
				throw new ArithmeticException("kernel degree wrong!");
			}
			T pointX = kernelPoint.getDehomogenisedCoord(1, 3);
			T pointY = kernelPoint.getDehomogenisedCoord(2, 3);
			T gxp = field.add(field.add(field.multiply(3, pointX, pointX), field.multiply(2, a2, pointX), a4),
					field.multiply(-1, a1, pointY));
			t = field.add(t, gxp);
			u = field.add(u, field.multiply(pointX, gxp));
		}
		T rangeA4 = field.subtract(a4, field.multiply(5, t));
		T rangeA6 = field.subtract(a6, field.add(field.multiply(7, u), field.multiply(domain.getB2(), t)));
		this.range = new EllipticCurve<T>(field, a1, a2, a3, rangeA4, rangeA6);
		this.range.setNumberOfPointsFrom(domain, this);
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof KernelPointIsogeny<?>) {
			KernelPointIsogeny<T> other = (KernelPointIsogeny<T>) o;
			int cmp = domain.compareTo(other.domain);
			if (cmp != 0) {
				return cmp;
			}
			return kernelGenerator.compareTo(other.kernelGenerator);
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public String toString() {
		return "Ker = <" + kernelGenerator.toString() + ">";
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
	
	@Override
	public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> ProjectivePoint<E> evaluateOverExtension(
			ProjectivePoint<E> point, Extension<T, B, E, FE> extension) {
		EllipticCurveExtensionResult<T, B, E, FE> extendedDomain = domain.extendBaseField(extension);
		EllipticCurveExtensionResult<T, B, E, FE> extendedRange = range.extendBaseField(extension);
		if (!extendedDomain.getCurve().hasRationalPoint(point)) {
			throw new ArithmeticException("Point not on curve");
		}
		if (point.equals(extendedDomain.getCurve().neutral())) {
			return extendedRange.getCurve().neutral();
		}
		FE field = extension.extension();
		E x = point.getDehomogenisedCoord(1, 3);
		E y = point.getDehomogenisedCoord(2, 3);
		for (int i = 1; i < kernelDegree; i++) {
			ProjectivePoint<E> kernelPoint = extendedDomain.getCurve().multiply(i,extendedDomain.getEmbedding().evaluate(kernelGenerator));
			if (kernelPoint.equals(point)) {
				return extendedRange.getCurve().neutral();
			}
			ProjectivePoint<E> added = extendedDomain.getCurve().add(point, kernelPoint);
			if (added.equals(domain.neutral())) {
				return extendedRange.getCurve().neutral();
			}
			x = field.add(x, added.getDehomogenisedCoord(1, 3));
			x = field.subtract(x, kernelPoint.getDehomogenisedCoord(1, 3));
			y = field.add(y, added.getDehomogenisedCoord(2, 3));
			y = field.subtract(y, kernelPoint.getDehomogenisedCoord(2, 3));
		}
		ProjectivePoint<E> result = new ProjectivePoint<>(field, x, y, field.one());
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
			asMorphism = new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), range.asGenericProjectiveScheme(),
					result);
		}
		return asMorphism;
	}

	public Isogeny<T> getDual() {
		ProjectivePoint<T> dualKernelPoint = null;
		for (ProjectivePoint<T> torsionPoint : domain.getTorsionPointBasis(kernelDegree)) {
			ProjectivePoint<T> evaluated = evaluate(torsionPoint);
			if (range.getOrder(evaluated).intValueExact() == kernelDegree) {
				dualKernelPoint = evaluated;
				break;
			}
		}
		if (dualKernelPoint == null) {
			throw new ArithmeticException("Dual isogeny not definied over this field. Please extend the field.");
		}
		Isogeny<T> dualUpToIsomorphism = new KernelPointIsogeny<T>(range, dualKernelPoint, kernelDegree);
		if (dualUpToIsomorphism.getRange().equals(domain)) {
			return dualUpToIsomorphism;
		}
		ProjectivePoint<T> testPoint;
		BigInteger order;
		int tries = 0;
		BigInteger degree = BigInteger.valueOf(kernelDegree);
		do {
			testPoint = domain.getRandomElement();
			order = domain.getOrder(testPoint);
			tries++;
		} while (tries < 10 && degree.mod(order).equals(degree.negate().mod(order)));
		ProjectivePoint<T> expected = domain.multiply(kernelDegree, testPoint);
		for (Isomorphism<T> isomorphism : dualUpToIsomorphism.getRange().getIsomorphisms(domain)) {
				ProjectivePoint<T> evaluated = evaluate(testPoint);
				evaluated = dualUpToIsomorphism.evaluate(evaluated);
				evaluated = isomorphism.evaluate(evaluated);
				if (evaluated.equals(expected)) {
					return new CompositionIsogeny<>(dualUpToIsomorphism, isomorphism);
				}
		}
		throw new ArithmeticException("Dual not computed correctly");
	}

	@Override
	public List<ProjectivePoint<T>> kernelGenerators() {
		return Collections.singletonList(kernelGenerator);
	}
}
