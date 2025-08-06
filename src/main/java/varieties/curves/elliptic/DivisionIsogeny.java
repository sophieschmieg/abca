package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.FieldExtension.SplittingFieldResult;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import varieties.curves.elliptic.EllipticCurve.EllipticCurveExtensionResult;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class DivisionIsogeny<T extends Element<T>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>>
		implements Isogeny<T> {
	private Fraction n;
	private EllipticCurve<T> curve;
	private Extension<T, B, E, FE> extension;
	private EllipticCurve<E> extendedCurve;
	private Isogeny<T> base;
	private Isogeny<E> equivalent;

	public static <T extends Element<T>> Isogeny<T> getDivisionIsogeny(Isogeny<T> base, BigInteger degree, BigInteger divisor) {
		Rationals q = Rationals.q();
		Field<T> field = base.getDomain().getField();
		return new DivisionIsogeny<>(q.inverse(q.getInteger(divisor)), base, degree, field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	
	public static <T extends Element<T>> Isogeny<T> getDivisionIsogeny(Fraction n, Isogeny<T> base) {
		Field<T> field = base.getDomain().getField();
		return new DivisionIsogeny<>(n, base, field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	public static <T extends Element<T>> Isogeny<T> getDivisionIsogeny(Fraction n, Isogeny<T> base,
			BigInteger degree) {
		Field<T> field = base.getDomain().getField();
		return new DivisionIsogeny<>(n, base, degree, field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	public DivisionIsogeny(Fraction n, Isogeny<T> base, Extension<T, B, E, FE> trivialExtension) {
		this(n, base, base.getDegree(), trivialExtension);
	}

	public DivisionIsogeny(Fraction n, Isogeny<T> base, BigInteger degree, Extension<T, B, E, FE> trivialExtension) {
		this.n = n;
		this.base = base;
		this.curve = base.getDomain();
		Polynomial<T> baseDivisionPolynomial = curve.getDivisionPolynomial(
				degree.multiply(n.getNumerator().getValue()).multiply(n.getDenominator().getValue()).intValueExact());
		if (baseDivisionPolynomial.degree(2) > 0) {
			baseDivisionPolynomial = baseDivisionPolynomial.getPolynomialRing().divide(baseDivisionPolynomial,
					baseDivisionPolynomial.getPolynomialRing().getVar(2));
		}
		baseDivisionPolynomial = curve.getField().getUnivariatePolynomialRing().getEmbedding(baseDivisionPolynomial,
				new int[1]);
		FE ext = trivialExtension.extension();
		UnivariatePolynomial<E> divisionPolynomial = ext.getUnivariatePolynomialRing()
				.getEmbedding(baseDivisionPolynomial, trivialExtension.embeddingMap());
		SplittingFieldResult<B, E, FE> splittingField = trivialExtension.extension().getSplittingField(divisionPolynomial);
		this.extension = trivialExtension.extendFurther(splittingField.getExtension());
		this.extendedCurve = curve.extendBaseField(extension).getCurve();
		BigInteger divisor = n.getDenominator().getValue();
		List<ProjectivePoint<E>> divisorBasis = extendedCurve.getTorsionPointBasis(divisor);
		if (divisorBasis.size() != 2) {
			throw new ArithmeticException("Not enough divisors to divide!");
		}
		for (ProjectivePoint<E> point : divisorBasis) {
			if (!base.evaluateOverExtension(point, extension).equals(extendedCurve.neutral())) {
				throw new ArithmeticException("isogeny not divisible!");
			}
		}
		List<ProjectivePoint<E>> degreeBasis = extendedCurve
				.getTorsionPointBasis(degree.multiply(n.getNumerator().getValue()));
		if (degreeBasis.isEmpty()) {
			throw new ArithmeticException("No degree torsion points!");
		}
		for (ProjectivePoint<E> point : degreeBasis) {
			if (!extendedCurve
					.multiply(n.getNumerator(), base.getDual()
							.evaluateOverExtension(base.evaluateOverExtension(point, extension), extension))
					.equals(extendedCurve.neutral())) {
				throw new ArithmeticException("Degree wrong!");
			}
		}
		if (!degree.multiply(n.getNumerator().getValue()).mod(divisor).equals(BigInteger.ZERO)) {
			throw new ArithmeticException("degree not divisible");
		}
		List<ProjectivePoint<E>> torsionPointBasis = extendedCurve.getTorsionPointBasis(degree.divide(divisor));
		if (torsionPointBasis.size() != 2) {
			throw new ArithmeticException("Not enough torsion points to divide!");
		}
		Set<ProjectivePoint<E>> kernelBasis = new TreeSet<>();
		for (ProjectivePoint<E> point : torsionPointBasis) {
			kernelBasis.add(extendedCurve.multiply(n.getNumerator(), base.getDual().evaluateOverExtension(point, extension)));
		}
		Isogeny<E> kernelIsogeny = new KernelIsogeny<>(extendedCurve, kernelBasis);
		List<Isomorphism<E>> isomorphisms = kernelIsogeny.getRange().getIsomorphisms(extendedCurve);
		ProjectivePoint<E> testPoint;
		int tries = 0;
		BigInteger order;
		do {
			testPoint = extendedCurve.getRandomElement();
			order = extendedCurve.getOrder(testPoint);
			tries++;
		} while (tries < 10 && order.mod(divisor).equals(BigInteger.ZERO)
				&& order.divide(divisor).equals(BigInteger.TWO));
		if (tries == 10) {
			throw new ArithmeticException("Could not find test point!");
		}
		for (Isomorphism<E> candidate : isomorphisms) {
			if (extendedCurve.multiply(divisor, candidate.evaluate(kernelIsogeny.evaluate(testPoint)))
					.equals(base.evaluateOverExtension(testPoint, extension))) {
				equivalent = new CompositionIsogeny<E>(kernelIsogeny, candidate);
			}
		}
	}
	
	private DivisionIsogeny(Fraction n, Isogeny<T> base, Isogeny<E> equivalent, Extension<T, B, E, FE> extension) {
	this.n = n;
		this.base = base;
		this.equivalent = equivalent;
		this.extension = extension;
		this.extendedCurve = base.getDomain().extendBaseField(extension).getCurve();
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof DivisionIsogeny<?, ?, ?, ?>) {
			DivisionIsogeny<T, B, E, FE> other = (DivisionIsogeny<T, B, E, FE>) o;
			int cmp = base.compareTo(other.base);
			if (cmp != 0) {
				return cmp;
			}
			return n.compareTo(other.n);
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public String toString() {
		return "[" + n + "]" + base;
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return curve;
	}

	@Override
	public EllipticCurve<T> getRange() {
		return base.getRange();
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		EllipticCurveExtensionResult<T, B, E, FE> ext = curve.extendBaseField(extension);
		ProjectivePoint<E> evaluated = equivalent.evaluate(ext.getEmbedding().evaluate(point));
		List<T> coords = new ArrayList<>();
		for (E resultCoord : evaluated.getCoords()) {
			coords.add(extension.retractionMap().evaluate(resultCoord));
		}
		return new ProjectivePoint<>(curve.getField(), coords);
	}
	
	@Override
	public <BC extends Element<BC>, EC extends AlgebraicExtensionElement<BC, EC>, FEC extends FieldExtension<BC, EC, FEC>> ProjectivePoint<EC> evaluateOverExtension(
			ProjectivePoint<EC> point, Extension<T, BC, EC, FEC> extension) {
	throw new ArithmeticException("Not implemented");
		//	return equivalent.evaluateOverExtension(point, extension);
	}

	@Override
	public BigInteger getDegree() {
		Rationals q = Rationals.q();
		return q.multiply(base.getDegree(), n, n).asInteger().getValue();
	}

	@Override
	public Isogeny<T> getDual() {
		return new DivisionIsogeny<>(n, base.getDual(), equivalent.getDual(), extension);
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		throw new ArithmeticException("not implemented");
//		if (asMorphism == null) {
//			asMorphism = curve.multiplicationMorphism(n.intValueExact());
//		}
//		return asMorphism;
	}

	@Override
	public List<ProjectivePoint<T>> kernelGenerators() {
		throw new ArithmeticException("not implemented");
	}

}
