package fields.floatingpoint;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.FieldAutomorphism;
import fields.helper.GaloisGroup;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.FloatingPointSet;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.ValueField;
import fields.vectors.Vector;
import util.Identity;

public class Complex extends AbstractFieldExtension<Real, ComplexNumber, Complex>
		implements ValueField<ComplexNumber>, FloatingPointSet<ComplexNumber, Complex> {
	private Reals r;
	private Complex doublePrecision;
	private static Map<Integer, Complex> complex = new TreeMap<>();

	public class ComplexNumber extends AbstractElement<ComplexNumber>
			implements AlgebraicExtensionElement<Real, ComplexNumber> {
		private Real realPart;
		private Real complexPart;

		public ComplexNumber(Real realPart, Real complexPart) {
			this.realPart = realPart;
			this.complexPart = complexPart;
		}

		@Override
		public int compareTo(ComplexNumber o) {
			int realComp = realPart.compareTo(o.realPart);
			if (realComp != 0) {
				return realComp;
			}
			return complexPart.compareTo(o.complexPart);
		}

		public Real realPart() {
			return realPart;
		}

		public Real complexPart() {
			return complexPart;
		}

		public String toString() {
			return realPart.toString() + " + " + complexPart.toString() + "i";
		}

		@Override
		public UnivariatePolynomial<Real> asPolynomial() {
			return r.getUnivariatePolynomialRing().getPolynomial(realPart, complexPart);
		}

	}

	private Complex(Reals reals) {
		super(reals.getUnivariatePolynomialRing().getPolynomial(reals.one(), reals.zero(), reals.one()), reals);
		this.r = reals;
	}

	private Complex(int precision) {
		this(Reals.r(precision));
	}
	
	public static Complex c(int precision) {
		if (!complex.containsKey(precision)) {
			complex.put(precision, new Complex(precision));
		}
		return complex.get(precision);
	}

	private Complex doublePrecision() {
		if (doublePrecision == null) {
			doublePrecision = withPrecision(2 * precision());
		}
		return doublePrecision;
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public boolean close(ComplexNumber t1, ComplexNumber t2) {
		int minRealExponent = Math.min(t1.realPart.exponent(), t2.realPart.exponent());
		int minImgExponent = Math.min(t1.complexPart.exponent(), t2.complexPart.exponent());
		int minExponent = Math.max(minRealExponent, minImgExponent);
		int precision = r.precision();
		return r.abs(r.subtract(t1.realPart, t2.realPart)).compareTo(r.getPowerOfTwo(-precision + minExponent + 1)) <= 0
				&& r.abs(r.subtract(t1.complexPart, t2.complexPart))
						.compareTo(r.getPowerOfTwo(-precision + minExponent + 1)) <= 0;
	}

	@Override
	public Complex withPrecision(int precision) {
		return new Complex(precision);
	}

	@Override
	public int precision() {
		return r.precision();
	}

	@Override
	public ComplexNumber roundToFixedPoint(ComplexNumber t, int precision) {
		return getNumber(r.roundToFixedPoint(t.realPart, precision), r.roundToFixedPoint(t.complexPart, precision));
	}

	public ComplexNumber getNumber(Real realPart, Real complexPart) {
		return new ComplexNumber(realPart, complexPart);
	}

	public ComplexNumber getDouble(double realPart, double complexPart) {
		return new ComplexNumber(r.getDouble(realPart), r.getDouble(complexPart));
	}

	public ComplexNumber getEmbedding(ComplexNumber t) {
		return getNumber(r.getEmbedding(t.realPart), r.getEmbedding(t.complexPart));
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	public String toString() {
		return "C";
	}

	@Override
	public Extension<ComplexNumber, Real, ComplexNumber, Complex> getExtension(
			UnivariatePolynomial<ComplexNumber> minimalPolynomial) {
		return new Extension<>(this, this, new Identity<>(), new MathMap<ComplexNumber, Vector<ComplexNumber>>() {
			@Override
			public Vector<ComplexNumber> evaluate(ComplexNumber t) {
				return new Vector<>(Collections.singletonList(t));
			}
		});

	}

	public ComplexNumber i() {
		return alpha();
	}

	public ComplexNumber e() {
		return getEmbedding(r.e());
	}

	public ComplexNumber pi() {
		return getEmbedding(r.pi());
	}

	public Real arg(ComplexNumber t) {
		return r.arctan2(t.complexPart, t.realPart);
	}

	@Override
	public ComplexNumber getInteger(BigInteger t) {
		return getEmbedding(r.getInteger(t));
	}

	public ComplexNumber getEmbedding(IntE t) {
		return getEmbedding(r.getEmbedding(t));
	}

	public ComplexNumber getEmbedding(Fraction t) {
		return getEmbedding(r.getEmbedding(t));
	}

	public ComplexNumber getEmbedding(Real t) {
		return new ComplexNumber(t, r.zero());
	}

	@Override
	public ComplexNumber add(ComplexNumber t1, ComplexNumber t2) {
		return new ComplexNumber(r.add(t1.realPart, t2.realPart), r.add(t1.complexPart, t2.complexPart));
	}

	@Override
	public ComplexNumber negative(ComplexNumber t) {
		return new ComplexNumber(r.negative(t.realPart), r.negative(t.complexPart));
	}

	@Override
	public ComplexNumber multiply(ComplexNumber t1, ComplexNumber t2) {
		return new ComplexNumber(
				r.subtract(r.multiply(t1.realPart, t2.realPart), r.multiply(t1.complexPart, t2.complexPart)),
				r.add(r.multiply(t1.realPart, t2.complexPart), r.multiply(t1.complexPart, t2.realPart)));
	}

	@Override
	public ComplexNumber inverse(ComplexNumber t) {
		ComplexNumber conjugate = conjugate(t);
		Real denom = norm(t);
		return new ComplexNumber(r.divide(conjugate.realPart, denom), r.divide(conjugate.complexPart, denom));
	}

	@Override
	public ComplexNumber divide(ComplexNumber dividend, ComplexNumber divisor) {
		Real denom = doublePrecision().norm(divisor);
		return getEmbedding(doublePrecision().multiply(doublePrecision().getEmbedding(doublePrecision().r.inverse(denom)), dividend, doublePrecision().conjugate(divisor)));
	}

	public ComplexNumber conjugate(ComplexNumber t) {
		return new ComplexNumber(t.realPart, r.negative(t.complexPart));
	}

	@Override
	public ComplexNumber getRandomElement() {
		return new ComplexNumber(r.getRandomElement(), r.getRandomElement());
	}

	@Override
	public Real value(ComplexNumber t) {
		return r.positiveSqrt(norm(t));
	}

	public ComplexNumber exp(ComplexNumber t) {
		return getEmbedding(doublePrecision().calculateExp(t, this));
	}

	private ComplexNumber calculateExp(ComplexNumber t, Complex halfPrecision) {
		ComplexNumber multiplier = one();
		if (r.abs(t.realPart()).compareTo(r.one()) > 0) {
			IntE offset = t.realPart().round();
			multiplier = power(halfPrecision.e(), offset.getValue());
			t = subtract(t, getInteger(offset));
		}
		if (r.abs(t.complexPart()).compareTo(r.getInteger(3)) > 0) {
			Real multiplesOfTwoPi = r.divide(t.complexPart(), r.multiply(2, halfPrecision.r.pi()));
			IntE rounded = multiplesOfTwoPi.round();
			t = subtract(t, multiply(Integers.z().multiply(2, rounded), halfPrecision.pi(), i()));
		}
		ComplexNumber result = one();
		ComplexNumber power = one();
		ComplexNumber prev;
		int i = 1;
		do {
			prev = result;
			power = multiply(power, divide(t, getInteger(i)));
			result = add(result, power);
			i++;
		} while (!close(result, prev));
		return multiply(multiplier, result);
	}

	public ComplexNumber log(ComplexNumber t) {
		return new ComplexNumber(r.log(r.positiveSqrt(norm(t))), arg(t));
	}

	@Override
	public Real norm(ComplexNumber t) {
		return r.add(r.multiply(t.realPart(), t.realPart()), r.multiply(t.complexPart(), t.complexPart()));
	}

	public ComplexNumber findRoot(Polynomial<ComplexNumber> t) {
		UnivariatePolynomialRing<ComplexNumber> ring = getUnivariatePolynomialRing();
//		UnivariatePolynomial<ComplexNumber> p = ring.toUnivariate(t);
		UnivariatePolynomial<ComplexNumber> derivative = ring.derivative(t);
//		Real scale = r.zero();
//		Real lcScale = value(t.leadingCoefficient());
//		for (int i = 0; i < t.degree(); i++) {
//			scale += Math.pow(value(p.univariateCoefficient(i)) / lcScale, 1.0 / (t.degree() - i));
//		}
		int maxNumIterations = 1000;
		while (true) {
			ComplexNumber x = getRandomElement();// multiply(getEmbedding(r.getDouble(scale)), getRandomElement());
			ComplexNumber prevX = zero();
			for (int i = 0; i < maxNumIterations; i++) {
				prevX = x;
				x = subtract(x, divide(ring.evaluate(t, x), ring.evaluate(derivative, x)));
				if (ring.evaluate(t, x).equals(zero()) || close(x, prevX)) {
					return x;
				}
			}
		}
	}

	public UnivariatePolynomial<IntE> findAlgebraicInteger(ComplexNumber t, int maxDegree, IntE bound) {
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, maxDegree + 3);
		Integers z = Integers.z();
//		Rationals q = Rationals.q();
		IntE sqrtBound = z.one();// z.getInteger(bound.getValue().sqrt());
		List<Vector<Real>> basis = new ArrayList<>();
		for (int i = 0; i <= maxDegree; i++) {
			List<Real> vector = new ArrayList<>();
			for (int j = 0; j <= maxDegree; j++) {
				if (i == j) {
					IntE divisor = i == maxDegree ? sqrtBound : z.one();
					vector.add(r.getInteger(divisor));
				} else {
					vector.add(r.zero());
				}
			}
			ComplexNumber power = multiply(getInteger(bound), power(t, i));
			vector.add(r.getInteger(power.realPart().round()));
			vector.add(r.getInteger(power.complexPart().round()));
			basis.add(new Vector<>(vector));
		}
		Vector<Real> lll = space.latticeReduction(basis).get(0);
		List<IntE> coefficients = new ArrayList<>();
		for (int i = 0; i <= maxDegree; i++) {
			IntE divisor = i == maxDegree ? sqrtBound : z.one();
			coefficients.add(z.divide(lll.get(i + 1).round(), z.getInteger(divisor)));
		}
		UnivariatePolynomialRing<IntE> integerPolynomials = z.getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> polynomial = integerPolynomials.getPolynomial(coefficients);
		FactorizationResult<Polynomial<IntE>> factors = z.factorization(polynomial);
		UnivariatePolynomialRing<ComplexNumber> polynomials = getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> result = null;
		Real minValue = null;
		for (Polynomial<IntE> factor : factors.primeFactors()) {
			UnivariatePolynomial<ComplexNumber> embedded = polynomials.getEmbedding(factor, new MathMap<>() {
				@Override
				public ComplexNumber evaluate(IntE t) {
					return getInteger(t);
				}
			});
			Real value = value(polynomials.evaluate(embedded, t));
			if (minValue == null || value.compareTo(minValue) < 0) {
				minValue = value;
				result = integerPolynomials.toUnivariate(factor);
			}
		}
		return result;
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<ComplexNumber> t) {
		return t.degree() == 1;
	}

	@Override
	public FactorizationResult<Polynomial<ComplexNumber>> factorization(UnivariatePolynomial<ComplexNumber> t) {
		Map<Polynomial<ComplexNumber>, Integer> squareFree = getUnivariatePolynomialRing().squareFreeFactorization(t);
		SortedMap<Polynomial<ComplexNumber>, Integer> result = new TreeMap<>();
		for (Polynomial<ComplexNumber> squareFreeFactor : squareFree.keySet()) {
			for (Polynomial<ComplexNumber> factor : squareFreeFactorization(squareFreeFactor)) {
				result.put(factor, squareFree.get(squareFreeFactor));
			}
		}
		return new FactorizationResult<>(getUnivariatePolynomialRing().getEmbedding(t.leadingCoefficient()), result);
	}

	public List<Polynomial<ComplexNumber>> squareFreeFactorization(Polynomial<ComplexNumber> t) {
		Polynomial<ComplexNumber> f = t;
		List<Polynomial<ComplexNumber>> factors = new ArrayList<>();
		if (t.degree() <= 0) {
			return Collections.singletonList(t);
		}
		PolynomialRing<ComplexNumber> ring = getUnivariatePolynomialRing();
		while (f.degree() > 0) {
			ComplexNumber root = findRoot(f);
			Polynomial<ComplexNumber> factor = ring.subtract(ring.getVar(1), ring.getEmbedding(root));
			factors.add(factor);
			f = ring.divide(f, factor);
		}
		return factors;
	}

	@Override
	public boolean isCyclic() {
		return true;
	}

	@Override
	public GaloisGroup<Real, ComplexNumber, Complex> galoisGroup() {
		return new GaloisGroup<>(this, new FieldAutomorphism<>(this, new int[] { 1, 0 }));
	}

	@Override
	public FieldAutomorphism<Real, ComplexNumber, Complex> cyclicGenerator() {
		return new FieldAutomorphism<>(this, new int[] { 1, 0 });
	}

	@Override
	public ComplexNumber fromSmallDegreePolynomial(UnivariatePolynomial<Real> t) {
		return new ComplexNumber(t.univariateCoefficient(0), t.univariateCoefficient(1));
	}

	@Override
	public List<ComplexNumber> conjugates(ComplexNumber s) {
		List<ComplexNumber> conjugates = new ArrayList<>();
		conjugates.add(s);
		conjugates.add(conjugate(s));
		return conjugates;
	}

	@Override
	public Vector<Real> asVector(ComplexNumber s) {
		List<Real> l = new ArrayList<>();
		l.add(s.realPart);
		l.add(s.complexPart);
		return new Vector<>(l);
	}

	@Override
	public Complex makeExtension(UnivariatePolynomial<Real> minimalPolynomial) {
		return this;
	}

	@Override
	protected Complex asExtensionType() {
		return this;
	}
}
