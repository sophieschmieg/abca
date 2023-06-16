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
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.FloatingPointSet;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.ValueField;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.Identity;
import util.MiscAlgorithms;
import util.SingletonSortedMap;

public class Complex extends AbstractFieldExtension<Real, ComplexNumber, Complex>
		implements ValueField<ComplexNumber>, FloatingPointSet<ComplexNumber, Complex> {
	private Reals r;
	private Complex highPrecision;
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
		super(reals.getUnivariatePolynomialRing().getPolynomial(reals.one(), reals.zero(), reals.one()), reals, "i");
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

	private Complex highPrecision() {
		if (highPrecision == null) {
			highPrecision = withPrecision(precision() + 4);
		}
		return highPrecision;
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public boolean close(ComplexNumber t1, ComplexNumber t2) {
		ComplexNumber diff = subtract(t1, t2);
		Real norm = norm(diff);
		return norm.compareTo(r.getPowerOfTwo(-r.precision() + 1)) < 0;
//		int minRealExponent = Math.min(t1.realPart.exponent(), t2.realPart.exponent());
//		int minImgExponent = Math.min(t1.complexPart.exponent(), t2.complexPart.exponent());
//		int minExponent = Math.max(minRealExponent, minImgExponent);
//		int precision = r.precision();
//		return r.abs(r.subtract(t1.realPart, t2.realPart)).compareTo(r.getPowerOfTwo(-precision + minExponent + 1)) <= 0
//				&& r.abs(r.subtract(t1.complexPart, t2.complexPart))
//						.compareTo(r.getPowerOfTwo(-precision + minExponent + 1)) <= 0;
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

	public boolean isReal(ComplexNumber t) {
		return r.abs(t.complexPart()).compareTo(r.getPowerOfTwo(-r.precision() + t.realPart().exponent())) < 0;
	}

	public Real asReal(ComplexNumber t) {
		if (!isReal(t)) {
			throw new ArithmeticException("Not a real number!");
		}
		return t.realPart();
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
		Real denom = highPrecision().norm(divisor);
		return getEmbedding(highPrecision().multiply(highPrecision().getEmbedding(highPrecision().r.inverse(denom)),
				dividend, highPrecision().conjugate(divisor)));
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

	@Override
	public Reals getReals() {
		return r;
	}

	@Override
	public Reals getBaseField() {
		return r;
	}

	public ComplexNumber exp(ComplexNumber t) {
		return getEmbedding(highPrecision().calculateExp(t, this));
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
		UnivariatePolynomial<ComplexNumber> p = getUnivariatePolynomialRing().toUnivariate(t);
		int scale = 0;
		int lcScale = value(t.leadingCoefficient()).exponent() - 1;
		for (int i = 0; i < t.degree(); i++) {
			int coeffScale = value(p.univariateCoefficient(i)).exponent() - lcScale;
			coeffScale = MiscAlgorithms.DivRoundUp(coeffScale, t.degree() - i);
			scale = Math.max(coeffScale, scale);
		}
		Complex highPrecision = withPrecision(precision() + scale + 2);
		UnivariatePolynomialRing<ComplexNumber> ring = highPrecision.getUnivariatePolynomialRing();
		UnivariatePolynomial<ComplexNumber> derivative = ring.derivative(t);
		int maxNumIterations = 1000;
		while (true) {
			ComplexNumber x = highPrecision.getRandomElement();
			ComplexNumber prevX = highPrecision.zero();
			for (int i = 0; i < maxNumIterations; i++) {
				prevX = x;
				x = highPrecision.subtract(x, highPrecision.divide(ring.evaluate(t, x), ring.evaluate(derivative, x)));
				if (ring.evaluate(t, x).equals(highPrecision.zero()) || highPrecision.close(x, prevX)) {
					return getEmbedding(x);
				}
			}
		}
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<ComplexNumber> t) {
		return t.degree() == 1;
	}

	@Override
	public FactorizationResult<Polynomial<ComplexNumber>, ComplexNumber> factorization(
			UnivariatePolynomial<ComplexNumber> t) {
		if (t.degree() == 1) {
			return new FactorizationResult<>(t.leadingCoefficient(),
					SingletonSortedMap.map(getUnivariatePolynomialRing().normalize(t), 1));
		}
		Rationals q = Rationals.q();
		NumberField gauss = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
		UnivariatePolynomial<NFE> rationalPolynomial = gauss.getUnivariatePolynomialRing().getEmbedding(t,
				new MathMap<>() {
					@Override
					public NFE evaluate(ComplexNumber t) {
						return gauss.add(gauss.getEmbedding(r.roundToFraction(t.realPart(), r.precision())),
								gauss.scalarMultiply(r.roundToFraction(t.complexPart(), r.precision()), gauss.alpha()));
					}
				});
		FactorizationResult<Polynomial<NFE>, NFE> rationalSquareFree = gauss.getUnivariatePolynomialRing()
				.squareFreeFactorization(rationalPolynomial);
		FactorizationResult<Polynomial<ComplexNumber>, ComplexNumber> squareFree;
		if (rationalSquareFree.isIrreducible()) {
			squareFree = new FactorizationResult<Polynomial<ComplexNumber>, ComplexNumber>(t.leadingCoefficient(),
					SingletonSortedMap.map(getUnivariatePolynomialRing().normalize(t), 1));
		} else {
			SortedMap<Polynomial<ComplexNumber>, Integer> squareFreeResult = new TreeMap<>();
			for (Polynomial<NFE> rationalSquareFreeFactor : rationalSquareFree.primeFactors()) {
				squareFreeResult
						.put(getUnivariatePolynomialRing().getEmbedding(rationalSquareFreeFactor, new MathMap<>() {
							@Override
							public ComplexNumber evaluate(NFE t) {
								return getNumber(r.getEmbedding(t.asPolynomial().univariateCoefficient(0)),
										r.getEmbedding(t.asPolynomial().univariateCoefficient(1)));
							}
						}), rationalSquareFree.multiplicity(rationalSquareFreeFactor));
			}
			squareFree = new FactorizationResult<Polynomial<ComplexNumber>, ComplexNumber>(t.leadingCoefficient(),
					squareFreeResult);
		}
		SortedMap<Polynomial<ComplexNumber>, Integer> result = new TreeMap<>();
		for (Polynomial<ComplexNumber> squareFreeFactor : squareFree.primeFactors()) {
			for (Polynomial<ComplexNumber> factor : squareFreeFactorization(squareFreeFactor)) {
				result.put(factor, squareFree.multiplicity(squareFreeFactor));
			}
		}
		return new FactorizationResult<>(t.leadingCoefficient(), result);
	}

	public List<Polynomial<ComplexNumber>> squareFreeFactorization(Polynomial<ComplexNumber> t) {
		Polynomial<ComplexNumber> f = t;
		List<Polynomial<ComplexNumber>> factors = new ArrayList<>();
		if (t.degree() <= 0) {
			return Collections.singletonList(t);
		}
		PolynomialRing<ComplexNumber> ring = highPrecision().getUnivariatePolynomialRing();
		UnivariatePolynomialRing<ComplexNumber> r = getUnivariatePolynomialRing();
		while (f.degree() > 0) {
			ComplexNumber root = highPrecision().findRoot(f);
			Polynomial<ComplexNumber> factor = ring.subtract(ring.getVar(1), ring.getEmbedding(root));
			factors.add(r.getEmbedding(factor, new MathMap<>() {
				@Override
				public ComplexNumber evaluate(ComplexNumber t) {
					return getEmbedding(t);
				}
			}));
			f = ring.divide(f, factor);
		}
		return factors;
	}

	@Override
	public boolean isSubModuleMember(MatrixModule<ComplexNumber> module, Matrix<ComplexNumber> m,
			Vector<ComplexNumber> b) {
		FiniteComplexVectorSpace space = new FiniteComplexVectorSpace(this, b.dimension());
		return space.isSubModuleMemberModule(module, m, b);
	}

	@Override
	public Vector<ComplexNumber> asSubModuleMember(MatrixModule<ComplexNumber> module, Matrix<ComplexNumber> m,
			Vector<ComplexNumber> b) {
		FiniteComplexVectorSpace space = new FiniteComplexVectorSpace(this, b.dimension());
		return space.asSubModuleMemberModule(module, m, b);
	}

	@Override
	public List<Vector<ComplexNumber>> syzygyProblem(MatrixModule<ComplexNumber> module, Matrix<ComplexNumber> m) {
		FiniteComplexVectorSpace space = new FiniteComplexVectorSpace(this, m.rows());
		return space.syzygyProblemModule(module, m);
	}

	@Override
	public List<Vector<ComplexNumber>> simplifySubModuleGenerators(MatrixModule<ComplexNumber> module,
			Matrix<ComplexNumber> m) {
		FiniteComplexVectorSpace space = new FiniteComplexVectorSpace(this, m.columns());
		return space.simplifySubModuleGeneratorsModule(module, m);
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
