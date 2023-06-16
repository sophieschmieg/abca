package fields.floatingpoint;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.FloatingPointSet;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.TotalOrder;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.ValueField;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class Reals extends AbstractField<Real> implements ValueField<Real>, FloatingPointSet<Real, Reals>, TotalOrder<Real> {
	private final int precision;
	private Reals highPrecision;
	private Real pi;
	private Real e;
	private Real zero;
	private Real one;
	private static Map<Integer, Reals> reals = new TreeMap<>();

	public class Real extends AbstractElement<Real> {
		private boolean positive;
		private BigInteger value;
		private int scale;

		private Real(BigInteger value, int scale, boolean positive) {
			if (value.compareTo(BigInteger.ZERO) < 0) {
				value = value.negate();
				positive = !positive;
			}
			if (value.bitLength() > precision) {
				boolean roundUp = value.testBit(value.bitLength() - precision - 1);
				scale += value.bitLength() - precision;
				value = value.shiftRight(value.bitLength() - precision);
				if (roundUp) {
					value = value.add(BigInteger.ONE);
				}
			}
			if (value.equals(BigInteger.ZERO)) {
				this.positive = true;
				this.value = value;
				this.scale = 0;
				return;
			}
			while (!value.testBit(0)) {
				value = value.shiftRight(1);
				scale++;
			}
			this.positive = positive;
			this.value = value;
			this.scale = scale;
		}

		@Override
		public int compareTo(Real o) {
			if (!this.positive && o.positive) {
				return -1;
			}
			if (this.positive && !o.positive) {
				return 1;
			}
			if (this.value.equals(BigInteger.ZERO) && o.value.equals(BigInteger.ZERO)) {
				return 0;
			}
			if (this.value.equals(BigInteger.ZERO) && o.positive) {
				return -1;
			}
			if (o.value.equals(BigInteger.ZERO) && this.positive) {
				return 1;
			}
			int expDiff = this.exponent() - o.exponent();
			if (expDiff != 0) {
				return positive ? expDiff : -expDiff;
			}
			BigInteger thisValue = value.shiftLeft(Math.max(0, this.scale - o.scale));
			BigInteger otherValue = o.value.shiftLeft(Math.max(0, o.scale - this.scale));
			int cmp = thisValue.compareTo(otherValue);
			return positive ? cmp : -cmp;
		}

		public IntE roundDown() {
			Integers z = Integers.z();
			if (!positive) {
				return z.negative(Reals.this.negative(this).roundUp());
			}
			return z.getInteger(value.shiftLeft(scale));
		}

		public IntE roundUp() {
			Integers z = Integers.z();
			if (!positive) {
				return z.negative(Reals.this.negative(this).roundDown());
			}
			if (scale >= 0) {
				return z.getInteger(value.shiftLeft(scale));
			}
			return z.add(z.getInteger(value.shiftLeft(scale)), z.one());
		}

		public IntE round() {
			if (scale >= 0) {
				return positive ? roundDown() : roundUp();
			}
			if (value.testBit(-scale - 1)) {
				return positive ? roundUp() : roundDown();
			}
			return positive ? roundDown() : roundUp();
		}

		public Real fractionalPart() {
			return subtract(this, getInteger(this.roundDown()));
		}

		public int exponent() {
			return scale + value.bitLength();
		}

		public String toString() {
			StringBuilder buf = new StringBuilder();
			if (!positive) {
				buf.append("-");
			}
			buf.append(value.shiftLeft(scale));
			buf.append(".");
			BigInteger fractionalPart = value.subtract(value.shiftLeft(scale).shiftRight(scale));
			if (fractionalPart.equals(BigInteger.ZERO)) {
				buf.append(BigInteger.ZERO);
				return buf.toString();
			}
			fractionalPart = fractionalPart.multiply(BigInteger.TEN.pow(-scale));
			fractionalPart = fractionalPart.shiftRight(-scale);
			String fractionalPartString = fractionalPart.toString();
			for (int i = 0; i < -scale - fractionalPartString.length(); i++) {
				buf.append("0");
			}
			buf.append(fractionalPartString.substring(0, Math.min(fractionalPartString.length(),
					(int) Math.round(Math.floor(precision * Math.log(2) / Math.log(10))))));
			return buf.toString();
		}

		public double doubleValue() {
			Reals r = withPrecision(64);
			Real t = r.getEmbedding(this);
			if (t.scale > 1024) {
				return t.positive ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
			}
			if (t.scale < -1024) {
				return t.positive ? 0.0 : -0.0;
			}
			double abs = t.value.doubleValue() * Math.pow(2.0, t.scale);
			return t.positive ? abs : -abs;
		}
	}

	private Reals(int precision) {
		if (precision <= 0) {
			throw new ArithmeticException("Precision too low");
		}
		this.precision = precision;
		this.zero = new Real(BigInteger.ZERO, 0, true);
		this.one = new Real(BigInteger.ONE, 0, true);
	}

	public static Reals r(int precision) {
		if (!reals.containsKey(precision)) {
			reals.put(precision, new Reals(precision));
		}
		return reals.get(precision);
	}

	public String toString() {
		return "R";
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public boolean close(Real t1, Real t2) {
		int minExponent = Math.min(t1.exponent(), t2.exponent());
		return abs(subtract(t1, t2)).compareTo(new Real(BigInteger.ONE, -precision + minExponent + 1, true)) <= 0;
	}

	@Override
	public Real roundToFixedPoint(Real t, int precision) {
		BigInteger value = t.value.shiftLeft(precision);
		value = value.shiftLeft(t.scale);
		return new Real(value, -precision, t.positive);
	}

	@Override
	public int precision() {
		return precision;
	}

	@Override
	public Reals withPrecision(int precision) {
		return r(precision);
	}

	@Override
	public Real zero() {
		return zero;
	}

	@Override
	public Real one() {
		return one;
	}

	public Real fromString(String string) {
		boolean positive = string.charAt(0) != '-';
		int dotIndex = string.indexOf(".");
		String intValueString;
		String fractionalValueString = "0";
		if (dotIndex < 0) {
			intValueString = string.substring(positive ? 0 : 1);
		} else {
			intValueString = string.substring(positive ? 0 : 1, dotIndex);
			if (dotIndex != string.length()) {
				fractionalValueString = string.substring(dotIndex + 1);
			}
		}
		if (intValueString.length() == 0) {
			intValueString = "0";
		}
		int fractionalValueLength = fractionalValueString.length();
		BigInteger intValue = new BigInteger(intValueString);
		BigInteger fractionalValue = new BigInteger(fractionalValueString);
		BigInteger scaledValue = intValue.multiply(BigInteger.TEN.pow(fractionalValueLength)).add(fractionalValue);
		Real value = divide(getInteger(scaledValue), getInteger(BigInteger.TEN.pow(fractionalValueLength)));
		return positive ? value : negative(value);
	}

	@Override
	public Real getInteger(BigInteger t) {
		return new Real(t, 0, true);
	}

	public Real getPowerOfTwo(int power) {
		return new Real(BigInteger.ONE, power, true);
	}

	public Real getEmbedding(Real t) {
		return new Real(t.value, t.scale, t.positive);
	}

	public Real getEmbedding(IntE t) {
		return getInteger(t.getValue());
	}

	public Real getEmbedding(Fraction t) {
		return divide(getEmbedding(t.getNumerator()), getEmbedding(t.getDenominator()));
	}

	@Override
	public Extension<Real, Real, ComplexNumber, Complex> getExtension(UnivariatePolynomial<Real> minimalPolynomial) {
		if (minimalPolynomial.degree() == 1) {
			throw new ArithmeticException("no degree 1 reals!");
		}
		if (minimalPolynomial.degree() != 2) {
			throw new ArithmeticException("not irreducible!");
		}
		Real discriminant = subtract(
				multiply(minimalPolynomial.univariateCoefficient(1), minimalPolynomial.univariateCoefficient(1)),
				multiply(4, minimalPolynomial.univariateCoefficient(2), minimalPolynomial.univariateCoefficient(0)));
		if (discriminant.compareTo(zero()) >= 0) {
			throw new ArithmeticException("not irreducible!");
		}
		Complex c = Complex.c(this.precision);
		return new Extension<>(c, this, c.getEmbeddingMap(), c.asVectorMap());
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	private Reals highPrecision() {
		if (highPrecision == null) {
			highPrecision = withPrecision(precision + 2);
		}
		return highPrecision;
	}

	@Override
	public Real add(Real t1, Real t2) {
		int minScale = Math.min(t1.scale, t2.scale);
		BigInteger value1Adjusted = t1.value.shiftLeft(t1.scale - minScale);
		BigInteger value2Adjusted = t2.value.shiftLeft(t2.scale - minScale);
		if (!t1.positive) {
			value1Adjusted = value1Adjusted.negate();
		}
		if (!t2.positive) {
			value2Adjusted = value2Adjusted.negate();
		}
		return new Real(value1Adjusted.add(value2Adjusted), minScale, true);
	}

	@Override
	public Real negative(Real t) {
		return new Real(t.value, t.scale, !t.positive);
	}

	@Override
	public Real multiply(Real t1, Real t2) {
		return new Real(t1.value.multiply(t2.value), t1.scale + t2.scale, !(t1.positive ^ t2.positive));
	}

	@Override
	public Real inverse(Real t) {
		if (t.equals(zero())) {
			throw new ArithmeticException("Division by 0");
		}
		if (t.compareTo(zero()) < 0) {
			return negative(inverse(negative(t)));
		}
		return getEmbedding(highPrecision().calculateInverse(t));
	}

	private Real calculateInverse(Real t) {
		Real x = new Real(BigInteger.ONE, -t.exponent(), true);
		Real prev;
		Real two = getInteger(2);
		do {
			prev = x;
			x = multiply(x, subtract(two, multiply(x, t)));
		} while (!close(prev, x));
		return x;
	}

	public Real abs(Real t) {
		return new Real(t.value, t.scale, true);
	}

	public Real exp(Real t) {
		Complex c = Complex.c(precision);
		return c.exp(c.getEmbedding(t)).realPart();
	}

	public Real log(Real t) {
		if (t.compareTo(zero()) <= 0) {
			throw new ArithmeticException("log of non positive number!");
		}
		return getEmbedding(highPrecision().calculateLog(t));
	}

	private Real calculateLog(Real t) {
		int k = 0;
		while (t.compareTo(getDouble(0.875)) < 0) {
			k--;
			t = multiply(t, e());
		}
		while (t.compareTo(getDouble(1.125)) > 0) {
			k++;
			t = divide(t, e());
		}
		Real x = subtract(one(), t);
		Real result = zero();
		Real prevResult;
		Real power = one();
		int n = 1;
		do {
			prevResult = result;
			power = multiply(power, x);
			result = add(result, divide(power, getInteger(n)));
			n++;
		} while (!result.equals(prevResult));
		return subtract(getInteger(k), result);
	}

	public MathMap<Real, Real> exp() {
		return new MathMap<>() {
			@Override
			public Real evaluate(Real t) {
				return exp(t);
			}
		};
	}

	public Real e() {
		if (e == null) {
			e = exp(one());
		}
		return e;
	}

	public Real pi() {
		if (pi == null) {
			pi = getEmbedding(highPrecision().calculatePi());
		}
		return pi;
	}

	private Real calculatePi() {
		Real sqrtTwo = positiveSqrt(getInteger(2));
		BigInteger factorial = BigInteger.ONE;
		BigInteger fourFactorial = BigInteger.ONE;
		Real sum = zero();
		Real prevSum;
		int k = 0;
		do {
			prevSum = sum;
			sum = add(sum, divide(multiply(fourFactorial, add(getInteger(1103), getInteger(26390 * k))),
					power(multiply(getInteger(factorial), power(getInteger(396), k)), 4)));
			k++;
			factorial = factorial.multiply(BigInteger.valueOf(k));
			for (int j = 0; j < 4; j++) {
				fourFactorial = fourFactorial.multiply(BigInteger.valueOf(4 * k - j));
			}
		} while (!sum.equals(prevSum));
		Real pi = divide(getInteger(9801), multiply(2, sqrtTwo, sum));
		return pi;
	}

	public Real positiveSqrt(Real t) {
		return positiveRoot(t, 2);
	}

	public Real positiveRoot(Real t, int degree) {
		if (t.equals(zero())) {
			return zero();
		}
		return multiply(t, power(inverseRoot(t, degree), degree-1));
//		if (t.compareTo(zero()) < 0) {
//			throw new ArithmeticException("No positive root");
//		}
//		return getEmbedding(highPrecision().calculateRoot(t, degree));
	}

	private Real calculateRoot(Real t, int degree) {
		Real result = getPowerOfTwo(MiscAlgorithms.DivRoundUp(t.exponent(), degree));
		// Real prevResult;
		do {
			// prevResult = result;
			result = divide(add(multiply(degree - 1, result), divide(t, power(result, degree - 1))),
					getInteger(degree));
		} while (!close(power(result, degree), t));
		return result;
	}

	public Real inverseSqrt(Real t) {
		return inverseRoot(t, 2);
	}

	public Real inverseRoot(Real t, int degree) {
		if (t.compareTo(zero()) <= 0) {
			throw new ArithmeticException("No positive root");
		}
		return getEmbedding(highPrecision().calculateInverseRoot(t, degree));
	}

	private Real calculateInverseRoot(Real t, int degree) {
		Real rootEstimate = getPowerOfTwo(MiscAlgorithms.DivRoundUp(t.exponent(), degree));
		for (int i = 0; i < 2; i++) {
			rootEstimate = divide(add(multiply(degree - 1, rootEstimate), divide(t, power(rootEstimate, degree - 1))),
					getInteger(degree));
		}
		Real result = inverse(rootEstimate);
		Real prevResult;
		Real coeff1 = divide(getInteger(degree + 1), getInteger(degree));
		Real coeff2 = divide(negative(t), getInteger(degree));
		do {
			prevResult = result;
			result = multiply(result, add(coeff1, multiply(coeff2, power(result, degree))));
		} while (!close(result, prevResult));
		return result;
	}
	
	public Real sin(Real angle) {
		return cosineAndSine(angle).get(2);
	}
	
	public Real cos(Real angle) {
		return cosineAndSine(angle).get(1);
	}
	
	public Vector<Real> cosineAndSine(Real angle) {
		Complex c = Complex.c(precision);
		ComplexNumber result = c.exp(c.getNumber(zero(), angle));
		return new Vector<>(result.realPart(), result.complexPart());
	}

	public Real arctan(Real t) {
		/*
		 * if (t.compareTo(one()) > 0) { return subtract(divide(pi(), getInteger(2)),
		 * doublePrecision().calculateArctan(doublePrecision().inverse(t))); } else if
		 * (t.compareTo(negative(one())) < 0) { return subtract(divide(negative(pi()),
		 * getInteger(2)),
		 * doublePrecision().calculateArctan(doublePrecision().inverse(t))); }
		 */
		return getEmbedding(highPrecision().calculateArctan(t));
	}

	private Real calculateArctan(Real t) {
		if (abs(t).compareTo(new Real(BigInteger.ONE, -8, true)) > 0) {
			Real x = divide(t, add(one(), positiveSqrt(add(one(), multiply(t, t)))));
			return multiply(2, calculateArctan(x));
		}
		Real result = t;
		Real power = t;
		Real prev;
		int i = 1;
		do {
			prev = result;
			power = multiply(power, t, t);
			Real summand = divide(power, getInteger(2 * i + 1));
			if (i % 2 == 0) {
				result = add(result, summand);
			} else {
				result = subtract(result, summand);
			}
			i++;
		} while (!close(result, prev));
		return result;
	}

	public Real arctan2(Real t1, Real t2) {
		Real piOverTwo = divide(pi(), getInteger(2));
		if (t2.equals(zero())) {
			if (t1.positive) {
				return piOverTwo;
			} else {
				return negative(piOverTwo);
			}
		}
		Real arctan = arctan(divide(t1, t2));
		if (t2.positive) {
			return arctan;
		}
		if (t1.positive) {
			return add(arctan, piOverTwo);
		}
		return subtract(arctan, piOverTwo);
	}

	public Real findZero(MathMap<Real, Real> map, Real start, Real end) {
		if (start.compareTo(end) > 0) {
			return findZero(map, end, start);
		}
		Real atStart = map.evaluate(start);
		Real atEnd = map.evaluate(end);
		if (atStart.positive == atEnd.positive) {
			throw new ArithmeticException("Cannot use intermediate value theorem");
		}
		while (!atStart.equals(zero())) {
			Real mid2 = add(start, end);
			Real mid = new Real(mid2.value, mid2.scale - 1, mid2.positive);
			Real atMid = map.evaluate(mid);
			if (atStart.positive == atMid.positive) {
				start = mid;
				atStart = atMid;
			} else {
				end = mid;
				atEnd = atMid;
			}
		}
		return start;
	}

	public Real invertMonotonic(MathMap<Real, Real> map, Real t) {
		Real start = negative(one());
		Real end = one();
		boolean startSmaller;
		boolean endSmaller;
		do {
			start = multiply(2, start);
			end = multiply(2, end);
			Real atStart = map.evaluate(start);
			Real atEnd = map.evaluate(end);
			startSmaller = atStart.compareTo(t) < 0;
			endSmaller = atEnd.compareTo(t) < 0;
		} while (!(startSmaller ^ endSmaller));
		return invertMonotonic(map, t, start, end);
	}

	public Real invertMonotonic(MathMap<Real, Real> map, Real at, Real start, Real end) {
		return findZero(new MathMap<>() {
			@Override
			public Real evaluate(Real t) {
				return subtract(map.evaluate(t), at);
			}
		}, start, end);
	}

	public MathMap<Real, Real> invertMonotonicMap(MathMap<Real, Real> map) {
		return new MathMap<>() {
			@Override
			public Real evaluate(Real t) {
				return invertMonotonic(map, t);
			}
		};
	}

	public Real getDouble(double t) {
		if (t < 0) {
			return negative(getDouble(-t));
		}
		int exponent = Math.getExponent(t);
		int scale = exponent - 53;
		BigInteger value = BigInteger.valueOf(Math.round(t * Math.pow(2.0, -scale)));
		return new Real(value, scale, true);
	}

	@Override
	public Real getRandomElement() {
		return getDouble(new Random().nextGaussian());
	}

	public Real getRandomElement(Real from, Real to) {
		return add(multiply(getUniformRandomElement(), subtract(from, to)), from);
	}

	public Real getUniformRandomElement() {
		return getDouble(new Random().nextDouble());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Real> iterator() {
		throw new InfinityException();
	}

	@Override
	public Real value(Real t) {
		return abs(t);
	}

	@Override
	public Reals getReals() {
		return this;
	}

	public Iterator<IntE> continuedFraction(Real t) {
		return MiscAlgorithms.continuedFraction(this, t, new MathMap<>() {

			@Override
			public IntE evaluate(Real t) {
				return t.roundDown();
			}
		});
	}

	public Iterator<Fraction> continuedFractionApproximation(Real t) {
		return MiscAlgorithms.continuedFractionApproximation(this, t, new MathMap<>() {

			@Override
			public IntE evaluate(Real t) {
				return t.roundDown();
			}
		});
	}

	public Fraction roundToFraction(Real t, int precision) {
		if (close(t, zero())) {
			return Rationals.q().zero();
		}
		Iterator<Fraction> it = continuedFractionApproximation(t);
		Real epsilon = power(getInteger(2), -precision);
		while (it.hasNext()) {
			Fraction approximation = it.next();
			if (abs(subtract(t, getFraction(approximation))).compareTo(epsilon) < 0) {
				return approximation;
			}
		}
		throw new ArithmeticException("Should not reach end of continued fraction before finding exact match!");
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<Real> t) {
		if (t.degree() == 1) {
			return true;
		}
		if (t.degree() > 2) {
			return false;
		}
		if (t.degree() <= 0) {
			throw new ArithmeticException("t is unit or 0");
		}
		Real discriminant = subtract(multiply(t.univariateCoefficient(1), t.univariateCoefficient(1)),
				multiply(4, t.univariateCoefficient(2), t.univariateCoefficient(0)));
		return discriminant.compareTo(zero()) < 0;
	}

	@Override
	public FactorizationResult<Polynomial<Real>, Real> factorization(UnivariatePolynomial<Real> t) {
		Complex c = Complex.c(precision);
		Map<ComplexNumber, Integer> complexRoots = c
				.roots(c.getUnivariatePolynomialRing().getEmbedding(t, c.getEmbeddingMap()));
		SortedMap<Polynomial<Real>, Integer> factors = new TreeMap<>();
		PolynomialRing<Real> ring = getUnivariatePolynomialRing();
		for (ComplexNumber root : complexRoots.keySet()) {
			int mult = complexRoots.get(root);
			if (root.complexPart().equals(zero())) {
				factors.put(ring.subtract(ring.getVar(1), ring.getEmbedding(root.realPart())), mult);
			} else if (root.complexPart().compareTo(zero()) > 0) {
				factors.put(ring.add(ring.getVarPower(1, 2),
						ring.getEmbedding(multiply(-1, c.trace(root)), new int[] { 1 }),
						ring.getEmbedding(c.norm(root))), mult);
			}
		}
		return new FactorizationResult<>(t.leadingCoefficient(), factors);
	}
	
	@Override
	public boolean isSubModuleMember(MatrixModule<Real> module, Matrix<Real> m,
			Vector<Real> b) {
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(this, b.dimension());
		return space.isSubModuleMemberModule(module, m, b);
	}

	@Override
	public Vector<Real> asSubModuleMember(MatrixModule<Real> module, Matrix<Real> m,
			Vector<Real> b) {
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(this, b.dimension());
		return space.asSubModuleMemberModule(module, m, b);
	}

	@Override
	public List<Vector<Real>> syzygyProblem(MatrixModule<Real> module, Matrix<Real> m) {
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(this, m.columns());
		return space.syzygyProblemModule(module, m);
	}

	@Override
	public List<Vector<Real>> simplifySubModuleGenerators(MatrixModule<Real> module,
			Matrix<Real> m) {
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(this, m.rows());
		return space.simplifySubModuleGeneratorsModule(module, m);
	}
}
