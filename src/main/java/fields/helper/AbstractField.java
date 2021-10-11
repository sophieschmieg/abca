package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.LocalRing;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.LocalRingImplementation;
import fields.local.TrivialDiscreteValuationField;
import fields.local.Value;
import fields.polynomials.GenericUnivariatePolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.TrivialPivotStrategy;
import util.ConstantMap;
import util.Identity;
import util.Pair;
import util.SingletonSortedMap;

public abstract class AbstractField<T extends Element<T>> implements Field<T>, Ring<T>, DedekindRing<T, T, T>, LocalRing<T, T, T> {
	private GenericUnivariatePolynomialRing<T> ring = new GenericUnivariatePolynomialRing<T>(this);

	@Override
	public Group<T> getAdditiveGroup() {
		return new Group<T>() {

			@Override
			public Exactness exactness() {
				return AbstractField.this.exactness();
			}

			@Override
			public T getRandomElement() {
				return this.getRandomElement();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractField.this.iterator();
			}

			@Override
			public T neutral() {
				return AbstractField.this.zero();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.negative(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.add(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return this.isFinite();
			}

		};
	}

	@Override
	public Group<T> getMultiplicativeGroup() {
		return new Group<T>() {

			@Override
			public Exactness exactness() {
				return AbstractField.this.exactness();
			}

			@Override
			public T getRandomElement() {
				T t;
				do {
					t = this.getRandomElement();
				} while (t.equals(AbstractField.this.zero()));
				return t;
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements().subtract(BigInteger.ONE);
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractField.this.getNonZeroElements().iterator();
			}

			@Override
			public T neutral() {
				return AbstractField.this.one();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.inverse(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.multiply(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractField.this.isFinite();
			}
		};
	}

	@Override
	public Extension<T, ?, ?, ?> getExtension(UnivariatePolynomial<T> minimalPolynomial) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public T add(T t1, T t2, T t3) {
		return this.add(this.add(t1, t2), t3);
	}

	@Override
	public T add(T t1, T t2, T t3, T t4) {
		return this.add(this.add(t1, t2), this.add(t3, t4));
	}

	@Override
	public T subtract(T minuend, T subtrahend) {
		return this.add(minuend, this.negative(subtrahend));
	}

	@Override
	public T getInteger(int n) {
		return getInteger(BigInteger.valueOf(n));
	}

	@Override
	public T getInteger(BigInteger n) {
		return multiply(n, this.one());
	}

	@Override
	public T getInteger(IntE n) {
		return getInteger(n.getValue());
	}

	@Override
	public T getFraction(Fraction t) {
		return divide(getInteger(t.getNumerator()), getInteger(t.getDenominator()));
	}

	@Override
	public T multiply(int n, T t) {
		return this.multiply(BigInteger.valueOf(n), t);
	}

	@Override
	public T multiply(BigInteger n, T t) {
		T dbl;
		if (n.signum() < 0) {
			n = n.negate();
			dbl = this.negative(t);
		} else {
			dbl = t;
		}
		T result = zero();
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.add(result, dbl);
			dbl = this.add(dbl, dbl);
		}
		return result;
	}

	@Override
	public T multiply(IntE n, T t) {
		return this.multiply(n.getValue(), t);
	}

	@Override
	public T multiply(T t1, T t2, T t3) {
		return multiply(multiply(t1, t2), t3);
	}

	@Override
	public T multiply(int n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}

	@Override
	public T multiply(BigInteger n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}

	@Override
	public T multiply(IntE n, T t1, T t2) {
		return this.multiply(n.getValue(), this.multiply(t1, t2));
	}

	@Override
	public T multiply(int n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}

	@Override
	public T multiply(BigInteger n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}

	@Override
	public T multiply(IntE n, T t1, T t2, T t3) {
		return this.multiply(n.getValue(), this.multiply(t1, t2, t3));
	}

	@Override
	public T divide(T dividend, T divisor) {
		return multiply(dividend, inverse(divisor));
	}

	@Override
	public T power(T t, int n) {
		return this.power(t, BigInteger.valueOf(n));
	}

	@Override
	public T power(T t, BigInteger n) {
		T square;
		if (n.signum() < 0) {
			n = n.negate();
			square = this.inverse(t);
		} else {
			square = t;
		}
		T result = one();
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.multiply(result, square);
			square = this.multiply(square, square);
		}
		return result;
	}

	public BigInteger discreteLogarithm(T base, T t) {
		if (!this.isFinite()) {
			throw new UnsupportedOperationException("not implemented");
		}
		throw new UnsupportedOperationException("not implemented");
	}

	public Iterable<T> getNonZeroElements() throws InfinityException {
		return new Iterable<T>() {
			private Iterable<T> elements = AbstractField.this;

			@Override
			public Iterator<T> iterator() {
				return new Iterator<T>() {
					private Iterator<T> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public T next() {
						T result = null;
						do {
							if (it.hasNext())
								result = it.next();
						} while (result.equals(AbstractField.this.zero()));
						return result;
					}
				};
			}

		};
	}

	@Override
	public boolean isUnit(T t) {
		return !t.equals(this.zero());
	}

	@Override
	public T projectToUnit(T t) {
		if (t.equals(zero())) {
			return one();
		}
		return t;
	}

	@Override
	public T upToUnit(T t) {
		if (t.equals(zero())) {
			return zero();
		}
		return one();
	}

	@Override
	public T numerator(T t) {
		return t;
	}

	@Override
	public T denominator(T t) {
		return one();
	}

	@Override
	public boolean hasRoots(Polynomial<T> t) {
		FactorizationResult<Polynomial<T>, T> factors = factorization(ring.toUnivariate(t));
		for (Polynomial<T> factor : factors.primeFactors()) {
			if (factor.degree() == 1) {
				return true;
			}
		}
		return false;
	}

	@Override
	public Map<T, Integer> roots(Polynomial<T> t) {
		Map<T, Integer> result = new TreeMap<>();
		Monomial constant = ring.getMonomial(new int[] { 0 });
		Monomial linear = ring.getMonomial(new int[] { 1 });
		FactorizationResult<Polynomial<T>, T> factors = factorization(ring.toUnivariate(t));
		for (Polynomial<T> factor : factors.primeFactors()) {
			if (factor.degree() == 1) {
				result.put(divide(negative(factor.coefficient(constant)), factor.coefficient(linear)),
						factors.multiplicity(factor));
			}
		}
		return result;
	}

	@Override
	public boolean hasRoot(T t, int n) {
		if (t.equals(this.zero())) {
			return true;
		}
		Monomial m = ring.getMonomial(new int[] { n });
		Monomial c = ring.getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = ring.getPolynomial(coeff);
		return this.hasRoots(p);
	}

	@Override
	public Map<T, Integer> roots(T t, int n) {
		if (t.equals(this.zero())) {
			return Collections.singletonMap(this.zero(), n);
		}
		Monomial m = ring.getMonomial(new int[] { n });
		Monomial c = ring.getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = ring.getPolynomial(coeff);
		return this.roots(p);
	}

	@Override
	public boolean hasSqrt(T t) {
		return hasRoot(t, 2);
	}

	@Override
	public Map<T, Integer> sqrt(T t) {
		return roots(t, 2);
	}

	@Override
	public T primitiveRoot() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
		return true;
	}

	@Override
	public Value valuation(T t, Ideal<T> maximalIdeal) {
		return t.equals(zero()) ? Value.INFINITY : Value.ZERO;
	}

	@Override
	public FieldOfFractionsResult<T, T> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, this, new Identity<>(), new Identity<>(), new ConstantMap<>(one()),
				new Identity<>());
	}

	@Override
	public LocalizeResult<T, T, T, T> localizeAtIdeal(Ideal<T> primeIdeal) {
		if (!primeIdeal.isPrime()) {
			throw new ArithmeticException("Not a prime ideal!");
		}
		return new LocalizeResult<>(this, primeIdeal, this, new Identity<>(), new Identity<>(),
				new ConstantMap<>(one()), new Identity<>());
	}
	
	@Override
	public FieldIdeal<T> maximalIdeal() {
		return getZeroIdeal();
	}
	
	@Override
	public Field<T> reduction() {
		return this;
	}
	
	@Override
	public T reduce(T t) {
		return t;
	}
	
	@Override
	public T lift(T t) {
		return t;
	}

	@Override
	public boolean isInteger(T t) {
		return true;
	}

	@Override
	public T asInteger(T t) {
		return t;
	}

	@Override
	public Field<T> reduction(Ideal<T> maximalIdeal) {
		return this;
	}

	@Override
	public T reduce(T t, Ideal<T> maximalIdeal) {
		return t;
	}

	@Override
	public T lift(T s, Ideal<T> maximalIdeal) {
		return s;
	}

	@Override
	public DiscreteValuationRing<T, T> localize(Ideal<T> maximalIdeal) {
		return new LocalRingImplementation<>(new TrivialDiscreteValuationField<>(this), toString());
	}

	@Override
	public boolean isZeroDivisor(T t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isDivisible(T dividend, T divisor) {
		if (!divisor.equals(zero())) {
			return true;
		}
		return dividend.equals(zero());
	}

	@Override
	public PivotStrategy<T> preferredPivotStrategy() {
		return new TrivialPivotStrategy<>(this);
	}

	@Override
	public QuotientAndRemainderResult<T> quotientAndRemainder(T dividend, T divisor) {
		if (dividend.equals(zero()) && divisor.equals(zero())) {
			return new QuotientAndRemainderResult<>(zero(), zero());
		}
		return new QuotientAndRemainderResult<>(this.divide(dividend, divisor), zero());
	}

	@Override
	public T divideChecked(T dividend, T divisor) {
		if (divisor.equals(zero())) {
			throw new ArithmeticException("Division by zero");
		}
		return divide(dividend, divisor);
	}

	@Override
	public T remainder(T dividend, T divisor) {
		return divisor.equals(zero()) ? dividend : zero();
	}

	@Override
	public BigInteger euclidMeasure(T t) {
		return BigInteger.ZERO;
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return this.getNumberOfElements().subtract(BigInteger.ONE);
	}

	@Override
	public Iterable<T> getUnits() {
		return this.getNonZeroElements();
	}

	public T gcd(T t1, T t2) {
		return this.extendedEuclidean(t1, t2).getGcd();
	}

	public T lcm(T t1, T t2) {
		return divide(multiply(t1, t2), gcd(t1, t2));
	}

	@Override
	public ExtendedEuclideanListResult<T> extendedEuclidean(List<T> ts) {
		List<T> result = new ArrayList<>();
		if (ts.size() == 0) {
			return new ExtendedEuclideanListResult<>(zero(), result);
		} else if (ts.size() == 1) {
			result.add(one());
			return new ExtendedEuclideanListResult<>(ts.get(0), result);
		}
		ExtendedEuclideanListResult<T> end = extendedEuclidean(ts.subList(1, ts.size()));
		T gcd = end.getGcd();
		ExtendedEuclideanResult<T> start = extendedEuclidean(ts.get(0), gcd);
		result.add(start.getCoeff1());
		for (T coeff : end.getCoeffs()) {
			result.add(multiply(start.getCoeff2(), coeff));
		}
		return new ExtendedEuclideanListResult<>(start.getGcd(), result);
	}

	public ExtendedEuclideanResult<T> extendedEuclidean(T t1, T t2) {
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(zero(), zero(), zero());
		}
		if (t1.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(one(), zero(), inverse(t2));
		}
		if (t2.equals(this.zero()) || this.characteristic().equals(BigInteger.TWO)) {
			return new ExtendedEuclideanResult<>(one(), inverse(t1), zero());
		}
		return new ExtendedEuclideanResult<>(one(), this.inverse(this.multiply(2, t1)),
				this.inverse(this.multiply(2, t2)));
	}

	@Override
	public BezoutIdentityResult<T> bezoutIdentity(Ideal<T> t1, Ideal<T> t2) {
		if (!t1.contains(one()) && !t2.contains(one())) {
			throw new ArithmeticException("Ideals not coprime!");
		}
		if (!t1.contains(one())) {
			return new BezoutIdentityResult<>(Collections.singletonList(zero()), Collections.singletonList(one()));
		}
		if (!t2.contains(one()) || characteristic().equals(BigInteger.TWO)) {
			return new BezoutIdentityResult<>(Collections.singletonList(one()), Collections.singletonList(zero()));
		}
		T oneHalf = inverse(getInteger(2));
		return new BezoutIdentityResult<>(Collections.singletonList(oneHalf), Collections.singletonList(oneHalf));
	}

	@Override
	public Pair<T, T> bezoutIdentity(T t1, T t2) {
		if (t1.equals(zero()) && t2.equals(zero())) {
			throw new ArithmeticException("Elements not coprime!");
		}
		if (t1.equals(zero())) {
			return new Pair<>(zero(), inverse(t2));
		}
		if (t2.equals(zero()) || characteristic().equals(BigInteger.TWO)) {
			return new Pair<>(inverse(t1), zero());
		}
		return new Pair<>(inverse(multiply(2, t1)), inverse(multiply(2, t2)));
	}

	@Override
	public boolean coprime(Ideal<T> t1, Ideal<T> t2) {
		return t1.contains(one()) || t2.contains(one());
	}

	@Override
	public boolean coprime(T t1, T t2) {
		return !t1.equals(zero()) || !t2.equals(zero());
	}

	@Override
	public ChineseRemainderPreparation<T> prepareChineseRemainderTheorem(List<Ideal<T>> ideals) {
		if (ideals.size() > 1) {
			throw new ArithmeticException("Not coprime proper ideals!");
		}
		if (ideals.size() == 1 && ideals.get(0).contains(one())) {
			throw new ArithmeticException("Not proper ideals!");
		}
		if (ideals.size() == 1) {
			return new ChineseRemainderPreparation<>(ideals, getZeroIdeal(), Collections.singletonList(one()));
		}
		return new ChineseRemainderPreparation<>(ideals, getZeroIdeal(), Collections.emptyList());
	}

	@Override
	public T chineseRemainderTheorem(List<T> elements, ChineseRemainderPreparation<T> preparation) {
		if (elements.size() != preparation.getMultipliers().size()) {
			throw new ArithmeticException("Mismatched multipliers!");
		}
		if (elements.size() == 0) {
			return zero();
		}
		return elements.get(0);
	}

	@Override
	public T chineseRemainderTheorem(List<T> elements, List<Ideal<T>> ideals) {
		if (elements.size() != ideals.size()) {
			throw new ArithmeticException("Mismatched multipliers!");
		}
		if (elements.size() == 0) {
			return zero();
		}
		return elements.get(0);
	}

	@Override
	public UnivariatePolynomialRing<T> getUnivariatePolynomialRing() {
		return ring;
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<T> t) {
		return factorization(t).isIrreducible();
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public final boolean hasCharacteristicRoot(T t) {
		return hasCharacteristicRoot(t, 1);
	}

	@Override
	public boolean hasCharacteristicRoot(T t, int power) {
		if (characteristic().equals(BigInteger.ZERO)) {
			return false;
		}
		throw new UnsupportedOperationException();
	}

	@Override
	public final T characteristicRoot(T t) {
		return characteristicRoot(t, 1);
	}

	@Override
	public T characteristicRoot(T t, int power) {
		throw new UnsupportedOperationException();
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain() {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start) {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start, Ideal<T> end) {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return true;
	}

	@Override
	public FactorizationResult<T, T> uniqueFactorization(T t) {
		return new FactorizationResult<>(t, Collections.emptySortedMap());
	}

	@Override
	public boolean isIrreducible(T t) {
		throw new ArithmeticException("not a non zero non unit! (it's a field!)");
	}

	@Override
	public boolean isPrime(T t) {
		throw new ArithmeticException("not a non zero non unit! (it's a field!)");
	}

	@Override
	public boolean isSquareFree(T t) {
		throw new ArithmeticException("not a non zero non unit! (it's a field!)");
	}

	@Override
	public List<T> factors(T t) {
		return Collections.singletonList(one());
	}

	@Override
	public List<T> adicDevelopment(T t, T base) {
		throw new ArithmeticException("No adic development over fields!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public boolean isDedekindDomain() {
		return true;
	}

	@Override
	public FactorizationResult<Ideal<T>, Ideal<T>> idealFactorization(Ideal<T> t) {
		if (t.equals(getZeroIdeal())) {
			return new FactorizationResult<>(getUnitIdeal(), SingletonSortedMap.map(getZeroIdeal(), 1));
		}
		return new FactorizationResult<>(getUnitIdeal(), new TreeMap<>());
	}

	@Override
	public PrimaryDecompositionResult<T, FieldIdeal<T>> primaryDecomposition(Ideal<T> t) {
		if (t.equals(getZeroIdeal())) {
			return new PrimaryDecompositionResult<>(Collections.singletonList(getZeroIdeal()),
					Collections.singletonList(getZeroIdeal()));
		}
		return new PrimaryDecompositionResult<>(Collections.emptyList(), Collections.emptyList());
	}

	@Override
	public IdealResult<T, FieldIdeal<T>> getIdealWithTransforms(List<T> generators) {
		FieldIdeal<T> ideal = new FieldIdeal<>(generators, this);
		List<T> expression = new ArrayList<>();
		boolean found = false;
		for (T g : generators) {
			if (!found && !g.equals(zero())) {
				found = true;
				expression.add(inverse(g));
			} else {
				expression.add(zero());
			}
		}
		return new IdealResult<>(Collections.singletonList(expression), generators, ideal);
	}

	@SuppressWarnings("unchecked")
	@Override
	public IdealResult<T, FieldIdeal<T>> getIdealWithTransforms(T... generators) {
		return getIdealWithTransforms(Arrays.asList(generators));
	}

	@Override
	public FieldIdeal<T> getIdeal(List<T> generators) {
		return getIdealWithTransforms(generators).getIdeal();
	}

	@SuppressWarnings("unchecked")
	@Override
	public FieldIdeal<T> getIdeal(T... generators) {
		return getIdeal(Arrays.asList(generators));
	}

	@Override
	public FieldIdeal<T> getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

	@Override
	public FieldIdeal<T> getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	public FieldIdeal<T> getNilRadical() {
		return getZeroIdeal();
	}

	public FieldIdeal<T> add(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) && t2.generators().get(0).equals(zero()) ? zero() : one()));
	}

	public FieldIdeal<T> multiply(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
	}

	public FieldIdeal<T> intersect(Ideal<T> t1, Ideal<T> t2) {
		return getIdeal(Collections.singletonList(
				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
	}

	@Override
	public FieldIdeal<T> radical(Ideal<T> t) {
		return (FieldIdeal<T>) t;
	}

	@Override
	public FieldIdeal<T> power(Ideal<T> t, int power) {
		if (power == 0) {
			return getUnitIdeal();
		}
		return (FieldIdeal<T>) t;
	}

	@Override
	public ModuloMaximalIdealResult<T, T> moduloMaximalIdeal(Ideal<T> ideal) {
		if (!ideal.equals(getZeroIdeal())) {
			throw new ArithmeticException("Not a maximal ideal!");
		}
		return new ModuloMaximalIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
	}

	@Override
	public ModuloIdealResult<T, T> moduloIdeal(Ideal<T> ideal) {
		if (ideal.contains(one())) {
			throw new ArithmeticException("Not a proper ideal!");
		}
		return new ModuloIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
	}

	@Override
	public <S extends Element<S>> Ideal<T> getIdealEmbedding(Ideal<S> t, MathMap<S, T> map) {
		for (S generator : t.generators()) {
			if (!map.evaluate(generator).equals(zero())) {
				return getUnitIdeal();
			}
		}
		return getZeroIdeal();
	}

	public static class FieldIdeal<T extends Element<T>> extends AbstractIdeal<T> {
		private boolean isZero;
		private Field<T> field;

		public FieldIdeal(List<T> generators, Field<T> field) {
			super(field);
			this.field = field;
			isZero = true;
			for (T g : generators) {
				if (!g.equals(field.zero())) {
					isZero = false;
					break;
				}
			}
		}

		@Override
		public boolean isPrimary() {
			return isZero;
		}

		@Override
		public boolean isPrime() {
			return isZero;
		}

		@Override
		public boolean isMaximal() {
			return isZero;
		}

		@Override
		public boolean isRadical() {
			return true;
		}

		@Override
		public List<T> generators() {
			return Collections.singletonList(isZero ? field.zero() : field.one());
		}

		@Override
		public List<T> generate(T t) {
			return Collections.singletonList(isZero ? zero() : t);
		}

		public T residue(T t) {
			return isZero ? t : zero();
		}

		@Override
		public boolean isFinite() {
			return isZero || field.isFinite();
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			return isZero ? BigInteger.ONE : field.getNumberOfElements();
		}

		@Override
		public boolean contains(T t) {
			return !isZero || t.equals(zero());
		}

		@Override
		public Value maximumPowerContains(T t) {
			if (!isZero || t.equals(zero())) {
				return Value.INFINITY;
			}
			return Value.ZERO;
		}
	}
}
