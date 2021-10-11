package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.FieldOfFractions.Fraction;
import fields.integers.Integers.IntE;
import fields.interfaces.DedekindRing;
import fields.interfaces.Element;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.GenericUnivariatePolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.pivot.DivisionPivotStrategy;
import fields.vectors.pivot.EuclidPivotStrategy;
import fields.vectors.pivot.PivotStrategy;
import util.Pair;

public abstract class AbstractRing<T extends Element<T>> implements Ring<T> {
	private UnivariatePolynomialRing<T> polynomialRing;

	public AbstractRing() {
		this(true);
	}

	public AbstractRing(boolean generateUnivariatePolynomialRing) {
		if (generateUnivariatePolynomialRing) {
			this.polynomialRing = new GenericUnivariatePolynomialRing<T>(this);
		}
	}

	public Group<T> getAdditiveGroup() {
		return new Group<T>() {

			@Override
			public Exactness exactness() {
				return AbstractRing.this.exactness();
			}

			@Override
			public T getRandomElement() {
				return AbstractRing.this.getRandomElement();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractRing.this.getNumberOfElements();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractRing.this.iterator();
			}

			@Override
			public T neutral() {
				return AbstractRing.this.zero();
			}

			@Override
			public T inverse(T t) {
				return AbstractRing.this.negative(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractRing.this.add(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractRing.this.isFinite();
			}

		};
	}

	public Group<T> getMultiplicativeGroup() {
		return new Group<T>() {

			@Override
			public Exactness exactness() {
				return AbstractRing.this.exactness();
			}

			@Override
			public T getRandomElement() {
				T t;
				do {
					t = AbstractRing.this.getRandomElement();
				} while (!AbstractRing.this.isUnit(t));
				return t;
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractRing.this.getNumberOfUnits();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractRing.this.getUnits().iterator();
			}

			@Override
			public T neutral() {
				return AbstractRing.this.one();
			}

			@Override
			public T inverse(T t) {
				return AbstractRing.this.inverse(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractRing.this.multiply(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractRing.this.isFinite();
			}
		};
	}

	public T add(T t1, T t2, T t3) {
		return this.add(this.add(t1, t2), t3);
	}

	@Override
	public T add(T t1, T t2, T t3, T t4) {
		return this.add(this.add(t1, t2), this.add(t3, t4));
	}

	public T subtract(T minuend, T subtrahend) {
		return this.add(minuend, this.negative(subtrahend));
	}

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
		return this.multiply(this.multiply(t1, t2), t3);
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

	@Override
	public T upToUnit(T t) {
		return divide(t, projectToUnit(t));
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
	public T divideChecked(T dividend, T divisor) {
		QuotientAndRemainderResult<T> result = quotientAndRemainder(dividend, divisor);
		if (!result.getRemainder().equals(zero())) {
			throw new ArithmeticException(dividend + " not divisble by " + divisor + " in " + this + "!");
		}
		return result.getQuotient();
	}

	@Override
	public T divide(T dividend, T divisor) {
		return quotientAndRemainder(dividend, divisor).getQuotient();
	}

	@Override
	public T remainder(T dividend, T divisor) {
		return quotientAndRemainder(dividend, divisor).getRemainder();
	}

	@Override
	public BezoutIdentityResult<T> bezoutIdentity(Ideal<T> t1, Ideal<T> t2) {
		if (isEuclidean()) {
			if (t1.generators().size() != 1 || t2.generators().size() != 1) {
				throw new ArithmeticException("Euclidean rings should use principal ideals");
			}
			ExtendedEuclideanResult<T> ee = extendedEuclidean(t1.generators().get(0), t2.generators().get(0));
			T inverseGcd = inverse(ee.getGcd());
			return new BezoutIdentityResult<>(Collections.singletonList(multiply(inverseGcd, ee.getCoeff1())),
					Collections.singletonList(multiply(inverseGcd, ee.getCoeff2())));
		}
		List<T> generators1 = t1.generators();
		List<T> generators2 = t2.generators();
		List<T> generators = new ArrayList<>();
		generators.addAll(generators1);
		generators.addAll(generators2);
		IdealResult<T, ? extends Ideal<T>> idealResult = getIdealWithTransforms(generators);
		if (!idealResult.getIdeal().contains(one())) {
			throw new ArithmeticException("Ideals not coprime!");
		}
		List<T> generateOne = idealResult.getIdeal().generate(one());
		List<T> bezout1 = new ArrayList<>();
		for (int i = 0; i < generators1.size(); i++) {
			bezout1.add(zero());
		}
		List<T> bezout2 = new ArrayList<>();
		for (int i = 0; i < generators2.size(); i++) {
			bezout2.add(zero());
		}
		for (int i = 0; i < generateOne.size(); i++) {
			T generate = generateOne.get(i);
			List<T> expression = idealResult.getGeneratorExpressions().get(i);
			for (int j = 0; j < generators1.size(); j++) {
				bezout1.set(j, add(bezout1.get(j), multiply(generate, expression.get(j))));
			}
			for (int j = 0; j < generators2.size(); j++) {
				bezout2.set(j, add(bezout2.get(j), multiply(generate, expression.get(j + generators1.size()))));
			}
		}
		return new BezoutIdentityResult<>(bezout1, bezout2);
	}

	@Override
	public Pair<T, T> bezoutIdentity(T t1, T t2) {
		if (isEuclidean()) {
			ExtendedEuclideanResult<T> ee = extendedEuclidean(t1, t2);
			if (!isUnit(ee.getGcd())) {
				throw new ArithmeticException("Not coprime");
			}
			return new Pair<>(divideChecked(ee.getCoeff1(), ee.getGcd()), divideChecked(ee.getCoeff2(), ee.getGcd()));
		}
		IdealResult<T, ?> i1 = getIdealWithTransforms(Collections.singletonList(t1));
		IdealResult<T, ?> i2 = getIdealWithTransforms(Collections.singletonList(t2));
		BezoutIdentityResult<T> bezout = bezoutIdentity(i1.getIdeal(), i2.getIdeal());
		T coeff1 = zero();
		for (int i = 0; i < bezout.getCoeff1().size(); i++) {
			coeff1 = add(coeff1, multiply(i1.getGeneratorExpressions().get(i).get(0), bezout.getCoeff1().get(i)));
		}
		T coeff2 = zero();
		for (int i = 0; i < bezout.getCoeff2().size(); i++) {
			coeff2 = add(coeff2, multiply(i2.getGeneratorExpressions().get(i).get(0), bezout.getCoeff2().get(i)));
		}
		return new Pair<>(coeff1, coeff2);
	}

	@Override
	public boolean coprime(Ideal<T> t1, Ideal<T> t2) {
		return add(t1, t2).contains(one());
	}

	@Override
	public boolean coprime(T t1, T t2) {
		if (isUniqueFactorizationDomain()) {
			return isUnit(gcd(t1, t2));
		}
		return coprime(getIdeal(Collections.singletonList(t1)), getIdeal(Collections.singletonList(t2)));
	}

	@Override
	public ChineseRemainderPreparation<T> prepareChineseRemainderTheorem(List<Ideal<T>> ideals) {
		Ideal<T> product = getUnitIdeal();
		List<T> multipliers = new ArrayList<>();
		for (int i = 0; i < ideals.size(); i++) {
			product = multiply(product, ideals.get(i));
			Ideal<T> cofactor = getUnitIdeal();
			for (int j = 0; j < ideals.size(); j++) {
				if (i == j) {
					continue;
				}
				cofactor = multiply(cofactor, ideals.get(j));
			}
			BezoutIdentityResult<T> bezout = bezoutIdentity(ideals.get(i), cofactor);
			T multiplier = one();
			List<T> generators = ideals.get(i).generators();
			for (int j = 0; j < generators.size(); j++) {
				multiplier = subtract(multiplier, multiply(bezout.getCoeff1().get(j), generators.get(j)));
			}
			multipliers.add(multiplier);
		}
		return new ChineseRemainderPreparation<>(ideals, product, multipliers);
	}

	@Override
	public T chineseRemainderTheorem(List<T> elements, ChineseRemainderPreparation<T> preparation) {
		if (elements.size() != preparation.getMultipliers().size()) {
			throw new ArithmeticException("Unequal list sizes!");
		}
		T result = zero();
		for (int i = 0; i < elements.size(); i++) {
			result = preparation.getProduct()
					.residue(add(result, multiply(elements.get(i), preparation.getMultipliers().get(i))));
		}
		return result;
	}

	@Override
	public T chineseRemainderTheorem(List<T> elements, List<Ideal<T>> ideals) {
		return chineseRemainderTheorem(elements, prepareChineseRemainderTheorem(ideals));
	}

	@Override
	public DedekindRing<T, ?, ?> asDedekindRing() {
		throw new ArithmeticException("Not a Dedekind ring!");
	}

	@Override
	public PivotStrategy<T> preferredPivotStrategy() {
		if (isEuclidean()) {
			return new EuclidPivotStrategy<>(this);
		}
		return new DivisionPivotStrategy<>(this);
	}

	public T gcd(T t1, T t2) {
		if (this.isEuclidean()) {
			return this.extendedEuclidean(t1, t2).getGcd();
		}
		if (this.isUniqueFactorizationDomain()) {
			FactorizationResult<T, T> factorization1 = uniqueFactorization(t1);
			FactorizationResult<T, T> factorization2 = uniqueFactorization(t2);
			T gcd = one();
			Set<T> primes = new TreeSet<>();
			primes.addAll(factorization1.primeFactors());
			primes.addAll(factorization2.primeFactors());
			for (T prime : primes) {
				gcd = multiply(gcd,
						power(prime, factorization1.multiplicity(prime) + factorization2.multiplicity(prime)));
			}
			return gcd;
		}
		throw new RuntimeException("Not possible");
	}

	public T lcm(T t1, T t2) {
		return divideChecked(multiply(t1, t2), gcd(t1, t2));
	}

	@Override
	public FieldOfFractionsResult<T, ?> fieldOfFractions() {
		FieldOfFractions<T> fieldOfFractions = new FieldOfFractions<>(this);
		return new FieldOfFractionsResult<>(this, fieldOfFractions, new MathMap<>() {

			@Override
			public Fraction<T> evaluate(T t) {
				return fieldOfFractions.getEmbedding(t);
			}
		}, new MathMap<>() {

			@Override
			public T evaluate(Fraction<T> t) {
				return t.getNumerator();
			}
		}, new MathMap<>() {

			@Override
			public T evaluate(Fraction<T> t) {
				return t.getDenominator();
			}
		}, new MathMap<>() {

			@Override
			public T evaluate(Fraction<T> t) {
				return t.asInteger();
			}
		});
	}
	
	@Override
	public LocalizeResult<T, ?, ?, ?> localizeAtIdeal(Ideal<T> primeIdeal) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	public ExtendedEuclideanResult<T> extendedEuclidean(T t1, T t2) {
		if (!this.isEuclidean())
			throw new RuntimeException("Not possible");
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(zero(), zero(), zero());
		}
		if (t1.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(t2, zero(), one());
		}
		if (t2.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(t1, one(), zero());
		}
		// T unit1 = projectToUnit(t1);
		// T unit2 = projectToUnit(t2);
		// t1 = upToUnit(t1);
		// t2 = upToUnit(t2);
		T u, v, xp, x, yp, y, q, r;
		u = t1;
		v = t2;
		xp = one();
		yp = zero();
		x = zero();
		y = one();
		// x(0) * t1 + y(0) * t2 = u = r(0)
		// x(1) * t1 + y(1) * t2 = v = r(1)
		while (true) {
			QuotientAndRemainderResult<T> qr = this.quotientAndRemainder(u, v);
			q = qr.getQuotient();
			r = qr.getRemainder();
			if (r.equals(this.zero()))
				break;
			// u = q * v + r
			// r(n+1) = r(n-1) - q * r(n)
			// r(n+1) = (x(n-1) - q * x(n)) * t1 + (y(n-1) - q * y(n)) * t2
			// x(n+1) = x(n-1) - q * x(n)
			// y(n+1) = y(n-1) - q * y(n)
			u = v;
			v = r;
			r = x;
			x = this.subtract(xp, this.multiply(q, x));
			xp = r;
			r = y;
			y = this.subtract(yp, this.multiply(q, y));
			yp = r;
		}
		T unit = inverse(projectToUnit(v));
		return new ExtendedEuclideanResult<>(upToUnit(v), multiply(x, unit), multiply(y, unit));
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

	@Override
	public List<T> adicDevelopment(T t, T base) {
		if (!this.isEuclidean()) {
			throw new ArithmeticException("Adic development only makes sense for Euclidean rings");
		}
		List<T> coefficients = new ArrayList<>();
		while (!t.equals(zero())) {
			QuotientAndRemainderResult<T> qr = quotientAndRemainder(t, base);
			coefficients.add(qr.getRemainder());
			t = qr.getQuotient();
		}
		return coefficients;
	}

	@Override
	public BigInteger getNumberOfUnits() {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public FactorizationResult<Ideal<T>, Ideal<T>> idealFactorization(Ideal<T> t) {
		if (!isUniqueFactorizationDomain() || !isPrincipalIdealDomain()) {
			throw new UnsupportedOperationException("Default implementation only works for UFDs!");
		}
		T generator = t.generators().get(0);
		FactorizationResult<T, T> generatorFactorization = uniqueFactorization(generator);
		SortedMap<Ideal<T>, Integer> result = new TreeMap<>();
		for (T prime : generatorFactorization.primeFactors()) {
			result.put(getIdeal(Collections.singletonList(prime)), generatorFactorization.multiplicity(prime));
		}
		return new FactorizationResult<>(getIdeal(Collections.singletonList(one())), result);
	}

	@Override
	public PrimaryDecompositionResult<T, ? extends Ideal<T>> primaryDecomposition(Ideal<T> t) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public ModuloMaximalIdealResult<T, ?> moduloMaximalIdeal(Ideal<T> ideal) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public ModuloIdealResult<T, ?> moduloIdeal(Ideal<T> ideal) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@SuppressWarnings("unchecked")
	@Override
	public IdealResult<T, ?> getIdealWithTransforms(T... generators) {
		return getIdealWithTransforms(Arrays.asList(generators));
	}

	@Override
	public Ideal<T> getIdeal(List<T> generators) {
		return getIdealWithTransforms(generators).getIdeal();
	}

	@SuppressWarnings("unchecked")
	@Override
	public Ideal<T> getIdeal(T... generators) {
		return getIdeal(Arrays.asList(generators));
	}

	@Override
	public <S extends Element<S>> Ideal<T> getIdealEmbedding(Ideal<S> t, MathMap<S, T> map) {
		List<T> generators = new ArrayList<>();
		for (S generator : t.generators()) {
			generators.add(map.evaluate(generator));
		}
		return getIdeal(generators);
	}

	@Override
	public Ideal<T> getUnitIdeal() {
		return getIdeal(Collections.singletonList(one()));
	}

	@Override
	public Ideal<T> getZeroIdeal() {
		return getIdeal(Collections.emptyList());
	}

	@Override
	public Ideal<T> getNilRadical() {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public List<T> factors(T t) {
		FactorizationResult<T, T> factorization = uniqueFactorization(t);
		List<T> result = new ArrayList<>();
		Map<T, Integer> powers = new TreeMap<>();
		for (T prime : factorization.primeFactors()) {
			powers.put(prime, 0);
		}
		T activePrime;
		do {
			T factor = one();
			activePrime = null;
			List<T> nullPrimes = new ArrayList<>();
			for (T prime : powers.keySet()) {
				factor = multiply(factor, power(prime, powers.get(prime)));
				if (activePrime == null && powers.get(prime) < factorization.multiplicity(prime)) {
					activePrime = prime;
				} else if (activePrime == null) {
					nullPrimes.add(prime);
				}
			}
			if (activePrime != null) {
				for (T prime : nullPrimes) {
					powers.put(prime, 0);
				}
				powers.put(activePrime, powers.get(activePrime) + 1);
			}
			result.add(factor);
		} while (activePrime != null);
		return result;
	}

	public Ideal<T> add(Ideal<T> t1, Ideal<T> t2) {
		List<T> generators = new ArrayList<T>();
		for (T gen1 : t1.generators()) {
			generators.add(gen1);
		}
		for (T gen2 : t2.generators()) {
			generators.add(gen2);
		}
		return getIdeal(generators);
	}

	public Ideal<T> multiply(Ideal<T> t1, Ideal<T> t2) {
		List<T> generators = new ArrayList<T>();
		for (T gen1 : t1.generators()) {
			for (T gen2 : t2.generators()) {
				generators.add(multiply(gen1, gen2));
			}
		}
		return getIdeal(generators);
	}

	public Ideal<T> power(Ideal<T> t, int power) {
		Ideal<T> result = getUnitIdeal();
		Ideal<T> square = t;
		while (power > 0) {
			if ((power & 1) == 1) {
				result = multiply(result, square);
			}
			power >>= 1;
			square = multiply(square, square);
		}
		return result;
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain() {
		return maximalPrimeIdealChain(getNilRadical());
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public List<Ideal<T>> maximalPrimeIdealChain(Ideal<T> start, Ideal<T> end) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	public UnivariatePolynomialRing<T> getUnivariatePolynomialRing() {
		if (this.polynomialRing == null) {
			this.polynomialRing = new GenericUnivariatePolynomialRing<T>(this);
		}
		return this.polynomialRing;
	}

	@Override
	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isIrreducible(T t) {
		if (!isUniqueFactorizationDomain()) {
			throw new UnsupportedOperationException();
		}
		return uniqueFactorization(t).isIrreducible();
	}

	@Override
	public boolean isPrime(T t) {
		if (!isUniqueFactorizationDomain()) {
			throw new UnsupportedOperationException();
		}
		return isIrreducible(t);
	}

	@Override
	public boolean isSquareFree(T t) {
		if (!isUniqueFactorizationDomain()) {
			throw new UnsupportedOperationException();
		}
		return uniqueFactorization(t).squareFree();
	}

	@Override
	public boolean hasRoots(Polynomial<T> t) {
		FactorizationResult<Polynomial<T>, T> factors = factorization(getUnivariatePolynomialRing().toUnivariate(t));
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
		Monomial constant = getUnivariatePolynomialRing().getMonomial(new int[] { 0 });
		Monomial linear = getUnivariatePolynomialRing().getMonomial(new int[] { 1 });
		FactorizationResult<Polynomial<T>, T> factors = factorization(getUnivariatePolynomialRing().toUnivariate(t));
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
		Monomial m = getUnivariatePolynomialRing().getMonomial(new int[] { n });
		Monomial c = getUnivariatePolynomialRing().getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = getUnivariatePolynomialRing().getPolynomial(coeff);
		return this.hasRoots(p);
	}

	@Override
	public Map<T, Integer> roots(T t, int n) {
		if (t.equals(this.zero())) {
			return Collections.singletonMap(this.zero(), n);
		}
		Monomial m = getUnivariatePolynomialRing().getMonomial(new int[] { n });
		Monomial c = getUnivariatePolynomialRing().getMonomial(new int[] { 0 });
		SortedMap<Monomial, T> coeff = new TreeMap<>();
		coeff.put(m, this.one());
		coeff.put(c, this.negative(t));
		Polynomial<T> p = getUnivariatePolynomialRing().getPolynomial(coeff);
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

}
