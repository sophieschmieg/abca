package fields.integers;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DedekindRing;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.LocalRing;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.GenericUnivariatePolynomial;
import fields.polynomials.Monomial;
import util.MiscAlgorithms;
import util.SingletonSortedMap;

public class Integers extends AbstractRing<IntE> implements DedekindRing<IntE, Fraction, PFE> {
	private static Integers z = new Integers();
	private Map<IntE, FactorizationResult<IntE, IntE>> factorizatonCache;
	private Map<Polynomial<IntE>, SquareFreeFactorizationResult> polynomialFactorizationCache;
	private IntE zero;
	private IntE one;
	private List<IntE> smallPrimes;

	private Integers() {
		factorizatonCache = new TreeMap<>();
		polynomialFactorizationCache = new TreeMap<>();
		zero = new IntE(0);
		one = new IntE(1);
		smallPrimes = new ArrayList<>();
	}

	public static Integers z() {
		return z;
	}

	public static class IntE extends AbstractElement<IntE> {
		private BigInteger value;

		public IntE(int value) {
			this.value = BigInteger.valueOf(value);
		}

		public IntE(BigInteger value) {
			this.value = value;
		}

		@Override
		public int compareTo(IntE other) {
			return this.value.compareTo(((IntE) other).value);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}

		public int intValueExact() {
			return value.intValueExact();
		}
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	public Iterator<IntE> primes() {
		return new Iterator<>() {
			private IntE prime = zero();
			private Iterator<IntE> it = smallPrimes.iterator();

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public IntE next() {
				if (it != null && it.hasNext()) {
					try {
						prime = it.next();
					} catch (ConcurrentModificationException e) {
						it = smallPrimes.iterator();
						if (!prime.equals(zero())) {
							while (!it.next().equals(prime))
								;
							prime = it.next();
						}
					}
				} else {
					it = null;
					prime = new IntE(prime.value.nextProbablePrime());
					smallPrimes.add(prime);
				}
				return prime;
			}
		};
	}

	public Iterable<IntE> setOfPrimes() {
		return new Iterable<>() {

			@Override
			public Iterator<IntE> iterator() {
				return primes();
			}
		};
	}

	@Override
	public IntE getRandomElement() {
		return new IntE((int) Math.round(new Random().nextGaussian() * 10.0));
	}

	public IntE getRandomElement(IntE max) {
		return new IntE(MiscAlgorithms.randomBigInteger(new Random(), max.value));
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
	public Iterator<IntE> iterator() {
		return new Iterator<IntE>() {
			private BigInteger it = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public IntE next() {
				IntE next = new IntE(it);
				if (it.compareTo(BigInteger.ZERO) > 0) {
					it = it.negate();
				} else {
					it = it.negate().add(BigInteger.ONE);
				}
				return next;
			}
		};
	}

	@Override
	public IntE zero() {
		return zero;
	}

	@Override
	public IntE one() {
		return one;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public IntE add(IntE t1, IntE t2) {
		return new IntE(t1.value.add(t2.value));
	}

	@Override
	public IntE negative(IntE t) {
		return new IntE(t.value.negate());
	}

	@Override
	public IntE multiply(IntE t1, IntE t2) {
		return new IntE(t1.value.multiply(t2.value));
	}

	@Override
	public IntE multiply(IntE t1, IntE t2, IntE t3) {
		return new IntE(t1.value.multiply(t2.value).multiply(t3.value));
	}

	@Override
	public IntE getInteger(BigInteger t) {
		return new IntE(t);
	}

	@Override
	public IntE divide(IntE t1, IntE t2) {
		return new IntE(t1.value.divide(t2.value));
	}

	@Override
	public IntE remainder(IntE t1, IntE t2) {
		return new IntE(t1.value.mod(t2.value));
	}

	@Override
	public boolean isUnit(IntE t) {
		return t.value.equals(BigInteger.ONE) || t.value.equals(BigInteger.ONE.negate());
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return BigInteger.TWO;
	}

	@Override
	public Iterable<IntE> getUnits() {
		List<IntE> list = new ArrayList<>();
		list.add(new IntE(1));
		list.add(new IntE(-1));
		return list;
	}

	@Override
	public IntE inverse(IntE t) {
		return t;
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
	public boolean isZeroDivisor(IntE t) {
		return t.value.equals(BigInteger.ZERO);
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return true;
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
	public boolean isDivisible(IntE dividend, IntE divisor) {
		return dividend.value.mod(divisor.value.abs()).equals(BigInteger.ZERO);
	}

	@Override
	public BigInteger euclidMeasure(IntE t) {
		return t.value.abs();
	}

	@Override
	public QuotientAndRemainderResult<IntE> quotientAndRemainder(IntE dividend, IntE divisor) {
		BigInteger[] result = dividend.value.divideAndRemainder(divisor.value);
		return new QuotientAndRemainderResult<>(new IntE(result[0]), new IntE(result[1]));
	}

	@Override
	public IntE gcd(IntE t1, IntE t2) {
		return new IntE(t1.value.gcd(t2.value));
	}

	public boolean likelySmooth(IntE n, IntE max) {
		if (n.equals(zero())) {
			return true;
		}
		PartialFactorizationResult pollardRhoResult = pollardRhoFactorization(n, 1,
				max.getValue().sqrt().intValueExact() + 1);
		return !pollardRhoResult.cofactor.equals(n);
	}

//	public List<Vector<IntE>> polynomialSieve(Polynomial<IntE> f, IntE lowerBoundX, IntE upperBoundX, IntE lowerBoundY,
//			IntE upperBoundY, List<IntE> smoothnessBase) {
//		IntE xDiff = subtract(upperBoundX, lowerBoundX);
//		IntE yDiff = subtract(upperBoundY, lowerBoundY);
//		if (xDiff.compareTo(zero()) <= 0 || yDiff.compareTo(zero()) <= 0) {
//			return Collections.emptyList();
//		}
//		if (xDiff.compareTo(getInteger(Integer.MAX_VALUE)) > 0 || yDiff.compareTo(getInteger(Integer.MAX_VALUE)) > 0) {
//			throw new ArithmeticException("Too many candidates!");
//		}
//		Rationals q = Rationals.q();
//		UnivariatePolynomialRing<IntE> polynomials = getUnivariatePolynomialRing();
//		PolynomialRing<IntE> polynomials2 = AbstractPolynomialRing.getPolynomialRing(this, 2, Monomial.GREVLEX);
//		if (!polynomials2.isHomogeneous(f)) {
//			throw new ArithmeticException("Not homogenous");
//		}
//		UnivariatePolynomial<IntE> dehomogenized = polynomials.getEmbedding(polynomials2.dehomogenize(f, 2),
//				new int[] { 0 });
//		List<IntE> zeroList = new ArrayList<>();
//		zeroList.add(null);
//		zeroList.add(zero());
//		UnivariatePolynomial<IntE> zeroed = polynomials.getEmbedding(polynomials2.partiallyEvaluate(f, zeroList),
//				new int[] { 0 });
//		List<Vector<IntE>> candidates = new ArrayList<>();
//		List<List<IntE>> reconstructed = new ArrayList<>();
//		for (IntE i = lowerBoundX; i.compareTo(upperBoundX) < 0; i = add(i, one())) {
//			List<IntE> reconstructedY = new ArrayList<>();
//			reconstructed.add(reconstructedY);
//			for (IntE j = lowerBoundY; j.compareTo(upperBoundY) < 0; j = add(j, one())) {
//				candidates.add(new Vector<>(i, j));
//				reconstructedY.add(one());
//			}
//		}
//		for (IntE prime : smoothnessBase) {
//			PrimeField reduction = PrimeField.getPrimeField(prime);
//			int power = -1;
//			UnivariatePolynomial<PFE> reduced;
//			IntE primePower = one();
//			do {
//				power++;
//				reduced = reduceUnivariatePolynomial(polynomials.divideScalar(f, primePower), prime);
//				primePower = multiply(prime, primePower);
//			} while (reduced.equals(reduction.getUnivariatePolynomialRing().zero()));
//			primePower = divideChecked(primePower, prime);
//			if (power > 0) {
//				for (int i = 0; i < reconstructed.size(); i++) {
//					reconstructed.set(i, multiply(reconstructed.get(i), primePower));
//				}
//			}
//			Set<IntE> roots = new TreeSet<>();
//			for (PFE root : reduction.roots(reduced).keySet()) {
//				roots.add(lift(reduction.subtract(root, reduce(lowerBound, prime.getValue()))));
//			}
//			while (!roots.isEmpty()) {
//				Set<IntE> nextRoots = new TreeSet<>();
//				primePower = multiply(prime, primePower);
//				for (IntE root : roots) {
//					if (root.compareTo(getInteger(reconstructed.size())) >= 0) {
//						continue;
//					}
//					int increment = primePower.compareTo(z.subtract(upperBound, lowerBound)) < 0
//							? primePower.intValueExact()
//							: reconstructed.size();
//					for (int index = root.intValueExact(); index < reconstructed.size(); index += increment) {
//						reconstructed.set(index, multiply(reconstructed.get(index), prime));
//					}
//					IntE lifted = add(root, lowerBound);
//					IntE result = polynomials.evaluate(f, lifted);
//					result = divideChecked(result, prime);
//					if (result.equals(zero()) && primePower.compareTo(subtract(upperBound, lowerBound)) >= 0) {
//						continue;
//					}
//					UnivariatePolynomial<PFE> nextReduced = reduceUnivariatePolynomial(
//							polynomials
//									.divideScalar(
//											polynomials.substitute(f,
//													Collections.singletonList(
//															polynomials.getPolynomial(lifted, primePower))),
//											primePower),
//							prime);
//					if (nextReduced.equals(reduction.getUnivariatePolynomialRing().zero())) {
//						for (int i = 0; i < prime.intValueExact(); i++) {
//							nextRoots.add(add(root, multiply(i, primePower)));
//						}
//					} else {
//						for (PFE nextRoot : reduction.roots(nextReduced).keySet()) {
//							nextRoots.add(add(root, multiply(lift(nextRoot), primePower)));
//						}
//					}
//				}
//				power++;
//				roots = nextRoots;
//			}
//		}
//		List<Vector<IntE>> result = new ArrayList<>();
//		for (Vector<IntE> candidate : candidates) {
//			IntE eval = polynomials.evaluate(f, candidate);
//			if (eval.compareTo(zero()) < 0) {
//				eval = negative(eval);
//			}
//			if (eval.equals(reconstructed.get(subtract(candidate.get(1), lowerBoundX).intValueExact())
//					.get(subtract(candidate.get(2), lowerBoundY).intValueExact()))) {
//				result.add(candidate);
//			}
//		}
//		return result;
//
//	}

	public List<IntE> polynomialSieveSingleVariable(GenericUnivariatePolynomial<IntE> f, IntE lowerBound,
			IntE upperBound, List<IntE> smoothnessBase) {
		if (lowerBound.compareTo(upperBound) >= 0) {
			return Collections.emptyList();
		}
		if (subtract(upperBound, lowerBound).compareTo(getInteger(Integer.MAX_VALUE)) > 0) {
			throw new ArithmeticException("Too many candidates!");
		}
		UnivariatePolynomialRing<IntE> polynomials = getUnivariatePolynomialRing();
		List<IntE> candidates = new ArrayList<>();
		List<IntE> reconstructed = new ArrayList<>();
		for (IntE i = lowerBound; i.compareTo(upperBound) < 0; i = add(i, one())) {
			candidates.add(i);
			reconstructed.add(one());
		}
		for (IntE prime : smoothnessBase) {
			PrimeField reduction = PrimeField.getPrimeField(prime);
			int power = -1;
			UnivariatePolynomial<PFE> reduced;
			IntE primePower = one();
			do {
				power++;
				reduced = reduceUnivariatePolynomial(polynomials.divideScalar(f, primePower), prime);
				primePower = multiply(prime, primePower);
			} while (reduced.equals(reduction.getUnivariatePolynomialRing().zero()));
			primePower = divideChecked(primePower, prime);
			if (power > 0) {
				for (int i = 0; i < reconstructed.size(); i++) {
					reconstructed.set(i, multiply(reconstructed.get(i), primePower));
				}
			}
			Set<IntE> roots = new TreeSet<>();
			for (PFE root : reduction.roots(reduced).keySet()) {
				roots.add(lift(reduction.subtract(root, reduce(lowerBound, prime.getValue()))));
			}
			while (!roots.isEmpty()) {
				Set<IntE> nextRoots = new TreeSet<>();
				primePower = multiply(prime, primePower);
				for (IntE root : roots) {
					if (root.compareTo(getInteger(reconstructed.size())) >= 0) {
						continue;
					}
					int increment = primePower.compareTo(z.subtract(upperBound, lowerBound)) < 0
							? primePower.intValueExact()
							: reconstructed.size();
					for (int index = root.intValueExact(); index < reconstructed.size(); index += increment) {
						reconstructed.set(index, multiply(reconstructed.get(index), prime));
					}
					IntE lifted = add(root, lowerBound);
					IntE result = polynomials.evaluate(f, lifted);
					result = divideChecked(result, prime);
					if (result.equals(zero()) && primePower.compareTo(subtract(upperBound, lowerBound)) >= 0) {
						continue;
					}
					UnivariatePolynomial<PFE> nextReduced = reduceUnivariatePolynomial(
							polynomials
									.divideScalar(
											polynomials.substitute(f,
													Collections.singletonList(
															polynomials.getPolynomial(lifted, primePower))),
											primePower),
							prime);
					if (nextReduced.equals(reduction.getUnivariatePolynomialRing().zero())) {
						for (int i = 0; i < prime.intValueExact(); i++) {
							nextRoots.add(add(root, multiply(i, primePower)));
						}
					} else {
						for (PFE nextRoot : reduction.roots(nextReduced).keySet()) {
							nextRoots.add(add(root, multiply(lift(nextRoot), primePower)));
						}
					}
				}
				power++;
				roots = nextRoots;
			}
		}
		List<IntE> result = new ArrayList<>();
		for (int i = 0; i < candidates.size(); i++) {
			IntE eval = polynomials.evaluate(f, candidates.get(i));
			if (eval.compareTo(zero()) < 0) {
				eval = negative(eval);
			}
			if (eval.equals(reconstructed.get(i))) {
				result.add(candidates.get(i));
			}
		}
		return result;
	}

	public Optional<FactorizationResult<IntE, IntE>> uniqueFactorizationIfSmooth(IntE n, SortedSet<IntE> primes) {
		if (n.equals(zero())) {
			return Optional.empty();
		}
		if (!likelySmooth(n, primes.last())) {
			return Optional.empty();
		}
		SortedMap<IntE, Integer> factorization = new TreeMap<>();
		for (IntE prime : primes) {
			int multiplicity = 0;
			while (isDivisible(n, prime)) {
				n = divideChecked(n, prime);
				multiplicity++;
			}
			factorization.put(prime, multiplicity);
		}
		if (isUnit(n)) {
			return Optional.of(new FactorizationResult<>(n, factorization));
		}
		return Optional.empty();
	}

	@Override
	public boolean isIrreducible(IntE t) {
		return t.value.isProbablePrime(100);
	}

	@Override
	public FactorizationResult<IntE, IntE> uniqueFactorization(IntE n) {
		if (factorizatonCache.containsKey(n)) {
			return factorizatonCache.get(n);
		}
		IntE t = n;
		IntE unit = one();
		if (n.value.compareTo(BigInteger.ZERO) < 0) {
			unit = negative(unit);
			n = negative(n);
		}
		FactorizationResult<IntE, IntE> result = new FactorizationResult<>(unit, uniqueFactorization(n, false));
		factorizatonCache.put(t, result);
		return result;
	}

	private class PartialFactorizationResult {
		private Map<IntE, Integer> factors;
		private IntE cofactor;

		public PartialFactorizationResult(Map<IntE, Integer> factors, IntE cofactor) {
			this.factors = factors;
			this.cofactor = cofactor;
		}
	}

	private SortedMap<IntE, Integer> uniqueFactorization(IntE n, boolean skipNaive) {
		if (n.value.compareTo(BigInteger.ZERO) < 0) {
			n = negative(n);
		}
		if (n.value.equals(BigInteger.ZERO)) {
			throw new ArithmeticException("0 cannot be factorized!");
		}
		if (n.value.isProbablePrime(100)) {
			return SingletonSortedMap.map(n, 1);
		}
		SortedMap<IntE, Integer> decomposition = new TreeMap<>();
		if (!skipNaive) {
			PartialFactorizationResult naive = naivePrimeDecomposition(n, BigInteger.valueOf(7919));
			n = naive.cofactor;
			decomposition.putAll(naive.factors);
		}
		PartialFactorizationResult pollardRho = pollardRhoFactorization(n, 100, 1000);
		decomposition.putAll(pollardRho.factors);
		n = pollardRho.cofactor;
		if (!n.equals(one())) {
			throw new ArithmeticException("Could not factor number");
		}
		return decomposition;
	}

	private PartialFactorizationResult pollardRhoFactorization(IntE n, int numberOfTries, int numberOfIterations) {
		if (n.equals(one())) {
			return new PartialFactorizationResult(Collections.emptyMap(), n);
		}
		if (n.value.isProbablePrime(100)) {
			return new PartialFactorizationResult(Collections.singletonMap(n, 1), one());
		}
		for (int i = 0; i < numberOfTries; i++) {
			IntE a = getRandomElement(n);
			IntE x = getRandomElement(n);
			IntE y = getRandomElement(n);
			for (int j = 0; j < numberOfIterations; j++) {
				IntE gcd = gcd(n, subtract(x, y));
				if (!gcd.equals(one())) {
					if (gcd.equals(n)) {
						break;
					}
					Map<IntE, Integer> result = new TreeMap<>();
					result.putAll(uniqueFactorization(gcd, true));
					result.putAll(uniqueFactorization(divideChecked(n, gcd), true));
					return new PartialFactorizationResult(result, one());
				}
				x = pollardRhoG(x, n, a);
				y = pollardRhoG(pollardRhoG(y, n, a), n, a);
			}
		}
		return new PartialFactorizationResult(Collections.emptyMap(), n);
	}

	private IntE pollardRhoG(IntE input, IntE n, IntE a) {
		return remainder(add(multiply(input, input), a), n);
	}

	private PartialFactorizationResult naivePrimeDecomposition(IntE n, BigInteger max) {
		if (n.equals(one())) {
			return new PartialFactorizationResult(Collections.emptyMap(), n);
		}
		if (n.value.isProbablePrime(100)) {
			return new PartialFactorizationResult(Collections.singletonMap(n, 1), one());
		}
		Map<IntE, Integer> result = new TreeMap<>();
		IntE zero = zero();
		IntE one = one();

		Iterator<IntE> primes = primes();
		IntE prime = primes.next();
		int power = 0;
		while (prime.value.compareTo(n.value.sqrt().add(one.value)) < 0 && prime.value.compareTo(max) < 0) {
			while (remainder(n, prime).equals(zero)) {
				power++;
				n = divideChecked(n, getInteger(prime));
				result.put(prime, power);
			}
			prime = primes.next();
			power = 0;
		}
		return new PartialFactorizationResult(result, n);
	}

	public IntE eulerToitent(IntE t) {
		FactorizationResult<IntE, IntE> factorization = uniqueFactorization(t);
		IntE result = one();
		for (IntE prime : factorization.primeFactors()) {
			result = multiply(power(prime, factorization.multiplicity(prime) - 1), subtract(prime, one()), result);
		}
		return result;
	}

	@Override
	public IntE projectToUnit(IntE t) {
		int signum = t.value.signum();
		if (signum == 0) {
			return one();
		}
		return getInteger(signum);
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	public Fraction round(Fraction t, IntE modulus) {
		IntE denominator = one();
		IntE invertible = t.getDenominator();
		IntE gcd;
		do {
			gcd = gcd(invertible, modulus);
			invertible = divide(invertible, gcd);
			denominator = multiply(denominator, gcd);
		} while (!gcd.equals(one()));
		BigInteger mod = modulus.getValue().multiply(denominator.getValue());
		BigInteger numerator = t.getNumerator().getValue().mod(mod).multiply(invertible.getValue().modInverse(mod))
				.mod(mod);
		return Rationals.q().getFraction(numerator, denominator.getValue());
	}

	public LocalRing<Fraction, PFE> localize(BigInteger prime) {
		return Rationals.q().withValuation(prime).ringOfIntegers();
	}

	public LocalRing<Fraction, PFE> localize(IntE prime) {
		return localize(prime.value);
	}

	public LocalRing<Fraction, PFE> localize(Ideal<IntE> prime) {
		return localize(prime.generators().get(0));
	}

	public Value valuation(IntE t, BigInteger prime) {
		if (t.equals(zero())) {
			return Value.INFINITY;
		}
		int valuation = 0;
		BigInteger n = t.value;
		while (n.mod(prime).equals(BigInteger.ZERO)) {
			valuation++;
			n = n.divide(prime);
		}
		return new Value(valuation);
	}

	public Value valuation(IntE t, Ideal<IntE> maximalIdeal) {
		return valuation(t, maximalIdeal.generators().get(0).value);
	}

	public IntE archimedeanValue(IntE t) {
		if (t.compareTo(zero) < 0) {
			return negative(t);
		}
		return t;
	}

	@Override
	public boolean hasSqrt(IntE t) {
		if (t.equals(zero())) {
			return true;
		}
		if (t.compareTo(zero()) < 0) {
			return false;
		}
		BigInteger[] sqrt = t.value.sqrtAndRemainder();
		return sqrt[1].equals(BigInteger.ZERO);
	}

	@Override
	public Map<IntE, Integer> sqrt(IntE t) {
		if (t.equals(zero())) {
			return Collections.singletonMap(zero(), 2);
		}
		if (t.compareTo(zero()) < 0) {
			return Collections.emptyMap();
		}
		BigInteger[] sqrt = t.value.sqrtAndRemainder();
		if (!sqrt[1].equals(BigInteger.ZERO)) {
			return Collections.emptyMap();
		}
		Map<IntE, Integer> result = new TreeMap<>();
		result.put(z.getInteger(sqrt[0]), 1);
		result.put(z.negative(z.getInteger(sqrt[0])), 1);
		return result;
	}

	@Override
	public FactorizationResult<Polynomial<IntE>, IntE> factorization(UnivariatePolynomial<IntE> t) {
		UnivariatePolynomialRing<IntE> ring = this.getUnivariatePolynomialRing();
		SortedMap<Polynomial<IntE>, Integer> result = new TreeMap<>();
		IntE content = ring.content(t);
		FactorizationResult<IntE, IntE> contentFactorization = uniqueFactorization(content);
		for (IntE contentFactor : contentFactorization.primeFactors()) {
			result.put(ring.getEmbedding(contentFactor), contentFactorization.multiplicity(contentFactor));
		}
		t = ring.contentFree(t);
		FactorizationResult<Polynomial<IntE>, IntE> squareFreeFactors = ring.squareFreeFactorization(t);
		IntE unit = multiply(contentFactorization.getUnit(), squareFreeFactors.getUnit());
		for (Polynomial<IntE> sff : squareFreeFactors.primeFactors()) {
			int power = squareFreeFactors.multiplicity(sff);
			SquareFreeFactorizationResult factors = factorizeSquareFree(ring.toUnivariate(sff));
			unit = multiply(unit, factors.unit);
			for (Polynomial<IntE> factor : factors.factors) {
				result.put(factor, power);
			}
		}
		return new FactorizationResult<>(unit, result);
	}

	private class CombinedFactors {
		private SortedSet<Integer> usedIndeces;
		private Polynomial<PAdicNumber> combined;

		public CombinedFactors(SortedSet<Integer> usedIndeces, Polynomial<PAdicNumber> combined) {
			this.usedIndeces = usedIndeces;
			this.combined = combined;
		}
	}

	private enum Method {
		LINEAR, QUADRATIC, OKUTSU;
	}

	private class SquareFreeFactorizationResult {
		private List<Polynomial<IntE>> factors;
		private IntE unit;

		private SquareFreeFactorizationResult(List<Polynomial<IntE>> factors, IntE unit) {
			this.factors = factors;
			this.unit = unit;
		}
	}

	private SquareFreeFactorizationResult factorizeSquareFree(UnivariatePolynomial<IntE> t) {
		if (t.degree() == 0) {
			return new SquareFreeFactorizationResult(Collections.singletonList(t), one());
		}
		if (t.degree() < 0) {
			throw new ArithmeticException("Cannot factor 0");
		}
		if (polynomialFactorizationCache.containsKey(t)) {
			return polynomialFactorizationCache.get(t);
		}
		UnivariatePolynomial<IntE> f = t;
		BigInteger limit = BigInteger.ZERO;
		for (Monomial m : t.monomials()) {
			IntE c = t.coefficient(m);
			BigInteger l = c.getValue().abs();
			if (l.compareTo(limit) > 0) {
				limit = l;
			}
		}
		limit = limit.add(BigInteger.ONE);
		limit = limit.shiftLeft(1);
		for (IntE prime : setOfPrimes()) {
			if (BigInteger.valueOf(t.degree()).mod(prime.value).equals(BigInteger.ZERO)) {
				continue;
			}
			int requiredAccuracy = (int) (Math.log(limit.doubleValue()) / Math.log(prime.value.doubleValue())) + 2;
			PAdicField qp = new PAdicField(prime.value, requiredAccuracy);
			LocalRing<PAdicNumber, PFE> zp = qp.ringOfIntegers();
			UnivariatePolynomialRing<PAdicNumber> zpr = zp.getUnivariatePolynomialRing();
			UnivariatePolynomial<PAdicNumber> fp = zpr.getEmbedding(t, new MathMap<>() {
				@Override
				public PAdicNumber evaluate(IntE t) {
					return qp.getInteger(t);
				}
			});
			if (!zp.hasGoodReduction(fp)) {
				continue;
			}
			UnivariatePolynomialRing<PFE> reducedRing = zp.reduction().getUnivariatePolynomialRing();
			UnivariatePolynomial<PFE> reduced = zp.reduceUnivariatePolynomial(fp);
			List<Polynomial<PAdicNumber>> liftedFactors = new ArrayList<>();
			System.err.println("Preparation successful");
			Method method = Method.OKUTSU;
			switch (method) {
			case LINEAR:
				FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> padicFactorization = zp
						.henselLiftFactorization(zpr.normalize(fp), requiredAccuracy);
				for (Polynomial<PAdicNumber> factor : padicFactorization.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(factor);
					}
				}
				break;
			case OKUTSU:
				FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> padicFactorizationOkutsu = zp
						.factorization(zpr.normalize(fp));
				for (Polynomial<PAdicNumber> factor : padicFactorizationOkutsu.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(factor);
					}
				}
				break;
			case QUADRATIC:
				FactorizationResult<Polynomial<PFE>, PFE> factors = zp.reduction().factorization(reduced);
				System.err.println("Finite Field factorization successful");
				for (Polynomial<PFE> factor : factors.primeFactors()) {
					if (factor.degree() > 0) {
						liftedFactors.add(zp.henselLiftFactor(zpr.normalize(fp), reducedRing.normalize(factor)));
					}
				}
				break;
			}
			System.err.println("Lifts computed successful");

			Map<Integer, Polynomial<PAdicNumber>> padicFactors = new TreeMap<>();
			for (int i = 0; i < liftedFactors.size(); i++) {
				padicFactors.put(i, liftedFactors.get(i));
			}
			List<CombinedFactors> padicCombinedFactors = new ArrayList<>();
			padicCombinedFactors.add(new CombinedFactors(Collections.emptySortedSet(), zpr.one()));
			List<Polynomial<IntE>> intFactors = new ArrayList<>();
			while (padicFactors.size() > 0) {
				CheckCombinationsResult result = checkCombinations(t, padicFactors, padicCombinedFactors, qp);
				t = result.cofactor;
				intFactors.addAll(result.factors);
				for (int k : result.usedFactors) {
					padicFactors.remove(k);
				}
				padicCombinedFactors = result.combined;
			}
			if (t.degree() != 0) {
				throw new ArithmeticException("Factorization recombination failed!");
			}
			SquareFreeFactorizationResult result = new SquareFreeFactorizationResult(intFactors,
					t.leadingCoefficient());
			polynomialFactorizationCache.put(f, result);
			return result;
//			System.err.println("Fast check not successful");
//
//			liftedFactors.clear();
//			for (int index : padicFactors.keySet()) {
//				liftedFactors.add(padicFactors.get(index));
//				zpr.trace(padicFactors.get(index), 5);
//			}
//			// for (int i = 0; i< liftedFactors.size();i++) {
//			// List<Polynomial<Fraction<IntegerRingElement>>> factors =
//			// }
//
		}
		throw new ArithmeticException("ran out of primes");
	}

	private static class CheckCombinationsResult {
		private List<Polynomial<IntE>> factors = new ArrayList<>();
		private UnivariatePolynomial<IntE> cofactor;
		private Set<Integer> usedFactors = new TreeSet<>();
		private List<CombinedFactors> combined = new ArrayList<>();
	}

	private CheckCombinationsResult checkCombinations(UnivariatePolynomial<IntE> t,
			Map<Integer, Polynomial<PAdicNumber>> padicFactors, List<CombinedFactors> padicCombinedFactors,
			PAdicField qp) {
		CheckCombinationsResult result = new CheckCombinationsResult();
		UnivariatePolynomialRing<PAdicNumber> ring = qp.getUnivariatePolynomialRing();
		result.cofactor = t;
		for (int i : padicFactors.keySet()) {
			Polynomial<PAdicNumber> padicFactor = padicFactors.get(i);
			for (CombinedFactors padicCombinedFactor : padicCombinedFactors) {
				if (padicCombinedFactor.usedIndeces.size() != 0 && padicCombinedFactor.usedIndeces.first() <= i) {
					continue;
				}
				SortedSet<Integer> indeces = new TreeSet<>();
				indeces.addAll(padicCombinedFactor.usedIndeces);
				indeces.add(i);
				Polynomial<PAdicNumber> newCombined = ring.multiply(padicFactor, padicCombinedFactor.combined);
				CheckFactorResult cfr = checkFactor(result.cofactor, newCombined, qp);
				if (cfr.foundFactor) {
					result.factors.add(cfr.factor);
					result.cofactor = cfr.cofactor;
					result.usedFactors.addAll(indeces);
					break;
				} else {
					result.combined.add(new CombinedFactors(indeces, newCombined));
				}
			}
		}
		return result;
	}

	private static class CheckFactorResult {
		private boolean foundFactor = false;
		private Polynomial<IntE> factor = null;
		private UnivariatePolynomial<IntE> cofactor = null;
	}

	private CheckFactorResult checkFactor(UnivariatePolynomial<IntE> t, Polynomial<PAdicNumber> potentialFactor,
			PAdicField qp) {
		CheckFactorResult result = new CheckFactorResult();
		UnivariatePolynomialRing<IntE> ring = getUnivariatePolynomialRing();
		UnivariatePolynomialRing<PAdicNumber> qpRing = qp.getUnivariatePolynomialRing();
		potentialFactor = qpRing.multiply(
				qp.divide(qp.getInteger(t.leadingCoefficient()), potentialFactor.leadingCoefficient()),
				potentialFactor);
		Polynomial<IntE> factor = ring.getEmbedding(potentialFactor, new MathMap<>() {
			@Override
			public IntE evaluate(PAdicNumber number) {
				return qp.roundToInteger(number, qp.getAccuracy());
			}
		});
		factor = ring.contentFree(factor);
		QuotientAndRemainderResult<Polynomial<IntE>> qr = ring.quotientAndRemainder(t, factor);
		if (qr.getRemainder().equals(ring.zero())) {
			result.foundFactor = true;
			result.factor = factor;
			result.cofactor = ring.toUnivariate(qr.getQuotient());
		}
		return result;
	}

	@Override
	public String toString() {
		return "Z";
	}

	@Override
	public IdealResult<IntE, IntegerIdeal> getIdealWithTransforms(List<IntE> generators) {
		if (generators.size() == 0) {
			return new IdealResult<>(Collections.singletonList(Collections.emptyList()), generators,
					new IntegerIdeal(0));
		}
		ExtendedEuclideanListResult<IntE> extendedEuclidean = extendedEuclidean(generators);
		return new IdealResult<>(Collections.singletonList(extendedEuclidean.getCoeffs()), generators,
				new IntegerIdeal(extendedEuclidean.getGcd()));
	}

	@Override
	public Ideal<IntE> intersect(Ideal<IntE> t1, Ideal<IntE> t2) {
		return getIdeal(Collections.singletonList(lcm(t1.generators().get(0), t2.generators().get(0))));
	}

	@Override
	public Ideal<IntE> radical(Ideal<IntE> t) {
		IntE m = t.generators().get(0);
		FactorizationResult<IntE, IntE> factors = uniqueFactorization(m);
		IntE radical = one();
		for (IntE prime : factors.primeFactors()) {
			radical = multiply(radical, prime);
		}
		return getIdeal(Collections.singletonList(radical));
	}

	public IntE lift(ModularIntegerRingElement t) {
		return new IntE(t.getValue());
	}

	public IntE lift(PFE t) {
		return new IntE(t.getValue());
	}

	public IntE lift(PFE t, Ideal<IntE> maximalIdeal) {
		return lift(t);
	}

	public IntE centeredLift(PFE t, BigInteger prime) {
		IntE lift = lift(t);
		IntE alternative = subtract(lift, new IntE(prime));
		if (lift.value.compareTo(alternative.value.abs()) > 0) {
			return alternative;
		}
		return lift;
	}

	public IntE centeredLift(PFE t, Ideal<IntE> maximalIdeal) {
		return centeredLift(t, maximalIdeal.generators().get(0).value);
	}

	public Polynomial<IntE> liftPolynomial(Polynomial<PFE> t) {
		return AbstractPolynomialRing.getPolynomialRing(this, t.getPolynomialRing().numberOfVariables(),
				t.getPolynomialRing().getComparator()).getEmbedding(t, new MathMap<>() {
					@Override
					public IntE evaluate(PFE t) {
						return lift(t);
					}
				});
	}

	public UnivariatePolynomial<IntE> liftUnivariatePolynomial(Polynomial<PFE> t) {
		return getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public IntE evaluate(PFE t) {
				return lift(t);
			}
		});
	}

	public UnivariatePolynomial<IntE> centeredLiftUnivariatePolynomial(Polynomial<PFE> t, BigInteger prime) {
		return getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public IntE evaluate(PFE t) {
				return centeredLift(t, prime);
			}
		});
	}

	@Override
	public Field<PFE> reduction(Ideal<IntE> maximalIdeal) {
		return PrimeField.getPrimeField(maximalIdeal.generators().get(0).value);
	}

	@Override
	public Rationals fieldOfFractions() {
		return Rationals.q();
	}

	@Override
	public MathMap<IntE, Fraction> embedding() {
		return new MathMap<>() {
			@Override
			public Fraction evaluate(IntE t) {
				return Rationals.q().getInteger(t);
			}
		};
	}

	@Override
	public boolean isInteger(Fraction t) {
		return isUnit(t.getDenominator());
	}

	@Override
	public IntE asInteger(Fraction t) {
		return t.asInteger();
	}

	public PFE reduce(IntE t, BigInteger prime) {
		return PrimeField.getPrimeField(prime).reduce(t);
	}

	public PFE reduce(IntE t, Ideal<IntE> prime) {
		return reduce(t, prime.generators().get(0).value);
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, BigInteger prime) {
		return AbstractPolynomialRing.getPolynomialRing(PrimeField.getPrimeField(prime),
				t.getPolynomialRing().numberOfVariables(), t.getPolynomialRing().getComparator())
				.getEmbedding(t, new MathMap<>() {
					@Override
					public PFE evaluate(IntE t) {
						return reduce(t, prime);
					}
				});
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, IntE prime) {
		return reducePolynomial(t, prime.getValue());
	}

	public Polynomial<PFE> reducePolynomial(Polynomial<IntE> t, Ideal<IntE> prime) {
		return reducePolynomial(t, prime.generators().get(0));
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, BigInteger prime) {
		return PrimeField.getPrimeField(prime).getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public PFE evaluate(IntE t) {
				return reduce(t, prime);
			}
		});
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, IntE prime) {
		return reduceUnivariatePolynomial(t, prime.getValue());
	}

	public UnivariatePolynomial<PFE> reduceUnivariatePolynomial(UnivariatePolynomial<IntE> t, Ideal<IntE> prime) {
		return reduceUnivariatePolynomial(t, prime.generators().get(0));
	}

	public boolean isEisenstein(GenericUnivariatePolynomial<IntE> p, BigInteger prime) {
		if (!valuation(p.leadingCoefficient(), prime).equals(new Value(0))) {
			return false;
		}
		Value one = new Value(1);
		if (!valuation(p.univariateCoefficient(0), prime).equals(one)) {
			return false;
		}
		for (int i = 1; i < p.degree(); i++) {
			if (valuation(p.univariateCoefficient(i), prime).compareTo(one) < 0) {
				return false;
			}
		}
		return true;
	}

	public static class IntegerIdeal extends AbstractIdeal<IntE> {
		private IntE m;

		public IntegerIdeal(IntE m) {
			super(z);
			this.m = m;
		}

		public IntegerIdeal(int m) {
			this(BigInteger.valueOf(m));
		}

		public IntegerIdeal(BigInteger m) {
			this(new IntE(m));
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
		public List<IntE> generators() {
			return Collections.singletonList(m);
		}

		@Override
		public List<IntE> generate(IntE t) {
			return Collections.singletonList(z.getInteger(t.value.divide(m.value)));
		}

		@Override
		public IntE residue(IntE t) {
			return new IntE(t.value.mod(m.value));
		}

		@Override
		public boolean isPrime() {
			return m.equals(z.zero()) || isMaximal();
		}

		@Override
		public boolean isMaximal() {
			return m.value.isProbablePrime(100);
		}

		@Override
		public boolean contains(IntE t) {
			return t.value.mod(m.value).equals(BigInteger.ZERO);
		}

		@Override
		public Value maximumPowerContains(IntE t) {
			if (t.equals(zero()) || m.equals(z.one())) {
				return Value.INFINITY;
			}
			if (m.equals(zero())) {
				return Value.ZERO;
			}
			FactorizationResult<Ideal<IntE>,Ideal<IntE>> factors = z.idealFactorization(this);
			Value value = Value.INFINITY;
			for (Ideal<IntE> factor : factors.primeFactors()) {
				value = value.min(new Value(z.valuation(t, factor).value() / z.valuation(m, factor).value()));
			}
			return value;
		}
	}
}
